import argparse
import pandas as pd
import json
import logging
from pathlib import Path
import sys
import hashlib
import random

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

def compute_sha256(file_path: Path) -> str:
    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        for block in iter(lambda: f.read(65536), b""):
            sha256.update(block)
    return sha256.hexdigest()

def validate_fasta(file_path: Path) -> dict:
    if not file_path.exists():
        return {"valid": False, "reason": "not_found"}
    size = file_path.stat().st_size
    if size == 0:
        return {"valid": False, "reason": "empty"}
    
    with open(file_path, "rb") as f:
        content = f.read()
    
    if not content:
        return {"valid": False, "reason": "empty"}
        
    content = content.lstrip()
    if not content.startswith(b">"):
        return {"valid": False, "reason": "no_header"}
        
    header_end = content.find(b"\n")
    if header_end == -1:
        return {"valid": False, "reason": "no_sequence"}
        
    # count sequence characters by subtracting newline occurrences
    seq_len = len(content) - header_end - content.count(b"\n", header_end) - content.count(b"\r", header_end)
    
    if seq_len < 1000:
        return {"valid": False, "reason": "too_short"}

    return {"valid": True, "size": size, "seq_len": seq_len}

def run_import(args):
    tsv_path = Path(args.tsv)
    fasta_dir = Path(args.fasta_dir)
    reports_dir = Path(args.reports_dir)
    manifest_out = Path(args.manifest)
    labels_out = Path(args.labels)
    
    reports_dir.mkdir(parents=True, exist_ok=True)
    manifest_out.parent.mkdir(parents=True, exist_ok=True)
    labels_out.parent.mkdir(parents=True, exist_ok=True)

    if not tsv_path.exists():
        logger.error(f"TSV Datei {tsv_path} existiert nicht.")
        sys.exit(1)

    logger.info("Lese TSV...")
    df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    
    # 2. TSV Validieren
    logger.info(f"Anzahl Zeilen TSV: {len(df)}")
    
    required_cols = ["Genome ID", "Genome Name", "Antibiotic", "Resistant Phenotype", "Evidence"]
    for c in required_cols:
        if c not in df.columns:
            logger.error(f"Spalte '{c}' fehlt in TSV.")
            sys.exit(1)
            
    df = df.dropna(subset=["Genome ID", "Antibiotic", "Resistant Phenotype"])
    logger.info(f"Nach Drop NA: {len(df)}")
    df = df.drop_duplicates()
    logger.info(f"Nach Drop Duplicates: {len(df)}")
    
    # Filter Evidence
    computational_keywords = ["Computational Prediction", "Predicted", "Model Prediction", "Machine Learning", "In Silico", "computational"]
    def is_computational(row):
        ev = str(row.get("Evidence", "")).lower()
        cm = str(row.get("Computational Method", "")).lower()
        lm = str(row.get("Laboratory Typing Method", "")).lower()
        combined = ev + " " + cm + " " + lm
        for kw in computational_keywords:
            if kw.lower() in combined:
                return True
        return False
        
    mask_comp = df.apply(is_computational, axis=1)
    df = df[~mask_comp]
    logger.info(f"Nach Filterung Computational: {len(df)}")
    
    # Normalize phenotype
    def norm_pheno(p):
        p = p.strip().lower()
        if p == "resistant": return "R"
        if p == "susceptible": return "S"
        if p == "intermediate": return "I"
        if p == "nonsusceptible": return "NS"
        return "UNKNOWN"
        
    df["norm_label"] = df["Resistant Phenotype"].apply(norm_pheno)
    df = df[df["norm_label"].isin(["R", "S"])]
    logger.info(f"Nach Filterung R/S: {len(df)}")
    
    # Resolve conflicts
    clean_labels = []
    for (gid, ab), group in df.groupby(["Genome ID", "Antibiotic"]):
        labels = group["norm_label"].unique()
        if len(labels) == 1:
            row = group.iloc[0].copy()
            clean_labels.append(row)
            
    clean_df = pd.DataFrame(clean_labels)
    logger.info(f"Eindeutige konfliktfreie R/S Labels: {len(clean_df)}")
    
    # 3. FASTA Matching
    logger.info("Suche FASTA Dateien...")
    fasta_files = []
    for ext in [".fna", ".fa", ".fasta"]:
        fasta_files.extend(list(fasta_dir.rglob(f"*{ext}")))
        
    fasta_map = {}
    duplicate_fastas = set()
    for f in fasta_files:
        name = f.name
        for ext in [".fna", ".fa", ".fasta"]:
            if name.endswith(ext):
                gid = name[: -len(ext)]
                if gid in fasta_map:
                    duplicate_fastas.add(gid)
                else:
                    fasta_map[gid] = f
                break
                
    logger.info(f"Gefundene FASTAs: {len(fasta_files)}")
    
    target_antibiotics = ["cefoxitin", "ciprofloxacin", "erythromycin"]
    pilot_df = clean_df[clean_df["Antibiotic"].isin(target_antibiotics)].copy()
    
    valid_genomes = {}
    for gid in pilot_df["Genome ID"].unique():
        if gid in fasta_map and gid not in duplicate_fastas:
            fpath = fasta_map[gid]
            val = validate_fasta(fpath)
            if val["valid"]:
                valid_genomes[gid] = {
                    "path": fpath,
                    "size": val["size"],
                    "seq_len": val["seq_len"]
                }
                
    logger.info(f"Gültige und gematchte Genome für die 3 Antibiotika: {len(valid_genomes)}")
    
    if len(valid_genomes) < 100:
        logger.error("Weniger als 100 gültige Genome vorhanden!")
        sys.exit(1)
        
    # 5. Auswahl 100 Genome
    candidate_df = pilot_df[pilot_df["Genome ID"].isin(valid_genomes.keys())]
    
    # Stratified/Balanced selection (heuristic)
    # We want max coverage of R/S and multiple antibiotics
    genome_stats = []
    for gid, group in candidate_df.groupby("Genome ID"):
        labels = group.set_index("Antibiotic")["norm_label"].to_dict()
        r_count = list(labels.values()).count("R")
        s_count = list(labels.values()).count("S")
        genome_stats.append({
            "Genome ID": gid,
            "ab_count": len(labels),
            "labels": labels,
            "r_count": r_count,
            "s_count": s_count
        })
        
    genome_stats.sort(key=lambda x: (x["ab_count"], min(x["r_count"], x["s_count"])), reverse=True)
    
    # Just take top 100 to ensure we get those with most labels, breaking ties deterministically
    random.seed(args.seed)
    random.shuffle(genome_stats) # Shuffle first to avoid name bias
    genome_stats.sort(key=lambda x: (x["ab_count"], min(x["r_count"], x["s_count"])), reverse=True)
    
    selected_gids = [x["Genome ID"] for x in genome_stats[:100]]
    selected_gids.sort() # "abschließend nach genome_id sortieren"
    
    if len(set(selected_gids)) != 100:
        logger.error("Doppelte Genome IDs in der Auswahl!")
        sys.exit(1)
        
    # 6. Manifest und Labels
    manifest_rows = []
    for gid in selected_gids:
        info = valid_genomes[gid]
        ab_count = next(x["ab_count"] for x in genome_stats if x["Genome ID"] == gid)
        manifest_rows.append({
            "genome_id": gid,
            "fasta_path": str(info["path"].absolute()),
            "fasta_filename": info["path"].name,
            "fasta_size_bytes": info["size"],
            "fasta_sha256": compute_sha256(info["path"]),
            "sequence_length": info["seq_len"],
            "number_of_binary_labels": ab_count,
            "usable_antibiotics": ";".join(sorted(next(x["labels"].keys() for x in genome_stats if x["Genome ID"] == gid)))
        })
        
    man_df = pd.DataFrame(manifest_rows)
    man_df.to_csv(manifest_out, index=False)
    
    final_labels_df = candidate_df[candidate_df["Genome ID"].isin(selected_gids)]
    out_labels = []
    for _, row in final_labels_df.iterrows():
        out_labels.append({
            "genome_id": row["Genome ID"],
            "antibiotic": row["Antibiotic"],
            "label": row["norm_label"],
            "raw_phenotype": row["Resistant Phenotype"],
            "evidence": row.get("Evidence", "")
        })
        
    pd.DataFrame(out_labels).to_csv(labels_out, index=False, compression="gzip")
    
    # 7. Audit Reports
    audit = {
        "total_label_genomes": len(df["Genome ID"].unique()),
        "total_fasta_files": len(fasta_files),
        "matched_genomes": len(valid_genomes),
        "missing_fastas": len(set(pilot_df["Genome ID"]) - set(fasta_map.keys())),
        "invalid_fastas": sum(1 for gid in fasta_map if gid in pilot_df["Genome ID"].unique() and gid not in valid_genomes and gid not in duplicate_fastas),
        "duplicate_fastas": len(duplicate_fastas),
        "pilot_size": 100
    }
    with open(reports_dir / "input_audit.json", "w") as f:
        json.dump(audit, f, indent=4)
        
    logger.info("Import und Audit erfolgreich abgeschlossen.")
    return 0

def main(args=None):
    parser = argparse.ArgumentParser(description="Import local dataset and prepare 100 genome pilot.")
    parser.add_argument("--tsv", required=True)
    parser.add_argument("--fasta-dir", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--labels", required=True)
    parser.add_argument("--reports-dir", required=True)
    parser.add_argument("--seed", type=int, default=42)
    parsed = parser.parse_args(args)
    return run_import(parsed)

if __name__ == "__main__":
    raise SystemExit(main())
