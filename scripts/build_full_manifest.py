import argparse
import logging
import sys
from pathlib import Path
import pandas as pd

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

def is_computational(row):
    """Prüft, ob ein Label auf rein rechnerischer (computational) Evidenz basiert."""
    comp_method = str(row.get("Computational Method", "")).lower()
    evidence = str(row.get("Evidence", "")).lower()
    lab_method = str(row.get("Laboratory Typing Method", "")).lower()
    
    comp_keywords = ["computational", "predicted", "in silico", "amrfinder", "resfinder", "card"]
    
    # Wenn "Computational Method" einen Wert enthält (nicht NaN oder leer), ist es verdächtig
    if pd.notna(row.get("Computational Method")) and str(row.get("Computational Method")).strip():
        # Manche können valid sein, wir schließen aber alles aus, was in Computational Method steht
        return True
        
    for k in comp_keywords:
        if k in comp_method or k in evidence or k in lab_method:
            return True
            
    return False

def find_fasta(genome_id, fasta_dir, recursive=False):
    """Sucht nach der FASTA-Datei für eine exakte genome_id."""
    extensions = [".fna", ".fa", ".fasta"]
    
    if recursive:
        # Glob recursively
        for ext in extensions:
            # We match the exact filename stem
            target_name = f"{genome_id}{ext}"
            for path in fasta_dir.rglob(target_name):
                if path.is_file() and path.stat().st_size > 0:
                    return path
    else:
        for ext in extensions:
            target_name = f"{genome_id}{ext}"
            path = fasta_dir / target_name
            if path.exists() and path.is_file() and path.stat().st_size > 0:
                return path
                
    return None

def main(args=None):
    parser = argparse.ArgumentParser(description="Bereitet das Manifest und die Labels für die AMR-Pipeline vor.")
    
    parser.add_argument("--labels", required=True, type=str, help="Pfad zur rohen TSV-Labeldatei")
    parser.add_argument("--fasta-dir", required=True, type=str, help="Pfad zum FASTA-Verzeichnis")
    parser.add_argument("--manifest-out", default="data/manifests/aureus_full_manifest.csv", type=str, help="Ausgabepfad für das Manifest")
    parser.add_argument("--labels-out", default="data/processed/aureus_labels_long.csv.gz", type=str, help="Ausgabepfad für die Labels")
    parser.add_argument("--antibiotics", nargs="+", default=["cefoxitin", "ciprofloxacin", "erythromycin"], help="Liste der Antibiotika")
    parser.add_argument("--recursive", action="store_true", help="Suche rekursiv nach FASTA-Dateien")
    parser.add_argument("--strict", action="store_true", help="Strikte Validierung")
    parser.add_argument("--dry-run", action="store_true", help="Nur simulieren, nichts schreiben")
    
    parsed_args = parser.parse_args(args)
    setup_logging()
    
    labels_path = Path(parsed_args.labels).resolve()
    fasta_dir = Path(parsed_args.fasta_dir).resolve()
    manifest_out = Path(parsed_args.manifest_out).resolve()
    labels_out = Path(parsed_args.labels_out).resolve()
    target_abs = parsed_args.antibiotics
    
    if not labels_path.exists():
        logging.error(f"Labeldatei existiert nicht: {labels_path}")
        return 1
        
    if not fasta_dir.exists() or not fasta_dir.is_dir():
        logging.error(f"FASTA-Ordner existiert nicht: {fasta_dir}")
        return 1
        
    logging.info(f"Lese TSV: {labels_path}")
    df = pd.read_csv(labels_path, sep="\t", dtype=str)
    
    raw_rows = len(df)
    logging.info(f"Zahl Rohzeilen: {raw_rows}")
    
    # Spaltenvalidierung
    required_cols = ["Genome ID", "Antibiotic", "Resistant Phenotype"]
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        logging.error(f"Erwartete TSV-Spalten fehlen: {missing_cols}")
        return 1
        
    # Standardisiere Spalten
    df = df.rename(columns={
        "Genome ID": "genome_id",
        "Antibiotic": "antibiotic",
        "Resistant Phenotype": "resistant_phenotype"
    })
    
    # NA drop
    df = df.dropna(subset=["genome_id", "antibiotic", "resistant_phenotype"])
    
    # Ausschluss von Computational Predictions
    df["is_comp"] = df.apply(is_computational, axis=1)
    df = df[~df["is_comp"]]
    
    # Antibiotika filtern
    df = df[df["antibiotic"].isin(target_abs)]
    
    # Phänotyp Normalisierung (nur R und S)
    df["label"] = df["resistant_phenotype"].str.strip().str.upper()
    df["label"] = df["label"].replace({"RESISTANT": "R", "SUSCEPTIBLE": "S"})
    df = df[df["label"].isin(["R", "S"])]
    
    # R/S-Konflikte auflösen
    valid_labels = []
    conflict_count = 0
    
    for (gid, ab), group in df.groupby(["genome_id", "antibiotic"]):
        labels = group["label"].unique()
        if len(labels) == 1:
            # Konfliktfrei
            valid_labels.append(group.iloc[0].copy())
        else:
            conflict_count += 1
            
    clean_df = pd.DataFrame(valid_labels)
    filtered_rows = len(clean_df)
    
    logging.info(f"Zahl gefilterte Labelzeilen: {filtered_rows}")
    logging.info(f"Zahl ausgeschlossener Konflikte (pro Genom+AB): {conflict_count}")
    
    if filtered_rows == 0:
        logging.error("Keine gültigen Labels nach Filterung übrig!")
        return 1
        
    # R/S Verteilung loggen
    for ab in target_abs:
        ab_df = clean_df[clean_df["antibiotic"] == ab]
        counts = ab_df["label"].value_counts().to_dict()
        logging.info(f"R/S-Verteilung {ab}: {counts}")
        
    unique_gids = clean_df["genome_id"].unique()
    logging.info(f"Zahl eindeutige Genome: {len(unique_gids)}")
    
    # FASTA Matching
    manifest_records = []
    missing_fastas = 0
    matched_fastas = 0
    
    for gid in unique_gids:
        fasta_path = find_fasta(gid, fasta_dir, recursive=parsed_args.recursive)
        
        if fasta_path:
            manifest_records.append({
                "genome_id": str(gid),
                "fasta_path": str(fasta_path),
                "fasta_filename": fasta_path.name
            })
            matched_fastas += 1
        else:
            missing_fastas += 1
            
    manifest_df = pd.DataFrame(manifest_records)
    logging.info(f"Zahl gematchte FASTAs: {matched_fastas}")
    logging.info(f"Zahl fehlende FASTAs: {missing_fastas}")
    
    if matched_fastas == 0:
        logging.error("Keine einzige FASTA-Datei gefunden!")
        return 1
        
    if not parsed_args.dry_run:
        manifest_out.parent.mkdir(parents=True, exist_ok=True)
        manifest_df.to_csv(manifest_out, index=False)
        logging.info(f"Manifest gespeichert unter: {manifest_out}")
        
        # Labels speichern (nur für gematchte Genome)
        labels_out.parent.mkdir(parents=True, exist_ok=True)
        final_labels = clean_df[clean_df["genome_id"].isin(manifest_df["genome_id"])]
        final_labels.to_csv(labels_out, index=False, compression="gzip")
        logging.info(f"Labels gespeichert unter: {labels_out}")
        
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
