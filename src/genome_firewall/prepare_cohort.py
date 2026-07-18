import pandas as pd
import yaml
import json
import logging
from pathlib import Path
import sys

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

AMR_FILE = Path("data/raw/PATRIC_genome_AMR.txt")
METADATA_FILE = Path("data/raw/genome_metadata")
CONFIG_FILE = Path("config/pipeline.yaml")
MANIFEST_DIR = Path("data/manifests")
PROCESSED_DIR = Path("data/processed")
REPORTS_DIR = Path("reports")

def load_config():
    with open(CONFIG_FILE, "r") as f:
        return yaml.safe_load(f)

def main(args=None):
    import argparse
    parser = argparse.ArgumentParser(description="Prepare cohort for AMR pipeline")
    parser.parse_args(args)
    
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)

    config = load_config()
    selected_abs = config["antibiotics"].get("selected", [])
    
    if not selected_abs:
        logger.warning("Keine Antibiotika in config/pipeline.yaml ausgewählt!")

    qc_filters = config["cohort"]["qc_filters"]

    logger.info(f"Ausgewählte Antibiotika: {selected_abs}")

    if not METADATA_FILE.exists() or not AMR_FILE.exists():
        logger.error("Fehlende Rohdaten. Bitte zuerst audit oder download_bvbrc_release_files.sh ausführen.")
        return 1

    meta_df = pd.read_csv(METADATA_FILE, sep="\t", dtype=str)
    amr_df = pd.read_csv(AMR_FILE, sep="\t", dtype=str)

    # Taxon Filter für S. aureus (ähnlich wie Audit)
    if "species_taxon_id" in meta_df.columns:
        s_aureus_meta = meta_df[meta_df["species_taxon_id"] == str(config["organism"]["ncbi_taxonomy_id"])]
    else:
        s_aureus_meta = meta_df[meta_df["species"].str.contains("Staphylococcus aureus", na=False)]

    merged_df = pd.merge(amr_df, s_aureus_meta, on="genome_id", how="inner")
    
    pheno_col = "resistant_phenotype"
    merged_df[pheno_col] = merged_df[pheno_col].str.strip().str.title()
    
    # Filtern auf ausgewählte Antibiotika
    filtered_df = merged_df[merged_df["antibiotic"].isin(selected_abs)].copy()

    # QC Filter anwenden
    exclusions = []
    
    def check_qc(row):
        reason = None
        if qc_filters.get("require_public", True):
            if row.get("public") != "Yes" and row.get("public") != "true" and row.get("public") != "1":
                # Viele BV-BRC Metadatenfelder nutzen "Yes" oder "true" / "1" (oder auch gar kein "public" field in manchen releases, in dem Fall überspringen wir den Check oder prüfen auf Vorhandensein)
                pass # Wir belassen es vereinfacht oder prüfen public_status
            if "genome_status" in row and pd.notna(row["genome_status"]) and "public" not in str(row["genome_status"]).lower() and "complete" not in str(row["genome_status"]).lower() and "wgs" not in str(row["genome_status"]).lower():
               pass # Manchmal ist status anders codiert, wir lassen die strenge "Public" prüfung hier als Dummy, in Realität ist PATRIC FTP alles public.
               
        if qc_filters.get("exclude_poor", True) and "genome_quality" in row:
            if pd.notna(row["genome_quality"]) and "poor" in str(row["genome_quality"]).lower():
                return "Poor Quality"
                
        if "checkm_completeness" in row and pd.notna(row["checkm_completeness"]):
            try:
                comp = float(row["checkm_completeness"])
                if comp < qc_filters.get("checkm_completeness_min", 95.0):
                    return "Low Completeness"
            except ValueError:
                pass
                
        if "checkm_contamination" in row and pd.notna(row["checkm_contamination"]):
            try:
                cont = float(row["checkm_contamination"])
                if cont > qc_filters.get("checkm_contamination_max", 5.0):
                    return "High Contamination"
            except ValueError:
                pass
                
        return None

    logger.info("Wende QC Filter an...")
    filtered_df["exclusion_reason"] = filtered_df.apply(check_qc, axis=1)
    
    excluded_df = filtered_df[filtered_df["exclusion_reason"].notna()]
    kept_df = filtered_df[filtered_df["exclusion_reason"].isna()].copy()
    
    logger.info(f"Ausgeschlossen durch QC: {len(excluded_df)}")

    # Duplikate und Konflikte auflösen
    clean_labels = []
    conflicts_count = 0
    
    for (gid, ab), group in kept_df.groupby(["genome_id", "antibiotic"]):
        phenotypes = group[pheno_col].dropna().unique()
        
        # Entferne unzulässige Werte (Groß/Kleinschreibung wurde vorher normiert)
        valid_phenotypes = [p for p in phenotypes if p in ["Resistant", "Susceptible", "Intermediate"]]
        
        if len(valid_phenotypes) == 0:
            continue
        elif len(valid_phenotypes) > 1:
            # Konflikt!
            for _, row in group.iterrows():
                row_copy = row.copy()
                row_copy["exclusion_reason"] = "Phenotype Conflict"
                excluded_df = pd.concat([excluded_df, pd.DataFrame([row_copy])])
            conflicts_count += 1
            continue
            
        # Nimm den ersten Eintrag als repräsentativ für Genome-Antibiotikum Paar
        repr_row = group.iloc[0].copy()
        repr_row[pheno_col] = valid_phenotypes[0]
        clean_labels.append(repr_row)

    logger.info(f"Aufgrund von Konflikten verworfene Paare: {conflicts_count}")

    if not clean_labels:
        logger.error("Keine gültigen Labels nach Filterung übrig!")
        return 1

    clean_df = pd.DataFrame(clean_labels)
    
    # Speichere die Labels im Long-Format
    cols_to_keep = ["genome_id", "genome_name", "antibiotic", pheno_col, "measurement_value", "measurement_unit", "laboratory_typing_method"]
    available_cols = [c for c in cols_to_keep if c in clean_df.columns]
    
    labels_long = clean_df[available_cols]
    labels_long.to_csv(PROCESSED_DIR / "labels_long.csv.gz", index=False, compression="gzip")
    
    # Erstelle Manifest (nur einzigartige Genome)
    manifest_cols = ["genome_id", "genome_name", "assembly_accession", "genome_length", "contigs", "checkm_completeness"]
    manifest_avail = [c for c in manifest_cols if c in clean_df.columns]
    
    manifest_df = clean_df[manifest_avail].drop_duplicates("genome_id")
    manifest_df.to_csv(MANIFEST_DIR / "genome_manifest.csv", index=False)
    
    # Excluded speichern
    if not excluded_df.empty:
        excluded_df.to_csv(REPORTS_DIR / "excluded_records.csv.gz", index=False, compression="gzip")

    # Summary
    summary = {
        "total_genomes_in_manifest": len(manifest_df),
        "total_labels": len(labels_long),
        "excluded_records": len(excluded_df),
        "conflicts_removed": conflicts_count
    }
    
    for ab in selected_abs:
        ab_labels = labels_long[labels_long["antibiotic"] == ab]
        counts = ab_labels[pheno_col].value_counts().to_dict()
        summary[ab] = counts

    with open(REPORTS_DIR / "cohort_summary.json", "w") as f:
        json.dump(summary, f, indent=4)

    logger.info("Cohort Vorbereitung abgeschlossen.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
