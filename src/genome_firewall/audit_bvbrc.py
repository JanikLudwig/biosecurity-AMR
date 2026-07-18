import pandas as pd
import json
import logging
from pathlib import Path
import sys

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

AMR_FILE = Path("data/raw/PATRIC_genome_AMR.txt")
METADATA_FILE = Path("data/raw/genome_metadata")
TARGET_TAXON_ID = 1280
REPORTS_DIR = Path("reports")

def main(args=None):
    import argparse
    parser = argparse.ArgumentParser(description="Audit BV-BRC Metadata for AMR")
    parser.parse_args(args)
    
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    
    if not AMR_FILE.exists() or not METADATA_FILE.exists():
        logger.error(f"Fehlende Rohdaten. Bitte 'make audit' oder den Download-Schritt ausführen.")
        return 1

    logger.info("Lade Genome Metadata...")
    # genome_metadata ist typischerweise tab-separated
    meta_df = pd.read_csv(METADATA_FILE, sep="\t", dtype=str)
    
    logger.info("Lade AMR Daten...")
    amr_df = pd.read_csv(AMR_FILE, sep="\t", dtype=str)

    schema_info = {
        "metadata_columns": list(meta_df.columns),
        "amr_columns": list(amr_df.columns)
    }

    # Filter S. aureus in Metadaten
    # Erlaube taxon_id 1280 oder lineage Felder (wenn lineage String vorhanden, suche nach '1280')
    if "taxon_id" in meta_df.columns:
        # Finde alle Genome, die zur Art 1280 gehören (auch Subspezies haben oft andere Taxon-IDs, 
        # aber in BV-BRC ist species_taxon_id oft verfügbar, oder wir filtern über species)
        if "species_taxon_id" in meta_df.columns:
            s_aureus_meta = meta_df[meta_df["species_taxon_id"] == str(TARGET_TAXON_ID)]
        else:
            # Fallback auf String matching in species
            s_aureus_meta = meta_df[meta_df["species"].str.contains("Staphylococcus aureus", na=False)]
    else:
        s_aureus_meta = meta_df[meta_df["species"].str.contains("Staphylococcus aureus", na=False)]

    logger.info(f"Gefundene S. aureus Genome in Metadaten: {len(s_aureus_meta)}")

    # Merge
    merged_df = pd.merge(amr_df, s_aureus_meta, on="genome_id", how="inner")
    logger.info(f"Gefundene AMR Einträge für S. aureus: {len(merged_df)}")

    # Filtere "Computational Prediction"
    if "evidence" in merged_df.columns:
        computational_mask = merged_df["evidence"].str.contains("Computational", case=False, na=False)
        merged_df = merged_df[~computational_mask]
    
    if "lab_typing_method" in merged_df.columns:
        computational_mask2 = merged_df["lab_typing_method"].str.contains("Computational", case=False, na=False)
        merged_df = merged_df[~computational_mask2]

    # Normalisiere Phenotypes
    # BV-BRC Spalte für Resistenz ist meist 'resistant_phenotype'
    pheno_col = "resistant_phenotype"
    if pheno_col not in merged_df.columns:
        logger.error(f"Spalte {pheno_col} nicht gefunden! Verfügbar: {list(merged_df.columns)}")
        return 1

    merged_df[pheno_col] = merged_df[pheno_col].str.strip().str.title()
    allowed_phenotypes = ["Resistant", "Susceptible", "Intermediate"]
    
    # Statistiken pro Antibiotikum
    antibiotic_col = "antibiotic"
    stats = []

    for ab, group in merged_df.groupby(antibiotic_col):
        # Betrachte eindeutige genome_id - phenotype Kombinationen
        unique_calls = group[["genome_id", pheno_col]].drop_duplicates()
        
        counts = unique_calls[pheno_col].value_counts()
        r = counts.get("Resistant", 0)
        s = counts.get("Susceptible", 0)
        i = counts.get("Intermediate", 0)
        
        # Finde Widersprüche (ein Genom hat R und S gleichzeitig für dasselbe AB)
        genomes_per_ab = unique_calls["genome_id"].value_counts()
        conflicts = (genomes_per_ab > 1).sum()

        total = r + s + i
        if total == 0:
            continue

        stats.append({
            "Antibiotikum": ab,
            "Genome (Eindeutig)": len(genomes_per_ab),
            "Resistant": r,
            "Susceptible": s,
            "Intermediate": i,
            "Widersprüche": conflicts,
            "Anteil Resistant (%)": round((r / total) * 100, 2) if total > 0 else 0,
            "Anteil Susceptible (%)": round((s / total) * 100, 2) if total > 0 else 0
        })

    stats_df = pd.DataFrame(stats).sort_values(by="Genome (Eindeutig)", ascending=False)
    
    stats_df.to_csv(REPORTS_DIR / "amr_availability.csv", index=False)
    
    with open(REPORTS_DIR / "amr_availability.md", "w") as f:
        f.write("# AMR Availability Report (Staphylococcus aureus)\n\n")
        f.write(stats_df.to_markdown(index=False))
        f.write("\n\n## Empfohlene Antibiotika für Prototyp\n")
        f.write("Wir empfehlen Antibiotika mit ausreichend Fällen (mindestens >50 R und S) und wenig Konflikten.\n")
        
        candidates = stats_df[(stats_df["Resistant"] >= 50) & (stats_df["Susceptible"] >= 50)].head(5)
        for _, row in candidates.iterrows():
            f.write(f"- **{row['Antibiotikum']}**: {row['Resistant']} R, {row['Susceptible']} S, {row['Widersprüche']} Konflikte\n")

    with open(REPORTS_DIR / "schema_report.json", "w") as f:
        json.dump(schema_info, f, indent=4)

    logger.info("Audit abgeschlossen. Reports in 'reports/' gespeichert.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
