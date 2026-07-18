import argparse
import pandas as pd
import yaml
import logging
from pathlib import Path
import sys
from .io_utils import robust_ftps_download, compute_sha256

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

CONFIG_FILE = Path("config/pipeline.yaml")
MANIFEST_FILE = Path("data/manifests/genome_manifest.csv")
RAW_DIR = Path("data/raw")

def load_config():
    with open(CONFIG_FILE, "r") as f:
        return yaml.safe_load(f)

def validate_fasta(file_path: Path) -> bool:
    """Prüft grob, ob die Datei eine gültige FASTA ist."""
    try:
        with open(file_path, "r") as f:
            first_line = f.readline().strip()
            if not first_line.startswith(">"):
                logger.error(f"Ungültige FASTA (kein > in der ersten Zeile): {file_path}")
                return False
            # Prüfe ein paar weitere Zeilen auf sinnvolle Zeichen
            seq_line = f.readline().strip()
            if not seq_line:
                logger.error(f"Leere Sequenz in FASTA: {file_path}")
                return False
        return True
    except Exception as e:
        logger.error(f"Fehler beim Validieren der FASTA {file_path}: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Download FASTA files from BV-BRC.")
    parser.add_argument("--mode", choices=["smoke", "pilot", "full"], required=True, help="Download mode")
    args = parser.parse_args()

    RAW_DIR.mkdir(parents=True, exist_ok=True)
    config = load_config()

    if args.mode == "smoke":
        # NCTC 8325
        ref = config["reference"]
        logger.info(f"Smoke Mode: Versuche Referenzgenom über Assembly {ref['assembly_accession']} zu finden...")
        # Lade temporär Metadaten um BV-BRC genome_id zu finden
        meta_file = RAW_DIR / "genome_metadata"
        if not meta_file.exists():
            logger.error("genome_metadata fehlt für Smoke-Test Auflösung.")
            sys.exit(1)
            
        meta_df = pd.read_csv(meta_file, sep="\t", dtype=str)
        ref_df = meta_df[meta_df["assembly_accession"] == ref["assembly_accession"]]
        if ref_df.empty:
            logger.error(f"Referenzgenom {ref['assembly_accession']} nicht in genome_metadata gefunden!")
            sys.exit(1)
            
        genome_id = ref_df.iloc[0]["genome_id"]
        target_genomes = pd.DataFrame([{"genome_id": genome_id}])
    else:
        if not MANIFEST_FILE.exists():
            logger.error("genome_manifest.csv fehlt. Bitte 'make audit' ausführen.")
            sys.exit(1)
        manifest_df = pd.read_csv(MANIFEST_FILE, dtype=str)
        
        if args.mode == "pilot":
            pilot_size = config["cohort"]["pilot_size"]
            seed = config["cohort"]["random_seed"]
            # Einfaches Sampling (könnte in Zukunft stratifiziert werden nach R/S)
            if len(manifest_df) > pilot_size:
                target_genomes = manifest_df.sample(n=pilot_size, random_state=seed)
            else:
                target_genomes = manifest_df
        else: # full
            target_genomes = manifest_df

    logger.info(f"Geplante Downloads: {len(target_genomes)} Genome.")

    hashes = {}
    for idx, row in target_genomes.iterrows():
        genome_id = row["genome_id"]
        url = f"ftps://ftp.bvbrc.org/genomes/{genome_id}/{genome_id}.fna"
        out_path = RAW_DIR / f"{genome_id}.fna"

        success = robust_ftps_download(url, out_path)
        if success:
            if validate_fasta(out_path):
                file_hash = compute_sha256(out_path)
                hashes[genome_id] = file_hash
            else:
                out_path.unlink()
        else:
            logger.error(f"Download fehlgeschlagen für {genome_id}")

    # Aktualisiere Manifest mit Hash-Werten wenn nicht im Smoke-Modus
    if args.mode != "smoke" and hashes:
        manifest_df["sha256"] = manifest_df["genome_id"].map(hashes)
        manifest_df.to_csv(MANIFEST_FILE, index=False)
        logger.info(f"Hashes im Manifest aktualisiert.")

if __name__ == "__main__":
    main()
