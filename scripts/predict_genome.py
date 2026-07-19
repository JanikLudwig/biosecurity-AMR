import argparse
import sys
import logging
from pathlib import Path

from genome_firewall.prediction import predict_genome

def main():
    parser = argparse.ArgumentParser(description="Predict AMR probabilities for a single Staphylococcus aureus genome.")
    parser.add_argument("fasta", type=str, help="Path to the input FASTA file")
    parser.add_argument("--model-dir", type=str, default="artifacts/models", help="Directory containing the .joblib models")
    parser.add_argument("--output-dir", type=str, default="artifacts/predictions", help="Directory for output predictions")
    parser.add_argument("--organism", type=str, default="Staphylococcus_aureus", help="Assumed organism for AMRFinderPlus")
    parser.add_argument("--docker-image", type=str, default="ncbi/amr:4.2.7-2026-05-15.1", help="AMRFinderPlus Docker Image")
    parser.add_argument("--backend", type=str, choices=["docker", "native"], default="docker", help="Execution backend for AMRFinderPlus")
    parser.add_argument("--keep-intermediates", action="store_true", help="Keep intermediate AMRFinderPlus TSV output")
    parser.add_argument("--force", action="store_true", help="Overwrite existing output if run ID matches or retry")
    parser.add_argument("--json", action="store_true", help="Ensure JSON output is generated (default behavior anyway)")
    parser.add_argument("--csv", action="store_true", help="Ensure CSV output is generated (default behavior anyway)")
    parser.add_argument("--quiet", action="store_true", help="Suppress non-error output")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose debug output")
    parser.add_argument("--genome-id", type=str, help="Optional genome ID (defaults to fasta filename stem)")

    args = parser.parse_args()

    logger = logging.getLogger("genome_firewall")
    if args.quiet:
        logger.setLevel(logging.ERROR)
    elif args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    if not logger.hasHandlers():
        handler = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    exit_code = predict_genome(
        fasta_path=Path(args.fasta),
        model_dir=Path(args.model_dir),
        output_dir=Path(args.output_dir),
        organism=args.organism,
        docker_image=args.docker_image,
        backend=args.backend,
        keep_intermediates=args.keep_intermediates,
        genome_id=args.genome_id,
        force=args.force
    )

    sys.exit(exit_code)

if __name__ == '__main__':
    main()
