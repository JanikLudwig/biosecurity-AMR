import json
import logging
from pathlib import Path
import sys

from genome_firewall.modeling.baseline import (
    build_training_dataset,
    assign_development_split,
    train_drug_model
)

import argparse

def main():
    parser = argparse.ArgumentParser(description="Train AMR models.")
    parser.add_argument("--features", type=Path, default=Path("data/local/full_run/features.csv.gz"))
    parser.add_argument("--labels", type=Path, default=Path("data/local/full_run/labels.csv.gz"))
    parser.add_argument("--groups", type=str, default=None)
    parser.add_argument("--calibration-method", type=str, choices=["sigmoid", "isotonic", "none"], default="sigmoid")
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/models"))
    parser.add_argument("--report-dir", type=Path, default=Path("reports/model_training/development"))
    parser.add_argument("--antibiotics", nargs="+", default=["cefoxitin", "ciprofloxacin", "erythromycin"])

    args = parser.parse_args()

    # We remove warnings.simplefilter('ignore') to help catch issues unless we really want them hidden.
    logging.basicConfig(level=logging.INFO, stream=sys.stdout)

    if not args.features.exists() or not args.labels.exists():
        logging.error("Local data not found. Please ensure features and labels exist.")
        return

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.report_dir.mkdir(parents=True, exist_ok=True)

    summaries = []

    for antibiotic in args.antibiotics:
        logging.info(f"Processing {antibiotic}...")
        table, feature_columns = build_training_dataset(args.features, args.labels, antibiotic)

        table, split_metadata = assign_development_split(
            table,
            seed=args.random_state,
            groups_path=args.groups
        )

        metadata, predictions = train_drug_model(
            table,
            feature_columns,
            antibiotic=antibiotic,
            output_directory=args.output_dir,
            seed=args.random_state,
            calibration_method=args.calibration_method,
            split_metadata=split_metadata
        )

        report_path = args.report_dir / f"{antibiotic}_report.json"
        with open(report_path, "w") as f:
            json.dump(metadata, f, indent=2)

        summaries.append(metadata)

        metrics = metadata.get("test_metrics", {})
        ba = metrics.get("balanced_accuracy")
        ba_str = f"{ba:.3f}" if ba is not None else "N/A"
        logging.info(f"Done training {antibiotic}. Test Balanced Accuracy: {ba_str}")

    summary_path = args.report_dir / "summary.json"
    with open(summary_path, "w") as f:
        json.dump(summaries, f, indent=2)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("CRASHED:", repr(e))
        import traceback
        traceback.print_exc()
        sys.exit(1)
