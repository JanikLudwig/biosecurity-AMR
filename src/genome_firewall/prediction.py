import csv
import json
import logging
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import pandas as pd

from genome_firewall.run_amrfinder import run_single_genome
from genome_firewall.build_feature_matrix import determine_feature
from genome_firewall.modeling.baseline import AMRModel

logger = logging.getLogger(__name__)

def validate_input_fasta(fasta_path: Path) -> None:
    """Validates the input FASTA file."""
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    if not fasta_path.is_file():
        raise ValueError(f"Not a regular file: {fasta_path}")
    if fasta_path.suffix.lower() not in ['.fna', '.fa', '.fasta']:
        raise ValueError(f"Invalid FASTA extension: {fasta_path.suffix}")
    if fasta_path.stat().st_size == 0:
        raise ValueError(f"FASTA file is empty: {fasta_path}")

    has_header = False
    has_sequence = False

    with open(fasta_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                has_header = True
            else:
                has_sequence = True
                # Check for valid IUPAC nucleotide characters + N
                if not re.match(r'^[ACGTURYSWKMBDHVNacgturyswkmbdhvn\-]+$', line):
                    raise ValueError("Sequence contains invalid nucleotide characters.")

            if has_header and has_sequence:
                break

    if not has_header:
        raise ValueError("No FASTA header found (must start with '>').")
    if not has_sequence:
        raise ValueError("No sequence data found in FASTA.")

def parse_amrfinder_output(tsv_path: Path, genome_id: str) -> Tuple[pd.DataFrame, int, int]:
    """
    Parses AMRFinderPlus TSV output into a single-row feature DataFrame.
    Returns: (DataFrame, recognized_count, unknown_count)
    (We don't actually know recognized vs unknown here until we align with models,
     so we just return all features found.)
    Wait, signature will just return (DataFrame, list of all features found)
    """
    if not tsv_path.exists():
        raise FileNotFoundError(f"AMRFinderPlus TSV not found: {tsv_path}")

    features_found = set()

    with open(tsv_path, 'r', encoding='utf-8') as f:
        content = f.read().strip()
        if content:
            lines = content.split('\n')
            if len(lines) > 1:
                reader = csv.DictReader(lines, delimiter='\t')
                for row in reader:
                    feat_id, f_type, gene_sym, mut = determine_feature(row)
                    features_found.add(feat_id)

    # Create single row DataFrame
    data = {'genome_id': genome_id}
    for feat in features_found:
        data[feat] = 1

    df = pd.DataFrame([data])
    return df, list(features_found)

def load_models(model_dir: Path) -> Dict[str, AMRModel]:
    models = {}
    expected_antibiotics = ['cefoxitin', 'ciprofloxacin', 'erythromycin']
    for ab in expected_antibiotics:
        model_path = model_dir / f"{ab}.joblib"
        if not model_path.exists():
            raise FileNotFoundError(f"Model file not found: {model_path}")
        try:
            model = AMRModel(model_path)
            if model.antibiotic != ab:
                raise ValueError(f"Model {ab} has incorrect antibiotic metadata: {model.antibiotic}")
            if len(model.feature_columns) != 158:
                logger.warning(f"Model {ab} has {len(model.feature_columns)} features instead of expected 158.")
            models[ab] = model
        except Exception as e:
            raise RuntimeError(f"Failed to load model {ab}: {e}")
    return models

def align_features(df: pd.DataFrame, models: Dict[str, AMRModel]) -> Tuple[pd.DataFrame, int, int]:
    # AMRModel.predict already aligns features. We just need to report recognized vs unknown.
    # Take the feature schema from one of the models (they should all be the same)
    ref_model = list(models.values())[0]
    schema_features = set(ref_model.feature_columns)

    found_features = set(df.columns) - {'genome_id'}

    recognized = found_features.intersection(schema_features)
    unknown = found_features - schema_features

    return df, len(recognized), len(unknown)

def predict_antibiotics(df: pd.DataFrame, models: Dict[str, AMRModel]) -> List[Dict[str, Any]]:
    predictions = []
    for ab, model in models.items():
        res = model.predict(df)
        row = res.iloc[0]

        # Determine string class based on standard threshold
        pred_class = "R" if row['probability_resistant'] >= 0.5 else "S"

        # Extract features
        evidence = row['evidence_features']
        evidence_list = evidence.split(',') if evidence else []

        pred_dict = {
            "antibiotic": ab,
            "probability_resistant": float(row['probability_resistant']),
            "prediction": pred_class,
            "confidence": float(row['confidence']),
            "evidence_features": evidence_list,
            "calibration_method": model.calibration_method,
            "evaluation_status": model.evaluation_status
        }
        predictions.append(pred_dict)
    return predictions

def write_prediction_report(
    prediction_json_path: Path,
    prediction_csv_path: Path,
    genome_id: str,
    fasta_path: Path,
    organism: str,
    image: str,
    backend: str,
    recognized_count: int,
    unknown_count: int,
    predictions: List[Dict[str, Any]],
    times: Dict[str, float]
):
    no_amr = (recognized_count == 0 and unknown_count == 0)

    # JSON
    report = {
        "genome_id": genome_id,
        "input_fasta": fasta_path.name,
        "assumed_organism": organism,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "pipeline": {
            "amrfinder_image": image,
            "backend": backend
        },
        "feature_summary": {
            "recognized_features": recognized_count,
            "unknown_features": unknown_count,
            "no_amr_features_detected": no_amr
        },
        "predictions": predictions,
        "timing_seconds": times,
        "warnings": [
            "Models are trained only for Staphylococcus aureus.",
            "Development-only evaluation.",
            "Not validated for clinical treatment decisions."
        ]
    }

    with open(prediction_json_path, 'w', encoding='utf-8') as f:
        json.dump(report, f, indent=4)

    # CSV
    with open(prediction_csv_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            "genome_id", "antibiotic", "probability_resistant", "prediction",
            "confidence", "calibration_method", "evaluation_status", "evidence_features"
        ])
        for p in predictions:
            writer.writerow([
                genome_id,
                p["antibiotic"],
                f"{p['probability_resistant']:.4f}",
                p["prediction"],
                f"{p['confidence']:.4f}",
                p["calibration_method"],
                p["evaluation_status"],
                ";".join(p["evidence_features"])
            ])

def predict_genome(
    fasta_path: Path,
    model_dir: Path,
    output_dir: Path,
    organism: str,
    docker_image: str,
    backend: str = "docker",
    keep_intermediates: bool = False,
    genome_id: Optional[str] = None,
    force: bool = False
) -> int:
    import time

    t0 = time.time()
    times = {}

    if not genome_id:
        genome_id = fasta_path.stem

    # Validation
    try:
        validate_input_fasta(fasta_path)
    except Exception as e:
        logger.error(f"Input validation failed: {e}")
        return 2

    t_val = time.time()
    times["validation_seconds"] = round(t_val - t0, 2)

    # Run setup
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_dir = output_dir / genome_id / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    amrfinder_dir = run_dir / "amrfinder"
    amrfinder_dir.mkdir(parents=True, exist_ok=True)

    tsv_path = amrfinder_dir / f"{genome_id}.tsv"

    # Run AMRFinderPlus
    try:
        res = run_single_genome(
            raw_id=genome_id,
            fasta_path=fasta_path,
            output_dir=amrfinder_dir,
            log_dir=run_dir,
            backend=backend,
            image=docker_image,
            organism=organism,
            threads=1,
            plus=True,
            force=force
        )
        if res['status'] == 'failed':
            logger.error(f"AMRFinderPlus failed: {res['error_message']}")
            return 3
    except Exception as e:
        logger.error(f"AMRFinderPlus execution failed: {e}")
        return 3

    t_amr = time.time()
    times["amrfinder_seconds"] = round(t_amr - t_val, 2)

    # Parsing
    try:
        from typing import Tuple
        df, found_features = parse_amrfinder_output(tsv_path, genome_id)
    except Exception as e:
        logger.error(f"Feature parsing failed: {e}")
        return 4

    t_parse = time.time()
    times["parsing_seconds"] = round(t_parse - t_amr, 2)

    # Load Models
    try:
        models = load_models(model_dir)
    except Exception as e:
        logger.error(f"Model loading failed: {e}")
        return 5

    # Align and predict
    try:
        df, rec_count, unk_count = align_features(df, models)
        predictions = predict_antibiotics(df, models)
    except Exception as e:
        logger.error(f"Inference failed: {e}")
        return 6

    t_inf = time.time()
    times["inference_seconds"] = round(t_inf - t_parse, 2)
    times["total_seconds"] = round(t_inf - t0, 2)

    # Write reports
    try:
        write_prediction_report(
            prediction_json_path=run_dir / "prediction.json",
            prediction_csv_path=run_dir / "prediction.csv",
            genome_id=genome_id,
            fasta_path=fasta_path,
            organism=organism,
            image=docker_image,
            backend=backend,
            recognized_count=rec_count,
            unknown_count=unk_count,
            predictions=predictions,
            times=times
        )
    except Exception as e:
        logger.error(f"Writing report failed: {e}")
        return 6

    # Cleanup
    if not keep_intermediates:
        try:
            import shutil
            shutil.rmtree(amrfinder_dir)
        except Exception as e:
            logger.warning(f"Could not remove intermediate files: {e}")

    # Console output
    print(f"\nGenome: {genome_id}")
    print(f"Assumed organism: {organism}")
    print(f"AMRFinderPlus: completed")
    print(f"Recognized AMR features: {rec_count}\n")
    print(f"{'Antibiotic':<16} {'P(resistant)':<14} {'Prediction':<12} {'Confidence'}")
    for p in predictions:
        print(f"{p['antibiotic'].capitalize():<16} {p['probability_resistant']:<14.2f} {p['prediction']:<12} {p['confidence']:.2f}")

    print("\nResearch-use warning:")
    print("Development-only models; not validated for clinical treatment decisions.\n")

    if logger.getEffectiveLevel() == logging.DEBUG:
        logger.debug(f"Run directory: {run_dir}")
        logger.debug(f"Docker Image: {docker_image}")
        logger.debug(f"Unknown Features: {unk_count}")
        logger.debug(f"Times: {times}")

    return 0
