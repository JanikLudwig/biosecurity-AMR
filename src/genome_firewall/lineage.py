from __future__ import annotations

from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd

from genome_firewall.splitting.homology import estimated_ani, sketch_fasta


def build_lineage_reference(
    *,
    splits_path: Path,
    fasta_directory: Path,
    output_path: Path,
    ksize: int,
    scaled: int,
    calibration_quantile: float,
) -> dict[str, Any]:
    splits = pd.read_csv(splits_path, dtype=object, keep_default_na=False)
    sketches = {}
    for row in splits.itertuples(index=False):
        fasta = fasta_directory / f"{row.genome_id}.fna"
        if not fasta.is_file():
            raise FileNotFoundError(f"Lineage reference FASTA is missing: {fasta}")
        sketches[row.genome_id] = sketch_fasta(fasta, ksize=ksize, scaled=scaled)
    train_ids = splits.loc[splits["split"].eq("train"), "genome_id"].tolist()
    calibration_ids = splits.loc[splits["split"].eq("calibration"), "genome_id"].tolist()
    train_sketches = [sketches[genome_id] for genome_id in train_ids]
    if not train_sketches or not calibration_ids:
        raise ValueError("Lineage reference requires train and calibration genomes")

    calibration_rows = []
    for genome_id in calibration_ids:
        maximum = max(
            estimated_ani(sketches[genome_id], reference, ksize=ksize)
            for reference in train_sketches
        )
        calibration_rows.append({"genome_id": genome_id, "maximum_training_ani": maximum})
    calibration = pd.DataFrame(calibration_rows)
    minimum_ani = float(
        np.quantile(calibration["maximum_training_ani"], calibration_quantile)
    )
    artifact = {
        "schema_version": "1.0",
        "ksize": ksize,
        "scaled": scaled,
        "calibration_quantile": calibration_quantile,
        "minimum_training_ani": minimum_ani,
        "training_genome_ids": train_ids,
        "training_sketches": train_sketches,
        "calibration_distances": calibration_rows,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(artifact, output_path, compress=3)
    calibration.to_csv(output_path.with_suffix(".calibration.csv"), index=False)
    return artifact


def evaluate_lineage(fasta: Path, artifact_path: Path) -> dict[str, Any]:
    if not artifact_path.is_file():
        return {
            "status": "unavailable",
            "maximum_training_ani": None,
            "minimum_training_ani": None,
            "nearest_training_genome": None,
        }
    artifact = joblib.load(artifact_path)
    query = sketch_fasta(
        fasta, ksize=int(artifact["ksize"]), scaled=int(artifact["scaled"])
    )
    similarities = [
        estimated_ani(query, reference, ksize=int(artifact["ksize"]))
        for reference in artifact["training_sketches"]
    ]
    best_index = int(np.argmax(similarities))
    maximum = float(similarities[best_index])
    minimum = float(artifact["minimum_training_ani"])
    return {
        "status": "in_distribution" if maximum >= minimum else "out_of_distribution",
        "maximum_training_ani": maximum,
        "minimum_training_ani": minimum,
        "nearest_training_genome": artifact["training_genome_ids"][best_index],
    }
