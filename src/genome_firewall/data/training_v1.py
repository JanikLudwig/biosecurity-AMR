from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd


LABEL_MAP = {"R": "Resistant", "S": "Susceptible"}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def prepare_training_v1(
    source_directory: Path,
    *,
    features_output: Path,
    phenotypes_output: Path,
    genomes_output: Path,
    manifest_output: Path,
) -> dict[str, Any]:
    """Normalize the supplied 3k tables into the current pipeline contracts."""
    feature_source = source_directory / "features.csv.gz"
    label_source = source_directory / "aureus_labels_long.csv.gz"
    dictionary_source = source_directory / "feature_dictionary.csv"
    for path in (feature_source, label_source, dictionary_source):
        if not path.is_file():
            raise FileNotFoundError(f"Required training-v1 input is missing: {path}")

    # BV-BRC accessions contain a dot and must never be parsed as floats.
    features = pd.read_csv(
        feature_source, dtype={"genome_id": "string"}, keep_default_na=False
    )
    labels = pd.read_csv(label_source, dtype="string", keep_default_na=False)
    dictionary = pd.read_csv(dictionary_source, dtype="string", keep_default_na=False)

    if features["genome_id"].duplicated().any():
        raise ValueError("features.csv.gz contains duplicate genome_id values")
    if not features["genome_id"].str.fullmatch(r"[0-9]+\.[0-9]+").all():
        raise ValueError("Unexpected BV-BRC genome_id format in features.csv.gz")

    feature_columns = [column for column in features.columns if column != "genome_id"]
    if set(feature_columns) != set(dictionary["feature_id"]):
        raise ValueError("Feature matrix and feature dictionary schemas differ")
    values = features[feature_columns].apply(pd.to_numeric, errors="raise")
    if not values.isin([0, 1]).all().all():
        raise ValueError("The training-v1 feature matrix must be binary")
    features[feature_columns] = values.astype("uint8")

    required = {"genome_id", "Genome Name", "antibiotic", "label", "Evidence"}
    missing = required.difference(labels.columns)
    if missing:
        raise ValueError(f"Training labels are missing columns: {sorted(missing)}")
    labels = labels.loc[labels["label"].isin(LABEL_MAP)].copy()
    labels["label"] = labels["label"].map(LABEL_MAP)
    if labels.duplicated(["genome_id", "antibiotic"]).any():
        raise ValueError("Training labels contain duplicate genome-antibiotic pairs")

    feature_ids = set(features["genome_id"])
    labels = labels.loc[labels["genome_id"].isin(feature_ids)].copy()
    provenance_columns = [
        column
        for column in ["Evidence", "Testing Standard", "Source"]
        if column in labels
    ]
    phenotypes = labels[["genome_id", "antibiotic", "label", *provenance_columns]].rename(
        columns={"Evidence": "evidence"}
    )
    genomes = (
        labels[["genome_id", "Genome Name"]]
        .drop_duplicates("genome_id")
        .rename(columns={"Genome Name": "genome_name"})
    )
    missing_names = sorted(feature_ids.difference(genomes["genome_id"]))
    if missing_names:
        genomes = pd.concat(
            [
                genomes,
                pd.DataFrame({"genome_id": missing_names, "genome_name": ""}),
            ],
            ignore_index=True,
        )

    features_output.parent.mkdir(parents=True, exist_ok=True)
    phenotypes_output.parent.mkdir(parents=True, exist_ok=True)
    genomes_output.parent.mkdir(parents=True, exist_ok=True)
    features.to_parquet(features_output, index=False)
    phenotypes.sort_values(["genome_id", "antibiotic"]).to_csv(
        phenotypes_output, index=False
    )
    genomes.sort_values("genome_id").to_csv(genomes_output, index=False)

    summary: dict[str, Any] = {
        "schema_version": "training-v1.0",
        "genomes": int(features["genome_id"].nunique()),
        "labels": int(len(phenotypes)),
        "antibiotics": sorted(phenotypes["antibiotic"].unique().tolist()),
        "features": len(feature_columns),
        "feature_schema_sha256": hashlib.sha256(
            "\n".join(sorted(feature_columns)).encode("utf-8")
        ).hexdigest(),
        "sources": {
            path.name: {"sha256": _sha256(path), "bytes": path.stat().st_size}
            for path in (feature_source, label_source, dictionary_source)
        },
        "outputs": {
            "features": str(features_output),
            "phenotypes": str(phenotypes_output),
            "genomes": str(genomes_output),
        },
    }
    manifest_output.parent.mkdir(parents=True, exist_ok=True)
    manifest_output.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary
