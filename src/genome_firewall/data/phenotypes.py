from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd

REQUIRED_COLUMNS = {
    "Taxon ID",
    "Genome ID",
    "Genome Name",
    "Antibiotic",
    "Resistant Phenotype",
    "Evidence",
}
VALID_LABELS = {"Resistant", "Susceptible"}


@dataclass(frozen=True)
class CleanResult:
    labels: pd.DataFrame
    conflicts: pd.DataFrame
    excluded_counts: dict[str, int]


@dataclass(frozen=True)
class PhenotypeAudit:
    summary: pd.DataFrame
    standards: pd.DataFrame
    methods: pd.DataFrame
    rows: pd.DataFrame


def _stable_rank(value: str, seed: int) -> str:
    return hashlib.sha256(f"{seed}:{value}".encode()).hexdigest()


def audit_source(
    source: Path,
    *,
    species: str,
    taxon_id: int,
    evidence: str,
    antibiotics: Iterable[str],
) -> PhenotypeAudit:
    """Profile laboratory AST provenance before reducing observations to binary labels."""
    frame = pd.read_csv(source, dtype=object, keep_default_na=False)
    missing = REQUIRED_COLUMNS.difference(frame.columns)
    if missing:
        raise ValueError(f"Source CSV is missing columns: {sorted(missing)}")

    supported = {drug.casefold() for drug in antibiotics}
    rows = frame.loc[
        frame["Genome Name"].str.startswith(f"{species} ")
        & frame["Taxon ID"].eq(str(taxon_id))
        & frame["Evidence"].eq(evidence)
    ].copy()
    rows["antibiotic"] = rows["Antibiotic"].str.strip().str.casefold()
    rows["label"] = rows["Resistant Phenotype"].str.strip()
    rows = rows.loc[rows["antibiotic"].isin(supported)].copy()

    optional = [
        "Measurement",
        "Measurement Sign",
        "Measurement Value",
        "Measurement Unit",
        "Laboratory Typing Method",
        "Laboratory Typing Method Version",
        "Laboratory Typing Platform",
        "Vendor",
        "Testing Standard",
        "Testing Standard Year",
        "Source",
        "PubMed",
    ]
    for column in optional:
        if column not in rows:
            rows[column] = ""

    rows["has_measurement"] = rows["Measurement"].str.strip().ne("")
    rows["has_standard"] = rows["Testing Standard"].str.strip().ne("")
    rows["is_binary_label"] = rows["label"].isin(VALID_LABELS)
    rows["genome_id"] = rows["Genome ID"]

    summary = (
        rows.groupby("antibiotic", dropna=False)
        .agg(
            observation_rows=("genome_id", "size"),
            genomes=("genome_id", "nunique"),
            binary_rows=("is_binary_label", "sum"),
            measurement_rows=("has_measurement", "sum"),
            rows_with_standard=("has_standard", "sum"),
        )
        .reset_index()
    )
    label_counts = (
        rows.groupby(["antibiotic", "label"])["genome_id"]
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    summary = summary.merge(label_counts, on="antibiotic", how="left")

    standards = (
        rows.assign(
            testing_standard=rows["Testing Standard"].replace("", "missing"),
            testing_standard_year=rows["Testing Standard Year"].replace("", "missing"),
        )
        .groupby(["antibiotic", "testing_standard", "testing_standard_year"])
        .agg(observation_rows=("genome_id", "size"), genomes=("genome_id", "nunique"))
        .reset_index()
    )
    methods = (
        rows.assign(
            laboratory_method=rows["Laboratory Typing Method"].replace("", "missing")
        )
        .groupby(["antibiotic", "laboratory_method"])
        .agg(observation_rows=("genome_id", "size"), genomes=("genome_id", "nunique"))
        .reset_index()
    )
    audit_columns = [
        "Genome ID",
        "Genome Name",
        "antibiotic",
        "label",
        *optional,
        "has_measurement",
        "has_standard",
        "is_binary_label",
    ]
    return PhenotypeAudit(
        summary=summary.sort_values("antibiotic").reset_index(drop=True),
        standards=standards.sort_values(
            ["antibiotic", "observation_rows"], ascending=[True, False]
        ).reset_index(drop=True),
        methods=methods.sort_values(
            ["antibiotic", "observation_rows"], ascending=[True, False]
        ).reset_index(drop=True),
        rows=rows[audit_columns].sort_values(["antibiotic", "Genome ID"]).reset_index(drop=True),
    )


def write_phenotype_audit(destination: Path, audit: PhenotypeAudit) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    audit.summary.to_csv(destination / "summary.csv", index=False)
    audit.standards.to_csv(destination / "testing-standards.csv", index=False)
    audit.methods.to_csv(destination / "laboratory-methods.csv", index=False)
    audit.rows.to_csv(destination / "observations.csv", index=False)


def load_and_clean(
    source: Path,
    *,
    species: str,
    taxon_id: int,
    evidence: str,
    antibiotics: Iterable[str],
) -> CleanResult:
    """Load explicit lab labels and remove ambiguous genome/drug observations."""
    # Object-backed strings avoid pandas/pyarrow string-array crashes during
    # Streamlit's threaded reruns (notably with pandas 3.0 on macOS).
    frame = pd.read_csv(source, dtype=object, keep_default_na=False)
    missing = REQUIRED_COLUMNS.difference(frame.columns)
    if missing:
        raise ValueError(f"Source CSV is missing columns: {sorted(missing)}")

    antibiotics = {drug.casefold() for drug in antibiotics}
    species_mask = frame["Genome Name"].str.startswith(f"{species} ")
    filtered = frame.loc[
        species_mask
        & frame["Taxon ID"].eq(str(taxon_id))
        & frame["Evidence"].eq(evidence)
    ].copy()
    filtered["antibiotic"] = filtered["Antibiotic"].str.strip().str.casefold()
    filtered["label"] = filtered["Resistant Phenotype"].str.strip()

    excluded_counts = {
        "outside_supported_drugs": int((~filtered["antibiotic"].isin(antibiotics)).sum()),
        "non_binary_or_missing_label": int((~filtered["label"].isin(VALID_LABELS)).sum()),
    }
    filtered = filtered.loc[
        filtered["antibiotic"].isin(antibiotics) & filtered["label"].isin(VALID_LABELS)
    ]

    grouped = filtered.groupby(["Genome ID", "antibiotic"])["label"].nunique()
    conflict_index = grouped[grouped > 1].index
    conflict_keys = pd.DataFrame(conflict_index.tolist(), columns=["genome_id", "antibiotic"])

    labels = filtered.rename(
        columns={"Genome ID": "genome_id", "Genome Name": "genome_name"}
    )[["genome_id", "genome_name", "antibiotic", "label"]].drop_duplicates()
    if not conflict_keys.empty:
        labels = labels.merge(
            conflict_keys.assign(_conflict=True),
            how="left",
            on=["genome_id", "antibiotic"],
        )
        labels = labels.loc[labels["_conflict"].isna()].drop(columns="_conflict")
    excluded_counts["conflicting_genome_drug_pairs"] = len(conflict_keys)
    return CleanResult(
        labels=labels.sort_values(["genome_id", "antibiotic"]).reset_index(drop=True),
        conflicts=conflict_keys,
        excluded_counts=excluded_counts,
    )


def summarize_labels(labels: pd.DataFrame) -> pd.DataFrame:
    """Return per-antibiotic label counts and resistance prevalence."""
    counts = (
        labels.groupby(["antibiotic", "label"])["genome_id"]
        .nunique()
        .unstack(fill_value=0)
        .reset_index()
    )
    for label in VALID_LABELS:
        if label not in counts:
            counts[label] = 0
    counts["total"] = counts["Resistant"] + counts["Susceptible"]
    counts["resistant_fraction"] = counts["Resistant"] / counts["total"]
    return counts.sort_values("total", ascending=False).reset_index(drop=True)


def select_genomes(
    labels: pd.DataFrame,
    *,
    max_genomes: int,
    minimum_per_class: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Select a deterministic, label-rich cohort while protecting minority classes."""
    strata = labels.groupby(["antibiotic", "label"], sort=True)
    # Genomes can satisfy several class floors, so the union is usually far smaller
    # than number-of-drugs x number-of-labels x minimum-per-class.
    protected: set[str] = set()
    for (_, _), stratum in strata:
        candidates = sorted(
            stratum["genome_id"].unique(), key=lambda value: _stable_rank(value, seed)
        )
        protected.update(candidates[:minimum_per_class])
    selected = protected

    coverage = labels.groupby("genome_id")["antibiotic"].nunique()
    remaining = [genome_id for genome_id in coverage.index if genome_id not in selected]
    remaining.sort(key=lambda value: (-int(coverage[value]), _stable_rank(value, seed)))
    selected.update(remaining[: max(0, max_genomes - len(selected))])

    if len(selected) > max_genomes:
        # This can happen only if the requested cap cannot satisfy all class floors.
        raise ValueError(
            f"Class-floor selection needs {len(selected)} genomes, above max_genomes={max_genomes}"
        )

    selected_labels = labels.loc[labels["genome_id"].isin(selected)].copy()
    manifest = (
        selected_labels.groupby(["genome_id", "genome_name"])
        .agg(label_count=("antibiotic", "nunique"))
        .reset_index()
    )
    manifest["selection_rank"] = manifest["genome_id"].map(
        lambda value: _stable_rank(value, seed)
    )
    manifest = manifest.sort_values(["label_count", "selection_rank"], ascending=[False, True])
    return manifest.reset_index(drop=True), selected_labels.reset_index(drop=True)


def write_selection(
    destination: Path,
    manifest: pd.DataFrame,
    labels: pd.DataFrame,
    metadata: dict[str, object],
) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(destination / "genomes.csv", index=False)
    labels.to_csv(destination / "phenotypes.csv", index=False)
    (destination / "selection.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
