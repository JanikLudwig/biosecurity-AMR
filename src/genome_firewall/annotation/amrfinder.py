from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd

REQUIRED_OUTPUT_COLUMNS = {
    "Element symbol",
    "Element name",
    "Type",
    "Subtype",
    "Class",
    "Subclass",
    "Method",
    "% Coverage of reference",
    "% Identity to reference",
}


def resolve_executable(configured: str) -> Path | None:
    candidate = Path(configured)
    if candidate.is_file():
        return candidate.resolve()
    discovered = shutil.which("amrfinder")
    return Path(discovered) if discovered else None


def executable_version(executable: Path) -> str:
    environment = {**os.environ, "CONDA_PREFIX": str(executable.parent.parent)}
    result = subprocess.run(
        [str(executable), "--version"],
        capture_output=True,
        text=True,
        check=True,
        env=environment,
    )
    return (result.stdout or result.stderr).strip()


def database_version(executable: Path) -> str | None:
    version_file = executable.parent.parent / "share/amrfinderplus/data/latest/version.txt"
    if not version_file.is_file():
        return None
    return version_file.read_text(encoding="utf-8").strip().splitlines()[0]


def run_nucleotide(
    executable: Path,
    fasta: Path,
    output_tsv: Path,
    *,
    organism: str,
    threads: int = 1,
) -> Path:
    """Run AMRFinderPlus defensively on an existing assembled nucleotide FASTA."""
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    command = [
        str(executable),
        "--nucleotide",
        str(fasta.resolve()),
        "--organism",
        organism,
        "--threads",
        str(threads),
        "--output",
        str(output_tsv.resolve()),
    ]
    environment = {**os.environ, "CONDA_PREFIX": str(executable.parent.parent)}
    result = subprocess.run(
        command, capture_output=True, text=True, check=True, env=environment
    )
    provenance = {
        "command": command,
        "executable_version": executable_version(executable),
        "database_version": database_version(executable),
        "stdout": result.stdout,
        "stderr": result.stderr,
    }
    output_tsv.with_suffix(".provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    return output_tsv


def parse_output(path: Path, *, genome_id: str) -> pd.DataFrame:
    """Parse AMRFinderPlus TSV into explicit known-gene/known-mutation evidence."""
    raw = pd.read_csv(path, sep="\t", dtype=object, keep_default_na=False)
    missing = REQUIRED_OUTPUT_COLUMNS.difference(raw.columns)
    if missing:
        raise ValueError(f"AMRFinder output is missing columns: {sorted(missing)}")

    point_mask = raw["Subtype"].eq("POINT") | raw["Method"].str.startswith("POINT")
    evidence = pd.DataFrame(
        {
            "genome_id": genome_id,
            "feature_key": raw["Element symbol"].map(lambda value: f"gene::{value}"),
            "element_symbol": raw["Element symbol"],
            "element_name": raw["Element name"],
            "evidence_category": "known_resistance_gene",
            "amr_class": raw["Class"],
            "amr_subclass": raw["Subclass"],
            "method": raw["Method"],
            "coverage": pd.to_numeric(raw["% Coverage of reference"], errors="coerce"),
            "identity": pd.to_numeric(raw["% Identity to reference"], errors="coerce"),
        }
    )
    def mutation_key(value: object) -> str:
        symbol = str(value)
        gene, separator, change = symbol.partition("_")
        return f"mutation::{gene}::{change}" if separator else f"mutation::{symbol}::unknown"

    evidence.loc[point_mask, "feature_key"] = raw.loc[point_mask, "Element symbol"].map(
        mutation_key
    )
    evidence.loc[point_mask, "evidence_category"] = "known_resistance_mutation"
    evidence["feature_value"] = 1
    return evidence.sort_values(["evidence_category", "feature_key"]).reset_index(drop=True)
