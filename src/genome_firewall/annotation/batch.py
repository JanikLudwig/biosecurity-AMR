from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd

from genome_firewall.annotation.amrfinder import (
    database_version,
    executable_version,
    parse_output,
    run_nucleotide,
)


def _provenance_path(output_tsv: Path) -> Path:
    return output_tsv.with_suffix(".provenance.json")


def _is_current(output_tsv: Path, *, software: str, database: str) -> bool:
    provenance_path = _provenance_path(output_tsv)
    if not output_tsv.is_file() or not provenance_path.is_file():
        return False
    try:
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return (
        provenance.get("executable_version") == software
        and provenance.get("database_version") == database
    )


def annotate_batch(
    fasta_paths: list[Path],
    output_directory: Path,
    *,
    executable: Path,
    organism: str,
    workers: int = 2,
    threads_per_worker: int = 1,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Annotate a resumable FASTA batch and return status, evidence, and features."""
    output_directory.mkdir(parents=True, exist_ok=True)
    software = executable_version(executable)
    database = database_version(executable)
    if database is None:
        raise RuntimeError("AMRFinderPlus database is not installed")

    def process(fasta_path: Path) -> dict[str, Any]:
        genome_id = fasta_path.stem
        output_tsv = output_directory / f"{genome_id}.tsv"
        cached = _is_current(output_tsv, software=software, database=database)
        try:
            if not cached:
                run_nucleotide(
                    executable,
                    fasta_path,
                    output_tsv,
                    organism=organism,
                    threads=threads_per_worker,
                )
            evidence = parse_output(output_tsv, genome_id=genome_id)
            evidence.to_csv(output_tsv.with_suffix(".evidence.csv"), index=False)
            return {
                "genome_id": genome_id,
                "fasta_path": str(fasta_path),
                "status": "cached" if cached else "annotated",
                "element_count": len(evidence),
                "error": "",
            }
        except Exception as error:
            return {
                "genome_id": genome_id,
                "fasta_path": str(fasta_path),
                "status": "failed",
                "element_count": 0,
                "error": f"{type(error).__name__}: {error}",
            }

    records: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(process, path): path for path in fasta_paths}
        for future in as_completed(futures):
            records.append(future.result())

    status = pd.DataFrame(records).sort_values("genome_id").reset_index(drop=True)
    evidence_frames: list[pd.DataFrame] = []
    for genome_id in status.loc[status["status"].ne("failed"), "genome_id"]:
        evidence_path = output_directory / f"{genome_id}.evidence.csv"
        evidence_frames.append(
            pd.read_csv(evidence_path, dtype=object, keep_default_na=False)
        )
    evidence = (
        pd.concat(evidence_frames, ignore_index=True)
        if evidence_frames
        else pd.DataFrame(columns=["genome_id", "feature_key", "feature_value"])
    )
    if evidence.empty:
        features = pd.DataFrame(index=status["genome_id"])
    else:
        features = evidence.pivot_table(
            index="genome_id",
            columns="feature_key",
            values="feature_value",
            aggfunc="max",
            fill_value=0,
        ).astype("uint8")
        features = features.reindex(
            status.loc[status["status"].ne("failed"), "genome_id"], fill_value=0
        )
    features.index.name = "genome_id"
    return status, evidence, features


def write_batch_outputs(
    destination: Path,
    status: pd.DataFrame,
    evidence: pd.DataFrame,
    features: pd.DataFrame,
) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    status.to_csv(destination / "annotation-status.csv", index=False)
    evidence.to_parquet(destination / "amr-evidence.parquet", index=False)
    flat_features = features.reset_index()
    flat_features.to_parquet(destination / "amr-features.parquet", index=False)
    # CSV is also emitted for UI tooling that should not require Arrow at runtime.
    flat_features.to_csv(destination / "amr-features.csv", index=False)
