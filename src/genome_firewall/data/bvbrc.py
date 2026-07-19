from __future__ import annotations

import asyncio
import json
from pathlib import Path
from typing import Any

import httpx
import pandas as pd

from genome_firewall.data.qc import evaluate_quality, inspect_fasta


class BvbrcClient:
    """Small public BV-BRC API client with bounded concurrent downloads."""

    def __init__(
        self,
        base_url: str,
        *,
        timeout_seconds: int,
        sequence_result_limit: int,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.sequence_result_limit = sequence_result_limit
        self.client = httpx.AsyncClient(
            timeout=httpx.Timeout(timeout_seconds),
            follow_redirects=True,
            headers={"User-Agent": "genome-firewall/0.1 (defensive AMR research)"},
        )

    async def __aenter__(self) -> "BvbrcClient":
        return self

    async def __aexit__(self, *_: object) -> None:
        await self.client.aclose()

    async def metadata(self, genome_id: str) -> dict[str, Any]:
        response = await self.client.get(
            f"{self.base_url}/genome/{genome_id}", headers={"Accept": "application/json"}
        )
        response.raise_for_status()
        return response.json()

    async def fasta(self, genome_id: str, destination: Path) -> Path:
        destination.parent.mkdir(parents=True, exist_ok=True)
        partial = destination.with_suffix(destination.suffix + ".part")
        url = (
            f"{self.base_url}/genome_sequence/"
            f"?eq(genome_id,{genome_id})&limit({self.sequence_result_limit})"
        )
        async with self.client.stream(
            "GET", url, headers={"Accept": "application/dna+fasta"}
        ) as response:
            response.raise_for_status()
            with partial.open("wb") as handle:
                async for chunk in response.aiter_bytes():
                    handle.write(chunk)
        partial.replace(destination)
        return destination


async def download_and_qc(
    manifest_path: Path,
    output_directory: Path,
    qc_output: Path,
    *,
    species: str,
    taxon_id: int,
    quality: dict[str, Any],
    bvbrc: dict[str, Any],
    limit: int | None = None,
    sample_seed: int | None = None,
) -> pd.DataFrame:
    """Download selected assemblies and record QC/provenance without silent truncation."""
    manifest = pd.read_csv(manifest_path, dtype=object)
    if limit is not None:
        if sample_seed is None:
            manifest = manifest.head(limit)
        else:
            manifest = manifest.sample(
                n=min(limit, len(manifest)), random_state=sample_seed
            ).sort_values("genome_id")
    semaphore = asyncio.Semaphore(bvbrc["download_concurrency"])

    async with BvbrcClient(
        bvbrc["api_base_url"],
        timeout_seconds=bvbrc["timeout_seconds"],
        sequence_result_limit=bvbrc["sequence_result_limit"],
    ) as client:

        async def process(row: dict[str, str]) -> dict[str, object]:
            genome_id = row["genome_id"]
            fasta_path = output_directory / f"{genome_id}.fna"
            metadata_path = output_directory / f"{genome_id}.metadata.json"
            async with semaphore:
                try:
                    metadata = await client.metadata(genome_id)
                    if metadata.get("species") != species or int(metadata.get("taxon_id", -1)) != taxon_id:
                        raise ValueError("BV-BRC metadata does not match the supported species")
                    if not fasta_path.exists():
                        await client.fasta(genome_id, fasta_path)
                    metadata_path.parent.mkdir(parents=True, exist_ok=True)
                    metadata_path.write_text(
                        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
                    )
                    metrics = inspect_fasta(fasta_path)
                    reasons = evaluate_quality(metrics, metadata, quality)
                    return {
                        "genome_id": genome_id,
                        "genome_name": row["genome_name"],
                        **metrics.as_dict(),
                        "passed_qc": not reasons,
                        "rejection_reasons": ";".join(reasons),
                        "download_error": "",
                    }
                except Exception as error:  # capture per-genome failures; do not lose the run
                    return {
                        "genome_id": genome_id,
                        "genome_name": row["genome_name"],
                        "passed_qc": False,
                        "rejection_reasons": "download_or_validation_error",
                        "download_error": f"{type(error).__name__}: {error}",
                    }

        records = await asyncio.gather(
            *(process(row) for row in manifest.to_dict(orient="records"))
        )

    qc = pd.DataFrame(records).sort_values("genome_id")
    qc_output.parent.mkdir(parents=True, exist_ok=True)
    qc.to_csv(qc_output, index=False)
    return qc
