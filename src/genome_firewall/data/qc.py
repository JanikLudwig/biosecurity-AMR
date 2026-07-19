from __future__ import annotations

import hashlib
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

DNA = frozenset("ACGTNRYKMSWBDHV")
DNA_BYTES = b"ACGTNRYKMSWBDHV"
CANONICAL_DNA_BYTES = b"ACGT"


@dataclass(frozen=True)
class FastaMetrics:
    genome_length: int
    contigs: int
    contig_n50: int
    ambiguous_bases: int
    ambiguous_fraction: float
    sha256: str

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def inspect_fasta(path: Path) -> FastaMetrics:
    """Validate an assembled FASTA and compute small, dependency-free QC metrics."""
    lengths: list[int] = []
    ambiguous = 0
    payload = path.read_bytes()
    digest = hashlib.sha256(payload)
    content = payload.strip()
    if not content:
        raise ValueError("No FASTA records found")
    if not content.startswith(b">"):
        raise ValueError("Sequence found before first FASTA header")
    for record_number, record in enumerate(content.split(b">")[1:], start=1):
        header, separator, body = record.partition(b"\n")
        if not separator or not header.strip():
            raise ValueError(f"FASTA record {record_number} has no header or sequence")
        sequence = b"".join(body.split()).upper()
        if not sequence:
            raise ValueError(f"FASTA record {record_number} is empty")
        invalid = set(sequence.translate(None, DNA_BYTES))
        if invalid:
            symbols = [chr(value) if value < 128 else repr(value) for value in sorted(invalid)]
            raise ValueError(f"Invalid DNA symbols {symbols} in FASTA record {record_number}")
        lengths.append(len(sequence))
        ambiguous += len(sequence.translate(None, CANONICAL_DNA_BYTES))

    genome_length = sum(lengths)
    halfway = genome_length / 2
    cumulative = 0
    n50 = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= halfway:
            n50 = length
            break
    return FastaMetrics(
        genome_length=genome_length,
        contigs=len(lengths),
        contig_n50=n50,
        ambiguous_bases=ambiguous,
        ambiguous_fraction=ambiguous / genome_length,
        sha256=digest.hexdigest(),
    )


def evaluate_quality(
    metrics: FastaMetrics,
    metadata: dict[str, Any],
    quality: dict[str, Any],
) -> list[str]:
    """Return explicit rejection reasons; an empty list means the assembly passes."""
    reasons: list[str] = []
    if metrics.genome_length < quality["min_genome_length"]:
        reasons.append("genome_too_short")
    if metrics.genome_length > quality["max_genome_length"]:
        reasons.append("genome_too_long")
    if metrics.contigs > quality["max_contigs"]:
        reasons.append("too_many_contigs")
    if metrics.contig_n50 < quality["min_n50"]:
        reasons.append("n50_too_low")
    if metrics.ambiguous_fraction > quality["max_ambiguous_fraction"]:
        reasons.append("too_many_ambiguous_bases")

    expected_length = metadata.get("genome_length")
    if expected_length and int(expected_length) != metrics.genome_length:
        reasons.append("download_length_mismatch")
    expected_contigs = metadata.get("contigs")
    if expected_contigs and int(expected_contigs) != metrics.contigs:
        reasons.append("download_contig_count_mismatch")

    completeness = metadata.get("checkm_completeness")
    if completeness is not None and float(completeness) < quality["min_checkm_completeness"]:
        reasons.append("checkm_completeness_too_low")
    contamination = metadata.get("checkm_contamination")
    if contamination is not None and float(contamination) > quality["max_checkm_contamination"]:
        reasons.append("checkm_contamination_too_high")
    allowed = quality.get("allowed_genome_quality", [])
    if allowed and metadata.get("genome_quality") not in allowed:
        reasons.append("unsupported_genome_quality")
    return reasons
