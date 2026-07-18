"""FASTA reading and lightweight genome quality checks (Module 01, input side).

Scope note (from the brief): the system starts *after* isolation, sequencing and
genome reconstruction. So we do not assemble or basecall; we only sanity-check a
reconstructed assembly and pass it to the annotator. QC failures become a
sample-level ``no-call`` in the predictor rather than a silent wrong answer.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Tuple

# Nucleotide characters we tolerate in an assembly (IUPAC ambiguity codes + gaps).
_NUC_ALPHABET = set("ACGTNRYSWKMBDHVU-acgtnryswkmbdhvu")


def _open_maybe_gzip(path: str):
    """Open a plain or gzip-compressed text file transparently."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_fasta(path: str) -> Iterator[Tuple[str, str]]:
    """Yield ``(header, sequence)`` records from a (optionally gzipped) FASTA file."""
    header = None
    chunks: List[str] = []
    with _open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


@dataclass
class GenomeQC:
    """Summary statistics + pass/fail flags for a reconstructed assembly."""

    n_contigs: int = 0
    total_length: int = 0
    longest_contig: int = 0
    n50: int = 0
    fraction_ambiguous: float = 0.0
    gc_fraction: float = 0.0
    invalid_char_fraction: float = 0.0
    passed: bool = True
    flags: List[str] = field(default_factory=list)

    def as_dict(self) -> Dict[str, object]:
        return {
            "n_contigs": self.n_contigs,
            "total_length": self.total_length,
            "longest_contig": self.longest_contig,
            "n50": self.n50,
            "fraction_ambiguous": round(self.fraction_ambiguous, 4),
            "gc_fraction": round(self.gc_fraction, 4),
            "invalid_char_fraction": round(self.invalid_char_fraction, 4),
            "passed": self.passed,
            "flags": list(self.flags),
        }


# Default QC thresholds tuned for a typical single bacterial isolate assembly.
# These are deliberately permissive: their only job is to catch obviously broken
# or non-bacterial input, which should become a sample-level no-call.
DEFAULT_QC = {
    "min_total_length": 1_000_000,   # bacterial genomes are megabases; below this = suspect
    "max_total_length": 15_000_000,  # far above this suggests contamination / multiple genomes
    "max_contigs": 2000,             # extreme fragmentation = poor assembly
    "max_fraction_ambiguous": 0.03,  # >3% Ns indicates low-quality reconstruction
    "max_invalid_char_fraction": 0.001,
}


def compute_qc(path: str, thresholds: Dict[str, float] | None = None) -> GenomeQC:
    """Compute assembly QC statistics and flag likely-unreliable input.

    Returns a :class:`GenomeQC`. ``passed=False`` means the predictor should treat
    the whole sample as out-of-distribution and return no-call for every drug.
    """
    t = dict(DEFAULT_QC)
    if thresholds:
        t.update(thresholds)

    lengths: List[int] = []
    n_ambiguous = 0
    n_invalid = 0
    n_gc = 0
    total = 0

    for _, seq in parse_fasta(path):
        lengths.append(len(seq))
        for ch in seq:
            total += 1
            up = ch.upper()
            if up in ("N", "-"):
                n_ambiguous += 1
            elif up in ("G", "C"):
                n_gc += 1
            elif up not in ("A", "T", "U"):
                if ch not in _NUC_ALPHABET:
                    n_invalid += 1
                else:
                    n_ambiguous += 1  # other IUPAC ambiguity codes

    qc = GenomeQC()
    qc.n_contigs = len(lengths)
    qc.total_length = total
    qc.longest_contig = max(lengths) if lengths else 0
    qc.n50 = _n50(lengths)
    qc.fraction_ambiguous = (n_ambiguous / total) if total else 1.0
    qc.invalid_char_fraction = (n_invalid / total) if total else 1.0
    # GC over unambiguous bases only.
    unambiguous = total - n_ambiguous - n_invalid
    qc.gc_fraction = (n_gc / unambiguous) if unambiguous else 0.0

    if qc.n_contigs == 0 or qc.total_length == 0:
        qc.flags.append("empty_or_unreadable_fasta")
    if qc.total_length < t["min_total_length"]:
        qc.flags.append("assembly_too_short_for_bacterial_genome")
    if qc.total_length > t["max_total_length"]:
        qc.flags.append("assembly_too_long_possible_contamination")
    if qc.n_contigs > t["max_contigs"]:
        qc.flags.append("excessive_fragmentation")
    if qc.fraction_ambiguous > t["max_fraction_ambiguous"]:
        qc.flags.append("high_ambiguous_base_fraction")
    if qc.invalid_char_fraction > t["max_invalid_char_fraction"]:
        qc.flags.append("non_nucleotide_characters")

    qc.passed = len(qc.flags) == 0
    return qc


def _n50(lengths: List[int]) -> int:
    """Standard N50: contig length at which half the assembly is in longer contigs."""
    if not lengths:
        return 0
    ordered = sorted(lengths, reverse=True)
    half = sum(ordered) / 2.0
    running = 0
    for length in ordered:
        running += length
        if running >= half:
            return length
    return ordered[-1]
