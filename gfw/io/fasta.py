"""FASTA reading + lightweight assembly QC.

The pipeline starts from a *reconstructed, quality-checked* assembly (the brief
places sample-to-genome processing out of scope). A S. aureus draft assembly is
normally many contigs — see ``# report.md`` — so multiple ``>`` records is
expected and is not an error; only mixed organisms would be.
"""

from __future__ import annotations

import gzip
import os
from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional, Tuple

from ..config import QC, QCConfig, GENOMES_DIR


@dataclass
class Contig:
    name: str
    seq: str

    def __len__(self) -> int:  # nucleotides
        return len(self.seq)


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def read_fasta(path: str) -> List[Contig]:
    """Parse a (optionally gzipped) FASTA into contigs. Stdlib only."""
    contigs: List[Contig] = []
    name: Optional[str] = None
    chunks: List[str] = []
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    contigs.append(Contig(name, "".join(chunks)))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        contigs.append(Contig(name, "".join(chunks)))
    return contigs


def genome_path(genome_id: str) -> str:
    """Resolve a genome_id to its assembly path in the repo ``genomes/`` dir."""
    p = os.path.join(GENOMES_DIR, f"{genome_id}.fna")
    if not os.path.exists(p):
        p_gz = p + ".gz"
        if os.path.exists(p_gz):
            return p_gz
    return p


@dataclass
class AssemblyStats:
    n_contigs: int
    total_length: int
    n50: int
    gc: float

    def as_dict(self) -> Dict[str, object]:
        return {"n_contigs": self.n_contigs, "total_length": self.total_length,
                "n50": self.n50, "gc": round(self.gc, 4)}


def assembly_stats(contigs: List[Contig]) -> AssemblyStats:
    lengths = sorted((len(c) for c in contigs), reverse=True)
    total = sum(lengths)
    gc = 0
    for c in contigs:
        s = c.seq.upper()
        gc += s.count("G") + s.count("C")
    n50 = 0
    acc = 0
    for L in lengths:
        acc += L
        if acc >= total / 2:
            n50 = L
            break
    return AssemblyStats(
        n_contigs=len(contigs), total_length=total, n50=n50,
        gc=(gc / total) if total else 0.0,
    )


@dataclass
class QCResult:
    passed: bool
    flags: List[str] = field(default_factory=list)
    stats: Optional[Dict[str, object]] = None
    checkm: Optional[Dict[str, object]] = None

    def as_dict(self) -> Dict[str, object]:
        return {"passed": self.passed, "flags": list(self.flags),
                "stats": self.stats, "checkm": self.checkm}


def qc_assembly(stats: AssemblyStats,
                checkm_completeness: Optional[float] = None,
                checkm_contamination: Optional[float] = None,
                cfg: QCConfig = QC) -> QCResult:
    """Gate an assembly on contiguity, size, and (if available) CheckM.

    A failed QC makes the sample out-of-distribution; the predictor returns
    ``no-call`` for every drug rather than guessing.
    """
    flags: List[str] = []
    if stats.n_contigs > cfg.max_contigs:
        flags.append(f"too_fragmented(n_contigs={stats.n_contigs}>{cfg.max_contigs})")
    if stats.total_length < cfg.min_length:
        flags.append(f"assembly_too_short({stats.total_length}<{cfg.min_length})")
    if stats.total_length > cfg.max_length:
        flags.append(f"assembly_too_long({stats.total_length}>{cfg.max_length})")
    if checkm_completeness is not None and checkm_completeness < cfg.min_completeness:
        flags.append(f"low_completeness({checkm_completeness}<{cfg.min_completeness})")
    if checkm_contamination is not None and checkm_contamination > cfg.max_contamination:
        flags.append(f"high_contamination({checkm_contamination}>{cfg.max_contamination})")
    checkm = None
    if checkm_completeness is not None or checkm_contamination is not None:
        checkm = {"completeness": checkm_completeness, "contamination": checkm_contamination}
    return QCResult(passed=not flags, flags=flags, stats=stats.as_dict(), checkm=checkm)
