"""Sequence-homology de-duplication for leakage-free *evaluation* (Module 02).

The zero-shot predictor does not train, so there is no train/test split and no
leakage to worry about. But when you evaluate the tool on a *collection* of
genomes, near-identical assemblies would over-count easy cases and inflate the
score (the brief's "weak submission" pattern). This module provides a light,
dependency-free MinHash k-mer sketch + greedy clustering so you can group genomes
by sequence similarity and evaluate one representative (or report grouped
metrics) at a chosen Jaccard threshold.

This is intentionally simple (pure stdlib). For production use a dedicated tool
(e.g. Mash/MinHash, dRep, or the brief's XTree) — the threshold to use is left to
the team to tune and justify.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Dict, List, Sequence, Set, Tuple

from .fasta import parse_fasta


def _kmers(seq: str, k: int) -> Set[str]:
    seq = seq.upper()
    return {seq[i:i + k] for i in range(len(seq) - k + 1)} if len(seq) >= k else set()


def minhash_sketch(fasta_path: str, k: int = 21, num_hashes: int = 200) -> Tuple[int, ...]:
    """Bottom-`num_hashes` MinHash sketch of a genome's k-mers (as sorted ints)."""
    mins: List[int] = []
    seen: Set[int] = set()
    for _, seq in parse_fasta(fasta_path):
        for kmer in _kmers(seq, k):
            h = int(hashlib.blake2b(kmer.encode(), digest_size=8).hexdigest(), 16)
            if h not in seen:
                seen.add(h)
    ordered = sorted(seen)[:num_hashes]
    return tuple(ordered)


def jaccard(a: Sequence[int], b: Sequence[int]) -> float:
    sa, sb = set(a), set(b)
    if not sa and not sb:
        return 1.0
    inter = len(sa & sb)
    union = len(sa | sb)
    return inter / union if union else 0.0


@dataclass
class Cluster:
    representative: str
    members: List[str]


def cluster_genomes(fasta_paths: List[str], k: int = 21, num_hashes: int = 200,
                    threshold: float = 0.9) -> List[Cluster]:
    """Greedy single-linkage-ish clustering of genomes by MinHash Jaccard.

    Two genomes with estimated Jaccard >= ``threshold`` are treated as the same
    (near-identical) and grouped. Returns clusters; use one representative per
    cluster to build a leakage-free evaluation set.
    """
    sketches: Dict[str, Tuple[int, ...]] = {p: minhash_sketch(p, k, num_hashes)
                                            for p in fasta_paths}
    clusters: List[Cluster] = []
    for path in fasta_paths:
        placed = False
        for cluster in clusters:
            if jaccard(sketches[path], sketches[cluster.representative]) >= threshold:
                cluster.members.append(path)
                placed = True
                break
        if not placed:
            clusters.append(Cluster(representative=path, members=[path]))
    return clusters
