from __future__ import annotations

import hashlib
import json
from itertools import combinations
from pathlib import Path
from typing import Iterable

import pandas as pd
import sourmash


def _sequences(path: Path) -> Iterable[str]:
    sequence: list[str] = []
    with path.open(encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if sequence:
                    yield "".join(sequence)
                    sequence = []
            else:
                sequence.append(line)
    if sequence:
        yield "".join(sequence)


def sketch_fasta(path: Path, *, ksize: int, scaled: int) -> sourmash.MinHash:
    sketch = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled)
    for sequence in _sequences(path):
        sketch.add_sequence(sequence, force=True)
    if not len(sketch):
        raise ValueError(f"No valid k-mers found in {path}")
    return sketch


def estimated_ani(left: sourmash.MinHash, right: sourmash.MinHash, *, ksize: int) -> float:
    """Estimate ANI from k-mer Jaccard using the standard Mash transform."""
    jaccard = left.jaccard(right)
    if jaccard <= 0:
        return 0.0
    return float((2 * jaccard / (1 + jaccard)) ** (1 / ksize))


class _DisjointSet:
    def __init__(self, values: Iterable[str]) -> None:
        self.parent = {value: value for value in values}

    def find(self, value: str) -> str:
        root = value
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[value] != value:
            parent = self.parent[value]
            self.parent[value] = root
            value = parent
        return root

    def union(self, left: str, right: str) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            keep, merge = sorted([left_root, right_root])
            self.parent[merge] = keep


def cluster_fastas(
    fasta_paths: list[Path],
    *,
    ksize: int,
    scaled: int,
    ani_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Single-linkage cluster genomes whose pairwise estimated ANI meets the threshold."""
    sketches = {path.stem: sketch_fasta(path, ksize=ksize, scaled=scaled) for path in fasta_paths}
    groups = _DisjointSet(sketches)
    edges: list[dict[str, object]] = []
    for left_id, right_id in combinations(sorted(sketches), 2):
        ani = estimated_ani(sketches[left_id], sketches[right_id], ksize=ksize)
        if ani >= ani_threshold:
            groups.union(left_id, right_id)
            edges.append({"left_genome_id": left_id, "right_genome_id": right_id, "estimated_ani": ani})

    roots = {genome_id: groups.find(genome_id) for genome_id in sketches}
    unique_roots = sorted(set(roots.values()))
    cluster_ids = {root: f"cluster_{index:05d}" for index, root in enumerate(unique_roots, start=1)}
    membership = pd.DataFrame(
        {
            "genome_id": sorted(sketches),
            "cluster_id": [cluster_ids[roots[genome_id]] for genome_id in sorted(sketches)],
        }
    )
    sizes = membership.groupby("cluster_id")["genome_id"].transform("size")
    membership["cluster_size"] = sizes
    return membership, pd.DataFrame(
        edges, columns=["left_genome_id", "right_genome_id", "estimated_ani"]
    )


def assign_grouped_splits(
    membership: pd.DataFrame,
    *,
    train_fraction: float,
    calibration_fraction: float,
    test_fraction: float,
    seed: int,
) -> pd.DataFrame:
    fractions = {
        "train": train_fraction,
        "calibration": calibration_fraction,
        "test": test_fraction,
    }
    if abs(sum(fractions.values()) - 1.0) > 1e-9:
        raise ValueError("Split fractions must sum to 1")

    counts = membership.groupby("cluster_id").size().rename("size").reset_index()
    counts["tie_break"] = counts["cluster_id"].map(
        lambda value: hashlib.sha256(f"{seed}:{value}".encode()).hexdigest()
    )
    counts = counts.sort_values(["size", "tie_break"], ascending=[False, True])
    targets = {name: fraction * len(membership) for name, fraction in fractions.items()}
    assigned = {name: 0 for name in fractions}
    cluster_split: dict[str, str] = {}
    for row in counts.itertuples(index=False):
        group_size = int(row.size)
        split = min(
            fractions,
            key=lambda name: (
                abs((assigned[name] + group_size) - targets[name])
                - abs(assigned[name] - targets[name]),
                {"train": 0, "calibration": 1, "test": 2}[name],
            ),
        )
        cluster_split[row.cluster_id] = split
        assigned[split] += group_size

    result = membership.copy()
    result["split"] = result["cluster_id"].map(cluster_split)
    if result.groupby("cluster_id")["split"].nunique().max() != 1:
        raise AssertionError("A homology cluster crosses dataset splits")
    return result.sort_values(["split", "cluster_id", "genome_id"]).reset_index(drop=True)


def phenotype_split_support(
    splits: pd.DataFrame,
    phenotypes: pd.DataFrame,
) -> pd.DataFrame:
    """Count susceptible/resistant labels and genetic groups in every drug partition."""
    merged = phenotypes.merge(
        splits[["genome_id", "cluster_id", "split"]],
        on="genome_id",
        how="inner",
        validate="many_to_one",
    )
    counts = (
        merged.groupby(["antibiotic", "split", "label"])["genome_id"]
        .nunique()
        .unstack(fill_value=0)
        .reset_index()
    )
    for label in ["Resistant", "Susceptible"]:
        if label not in counts:
            counts[label] = 0
    groups = (
        merged.groupby(["antibiotic", "split"])["cluster_id"]
        .nunique()
        .rename("clusters")
        .reset_index()
    )
    result = counts.merge(groups, on=["antibiotic", "split"], how="left")
    result["minimum_class_count"] = result[["Resistant", "Susceptible"]].min(axis=1)
    return result.sort_values(["antibiotic", "split"]).reset_index(drop=True)


def write_split_outputs(
    destination: Path,
    membership: pd.DataFrame,
    edges: pd.DataFrame,
    parameters: dict[str, object],
) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    membership.to_csv(destination / "genome-splits.csv", index=False)
    edges.to_parquet(destination / "near-duplicate-edges.parquet", index=False)
    summary = (
        membership.groupby("split")
        .agg(genomes=("genome_id", "size"), clusters=("cluster_id", "nunique"))
        .reset_index()
        .to_dict(orient="records")
    )
    payload = {**parameters, "summary": summary}
    (destination / "split-provenance.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
