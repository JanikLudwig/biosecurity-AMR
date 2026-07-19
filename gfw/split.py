"""Leakage-safe train / calibration / test split.

Two levels, matching the brief's "strong submission" definition:

1. **Dedup** — collapse near-identical genomes with the provided **hc10** cgMLST
   clusters (≤10-allele differences), keeping one representative each. Stops the
   model from memorising near-duplicate assemblies.
2. **Grouped split** — assign whole **MLST sequence-type** lineages to one
   partition only, so the hidden test contains clonal groups never seen in
   training. (hc10 is too fine here — nearly every genome is its own hc10 cluster
   — so MLST ST is the honest grouping unit for generalization.)

Deterministic: groups are ordered by a seeded shuffle then greedily packed to hit
the target proportions, so the same split is reproduced every run.
"""

from __future__ import annotations

import random
from typing import Dict, List, Optional

import pandas as pd

from .config import SPLIT, SplitConfig, SPLIT_CSV, ensure_dirs
from .io.labels import load_lab_labels, load_manifest


def _genome_table(labels: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """One row per labelled genome: genome_id, mlst_group, dedup cluster."""
    lab = labels if labels is not None else load_lab_labels()
    g = lab[["genome_id", "mlst_group"]].drop_duplicates("genome_id").copy()
    man = load_manifest()[["genome_id", "cluster"]]
    g = g.merge(man, on="genome_id", how="left")
    # Genomes with no MLST call become their own singleton group (never leak).
    missing = g["mlst_group"].isna()
    g.loc[missing, "mlst_group"] = "STsolo_" + g.loc[missing, "genome_id"].astype(str)
    # Genomes with no dedup cluster are treated as unique (kept).
    g["cluster"] = g["cluster"].fillna("clu_" + g["genome_id"].astype(str))
    return g


def dedup_representatives(g: pd.DataFrame) -> pd.DataFrame:
    """Keep one representative genome per hc10 cluster (deterministic)."""
    return (g.sort_values("genome_id")
             .drop_duplicates(subset="cluster", keep="first")
             .reset_index(drop=True))


def _greedy_group_assignment(sizes: Dict[str, int], cfg: SplitConfig) -> Dict[str, str]:
    """Assign whole groups to partitions to approach the target proportions."""
    total = sum(sizes.values())
    targets = {"train": cfg.train * total,
               "calibration": cfg.calibration * total,
               "test": cfg.test * total}
    filled = {k: 0.0 for k in targets}
    assignment: Dict[str, str] = {}

    groups = list(sizes.keys())
    random.Random(cfg.seed).shuffle(groups)
    # Largest groups first so they don't overshoot a small partition late.
    groups.sort(key=lambda gid: -sizes[gid])
    for gid in groups:
        # Partition with the largest remaining quota (target - filled).
        part = max(targets, key=lambda p: targets[p] - filled[p])
        assignment[gid] = part
        filled[part] += sizes[gid]
    return assignment


def make_split(cfg: SplitConfig = SPLIT,
               labels: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Return genome_id -> partition after dedup + MLST-grouped assignment.

    Columns: genome_id, mlst_group, cluster, partition, is_representative.
    Only representative genomes are assigned to a partition; deduped-out genomes
    are marked ``partition == "dropped_dedup"``.
    """
    g = _genome_table(labels)
    reps = dedup_representatives(g)
    sizes = reps.groupby("mlst_group").size().to_dict()
    assignment = _greedy_group_assignment(sizes, cfg)

    g["is_representative"] = g["genome_id"].isin(set(reps["genome_id"]))
    g["partition"] = g["mlst_group"].map(assignment)
    g.loc[~g["is_representative"], "partition"] = "dropped_dedup"
    return g.reset_index(drop=True)


def save_split(path: str = SPLIT_CSV, cfg: SplitConfig = SPLIT) -> pd.DataFrame:
    ensure_dirs()
    split = make_split(cfg)
    split.to_csv(path, index=False)
    return split


def load_split(path: str = SPLIT_CSV) -> pd.DataFrame:
    return pd.read_csv(path, dtype={"genome_id": str})


def summarize(split: pd.DataFrame) -> str:
    reps = split[split["is_representative"]]
    lines = [f"Genomes: {len(split)} labelled | {len(reps)} representatives "
             f"after hc10 dedup | {split['mlst_group'].nunique()} MLST groups"]
    for part in ("train", "calibration", "test"):
        sub = reps[reps["partition"] == part]
        lines.append(f"  {part:12s}: {len(sub):5d} genomes | "
                     f"{sub['mlst_group'].nunique():4d} groups")
    return "\n".join(lines)
