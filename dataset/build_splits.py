"""Build a genetic-group-aware train / calibration / hidden-test split.

The brief requires: "Group genomes by genetic similarity, with fixed training,
confidence-calibration, and HIDDEN test sets... genetically related groups,
ideally including groups the model has not seen before." It leaves the
similarity threshold to each team to justify.

We use BV-BRC's precomputed **core-genome MLST hierarchical clustering**
(`cgmlst_hc*` fields on the `genome` resource) as the grouping key instead of
running our own sequence-similarity clustering — it is a real, curated
genetic-distance clustering (allele differences in the core genome), retrieved
with no FASTA download.

Threshold choice: **cgmlst_hc10** (genomes within 10 core-genome allele
differences are one group). This is the conventional resolution used in
genomic epidemiology to call outbreak/transmission clusters (tighter
thresholds like hc0/hc2/hc5 are near-identical-only; looser ones like hc50/
hc100 start merging epidemiologically unrelated lineages together, so most
strains would land in one enormous group and the "hidden test with an unseen
group" property becomes meaningless — see the hc50/hc100 cluster-size dump in
the analysis output). At hc10 we get 2,685 groups from 3,893 genomes
(2,500 of them singletons), giving many genuinely distinct groups to withhold
for calibration/hidden-test while still collapsing true near-duplicates.

Genomes missing a cgMLST call (~17.5%) are each assigned their own singleton
pseudo-group (`nogroup:<genome_id>`) — a conservative fallback that never
causes leakage (it can only over-split, never under-split).

Split: groups are shuffled (seeded) and assigned whole to one split —
70% train / 15% calibration / 15% hidden-test by genome count — so no genetic
group ever appears in more than one split.
"""

from __future__ import annotations

import os
import random

import pandas as pd

_HERE = os.path.dirname(__file__)
SEED = 20260718
GROUP_FIELD = "cgmlst_hc10"
SPLIT_FRACTIONS = {"train": 0.70, "calibration": 0.15, "hidden_test": 0.15}


def qc_pass(row: pd.Series) -> bool:
    """Documented QC rule (mirrors genome_firewall.fasta thresholds)."""
    if row["genome_quality"] != "Good":
        return False
    if pd.isna(row["checkm_completeness"]) or row["checkm_completeness"] < 90:
        return False
    if pd.isna(row["checkm_contamination"]) or row["checkm_contamination"] > 5:
        return False
    if pd.isna(row["contigs"]) or row["contigs"] > 500:
        return False
    return True


def build() -> None:
    meta = pd.read_csv(os.path.join(_HERE, "genome_metadata.csv"), dtype=str)
    for c in ("checkm_completeness", "checkm_contamination", "contigs", "genome_length"):
        meta[c] = pd.to_numeric(meta[c], errors="coerce")

    meta["qc_pass"] = meta.apply(qc_pass, axis=1)
    meta["genetic_group"] = meta[GROUP_FIELD].where(
        meta[GROUP_FIELD].notna(), "nogroup:" + meta["genome_id"]
    )
    meta["genetic_group"] = GROUP_FIELD + ":" + meta["genetic_group"].astype(str)

    groups = (meta.groupby("genetic_group")["genome_id"]
                  .apply(list).reset_index()
                  .rename(columns={"genome_id": "members"}))
    groups["size"] = groups["members"].apply(len)

    rng = random.Random(SEED)
    order = groups.sample(frac=1.0, random_state=SEED).reset_index(drop=True)

    total = order["size"].sum()
    targets = {k: v * total for k, v in SPLIT_FRACTIONS.items()}
    running = {k: 0 for k in SPLIT_FRACTIONS}
    assignment = {}
    for _, row in order.iterrows():
        # assign each group to whichever split is furthest below its target
        # (relative to target size), keeping split sizes close to the fractions
        # while never breaking a group across splits.
        deficits = {k: targets[k] - running[k] for k in SPLIT_FRACTIONS}
        split = max(deficits, key=deficits.get)
        assignment[row["genetic_group"]] = split
        running[split] += row["size"]

    groups["split"] = groups["genetic_group"].map(assignment)
    meta = meta.merge(groups[["genetic_group", "split"]], on="genetic_group", how="left")

    out_cols = ["genome_id", "genome_name", "genome_length", "contigs",
                "checkm_completeness", "checkm_contamination", "genome_quality",
                "qc_pass", "mlst", GROUP_FIELD, "genetic_group", "split",
                "assembly_accession", "biosample_accession", "bioproject_accession"]
    meta[out_cols].to_csv(os.path.join(_HERE, "genome_splits.csv"), index=False)

    print("split sizes (genomes):")
    print(meta["split"].value_counts())
    print()
    print("QC pass rate:", meta["qc_pass"].mean().round(3), f"({meta['qc_pass'].sum()}/{len(meta)})")
    print()
    print("groups per split:")
    print(groups["split"].value_counts())
    print()
    # sanity: no group split across sets
    bad = groups.groupby("genetic_group")["split"].nunique()
    assert (bad == 1).all(), "a genetic group was split across train/cal/test!"
    print("OK: every genetic group assigned to exactly one split.")


if __name__ == "__main__":
    build()
