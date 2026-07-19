from pathlib import Path

import pandas as pd

from genome_firewall.splitting.homology import (
    assign_grouped_splits,
    cluster_fastas,
    phenotype_split_support,
)


def test_identical_genomes_never_cross_splits(tmp_path: Path) -> None:
    sequence = "ACGT" * 100
    paths = []
    for name, content in [("a", sequence), ("b", sequence), ("c", "TGCA" * 100)]:
        path = tmp_path / f"{name}.fna"
        path.write_text(f">contig\n{content}\n", encoding="ascii")
        paths.append(path)

    membership, edges = cluster_fastas(
        paths, ksize=21, scaled=1, ani_threshold=0.99
    )
    assert membership.set_index("genome_id").loc["a", "cluster_id"] == membership.set_index(
        "genome_id"
    ).loc["b", "cluster_id"]
    assert len(edges) == 1

    splits = assign_grouped_splits(
        membership,
        train_fraction=0.6,
        calibration_fraction=0.2,
        test_fraction=0.2,
        seed=42,
    )
    assert splits.groupby("cluster_id")["split"].nunique().max() == 1
    duplicate_cluster = splits.loc[splits["genome_id"].eq("a"), "cluster_id"].iloc[0]
    assert splits.loc[splits["cluster_id"].eq(duplicate_cluster), "split"].iloc[0] == "train"
    assert isinstance(splits, pd.DataFrame)


def test_phenotype_support_counts_classes_and_clusters() -> None:
    splits = pd.DataFrame(
        {
            "genome_id": ["a", "b", "c"],
            "cluster_id": ["one", "one", "two"],
            "split": ["train", "train", "test"],
        }
    )
    phenotypes = pd.DataFrame(
        {
            "genome_id": ["a", "b", "c"],
            "antibiotic": ["drug", "drug", "drug"],
            "label": ["Resistant", "Susceptible", "Susceptible"],
        }
    )
    support = phenotype_split_support(splits, phenotypes).set_index("split")
    assert support.loc["train", "minimum_class_count"] == 1
    assert support.loc["train", "clusters"] == 1
    assert support.loc["test", "Resistant"] == 0
