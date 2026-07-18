"""Build a model-ready AMR gene presence/absence feature matrix.

Input : bvbrc_data/staph_aureus_sp_gene_amr.csv (raw sp_gene "Antibiotic
        Resistance" hits — many rows per genome, since BV-BRC merges calls
        from multiple tools: AMRFinderPlus, RGI/CARD, PATRIC k-mer search).
Output: dataset/features_gene_presence.csv — one row per genome, one binary
        column per normalized gene family (Module 01 "AI features" output).
        dataset/model_table.csv — final table joining features + labels +
        QC + genetic-group split, ready for the Module 02 baseline model.

Gene normalization: RGI/CARD emits allele-numbered variants of the same gene
family (e.g. "tetM_1", "norB_3"); we strip the trailing "_<n>" so all alleles
of one family collapse to one presence feature (documented simplification —
a v1 could keep alleles separate for point-mutation-level resolution).
Rows with no usable gene symbol (~rare CARD "=>" product-only calls) are kept
under a synthetic id extracted from the product string, or dropped if none.

Rare-gene filtering: gene families detected in fewer than MIN_SUPPORT genomes
are dropped — at n≈3,900 genomes, a family seen in 1-4 genomes carries almost
no statistical signal for a per-antibiotic model and mainly adds noise/overfit
risk (documented threshold, tune freely).
"""

from __future__ import annotations

import os
import re

import pandas as pd

_HERE = os.path.dirname(__file__)
MIN_SUPPORT = 5


def _canonicalize(symbol: str) -> str:
    """Fold naming-convention differences across tools into one key.

    AMRFinderPlus, RGI/CARD and PATRIC's k-mer method disagree on formatting
    for the *same* gene — e.g. "ermA" vs "Erm(A)" vs "ERM(A)", or "BlaZ" vs
    "BlaZ family". Without this, ~5 gene columns fragment what should be one
    strong signal (verified: raw pull had "Erm(A)", "ErmA", "Erm(B)", "Erm(C)",
    "Erm(T)", "Erm(33)" as distinct columns). We lowercase, drop " family",
    strip allele suffixes, and remove punctuation that's purely a formatting
    choice (parentheses/apostrophes) while leaving the alphanumeric content
    that actually distinguishes genes (e.g. aac(6')-Ib stays distinguishable
    from aac(6')-Ie-APH(2'')-Ia after punctuation removal).
    """
    s = symbol.strip().lower()
    s = re.sub(r"\s+family$", "", s)
    s = re.sub(r"[()'’]", "", s)
    s = re.sub(r"_\d+$", "", s)
    s = re.sub(r"\s+", "", s)
    return s


def normalize_gene(row: pd.Series) -> str | None:
    gene = row.get("gene")
    if isinstance(gene, str) and gene.strip():
        return _canonicalize(gene)
    product = row.get("product") or ""
    m = re.search(r"=>\s*([A-Za-z0-9()'\-./]+)", str(product))
    if m:
        return _canonicalize(m.group(1))
    return None


def build() -> None:
    raw = pd.read_csv(os.path.join(_HERE, "..", "bvbrc_data",
                                    "staph_aureus_sp_gene_amr.csv"), dtype=str)
    raw["gene_norm"] = raw.apply(normalize_gene, axis=1)
    n_before = len(raw)
    raw = raw.dropna(subset=["gene_norm"])
    print(f"gene hits: {n_before} raw -> {len(raw)} with a resolvable gene symbol "
          f"({n_before - len(raw)} dropped, no gene/product pattern)")

    # collapse multi-tool duplicate calls: one presence flag per (genome, gene),
    # but keep which tool(s)/evidence supported it for auditability.
    per_pair = (raw.groupby(["genome_id", "gene_norm"])["evidence"]
                    .apply(lambda s: ";".join(sorted(set(s.dropna()))))
                    .reset_index())

    support = per_pair["gene_norm"].value_counts()
    keep_genes = support[support >= MIN_SUPPORT].index
    dropped = support[support < MIN_SUPPORT]
    print(f"gene families: {len(support)} distinct -> {len(keep_genes)} kept "
          f"(support >= {MIN_SUPPORT} genomes), {len(dropped)} rare families dropped")

    per_pair = per_pair[per_pair["gene_norm"].isin(keep_genes)]

    wide = (per_pair.assign(present=1)
                     .pivot_table(index="genome_id", columns="gene_norm",
                                  values="present", fill_value=0))
    wide.columns = [f"gene__{c}" for c in wide.columns]
    wide = wide.reset_index()

    out_features = os.path.join(_HERE, "features_gene_presence.csv")
    wide.to_csv(out_features, index=False)
    print(f"wrote {out_features}  shape={wide.shape}  "
          f"({wide.shape[0]} genomes x {wide.shape[1]-1} gene features)")

    # ---- final joined model table ----
    labels = pd.read_csv(os.path.join(_HERE, "labels_long.csv"), dtype=str)
    splits = pd.read_csv(os.path.join(_HERE, "genome_splits.csv"), dtype=str)
    splits["qc_pass"] = splits["qc_pass"].map({"True": True, "False": False})

    model_table = labels.merge(splits[["genome_id", "split", "qc_pass",
                                       "genetic_group", "genome_quality"]],
                               on="genome_id", how="left")
    model_table = model_table.merge(wide, on="genome_id", how="left")
    gene_cols = [c for c in wide.columns if c.startswith("gene__")]
    model_table[gene_cols] = model_table[gene_cols].fillna(0).astype(int)
    model_table["has_gene_features"] = model_table["genome_id"].isin(wide["genome_id"])

    out_model = os.path.join(_HERE, "model_table.csv")
    model_table.to_csv(out_model, index=False)
    print(f"wrote {out_model}  shape={model_table.shape}")
    print()
    print("rows missing gene features (genome had 0 AMR-gene hits of any kind):",
          (~model_table["has_gene_features"]).sum(), "/", len(model_table))
    print()
    print("split x label counts, per antibiotic:")
    for ab, g in model_table.groupby("antibiotic"):
        ct = pd.crosstab(g["split"], g["final_label"])
        print(f"\n-- {ab} --")
        print(ct)


if __name__ == "__main__":
    build()
