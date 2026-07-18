"""Build a clean, model-ready label dataset from the raw BV-BRC AMR pull.

Input : bvbrc_data/staph_aureus_amr_lab_method.tsv (raw genome_amr API export,
        evidence == "Laboratory Method" only, taxon_id == 1280 / S. aureus).
Output: dataset/labels_long.csv   — one row per (genome_id, antibiotic) with a
                                     single resolved final_label
        dataset/panel_summary.csv — per-drug genome/class-balance summary
        dataset/DATASET_CARD.md   — provenance + the labeling rule, documented

Labeling rule (documented, per the brief's "one final label for each
genome-antibiotic pair" requirement):
  1. Collapse duplicate lab records for the same (genome, antibiotic).
  2. If all non-blank phenotypes for the pair agree on Resistant/Susceptible
     -> that is the final_label.
  3. If Resistant and Susceptible both occur for the same pair -> "conflicting"
     (excluded from training; kept in the table for transparency).
  4. If the only phenotypes present are Intermediate/Nonsusceptible (CLSI/EUCAST
     categories between clearly-S and clearly-R) -> "intermediate" (excluded
     from a binary train set; may be useful for no-call studies later).
  5. Pairs with no usable phenotype (blank only) are dropped entirely.

This script only touches the raw label table — no genome sequences are
downloaded here (see DATASET_CARD.md "Status / Phase 2").
"""

from __future__ import annotations

import os

import pandas as pd

_HERE = os.path.dirname(__file__)
_RAW = os.path.join(_HERE, "..", "bvbrc_data", "staph_aureus_amr_lab_method.tsv")

# Known raw-string synonyms/typos in the BV-BRC export that refer to one drug.
_ANTIBIOTIC_SYNONYMS = {
    "phosphomycin": "fosfomycin",
    "ceftarolin": "ceftaroline",
}

# The v0 training panel: 5 antibiotics, 5 distinct resistance mechanisms,
# chosen from the candidates by genome coverage (>=1500 genomes each) and
# R:S class balance (see dataset/DATASET_CARD.md for the full comparison
# table and the rationale for each pick).
PANEL = {
    "cefoxitin": "BETA-LACTAM (mecA/mecC surrogate marker for MRSA)",
    "ciprofloxacin": "FLUOROQUINOLONE (gyrA/grlA)",
    "erythromycin": "MACROLIDE (ermA/ermC/msrA)",
    "gentamicin": "AMINOGLYCOSIDE (aac(6')-aph(2''))",
    "trimethoprim/sulfamethoxazole": "FOLATE-PATHWAY-INHIBITOR (dfrA/dfrG, sul)",
}


def load_raw() -> pd.DataFrame:
    df = pd.read_csv(_RAW, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    for c in df.columns:
        df[c] = df[c].astype(str).str.strip().str.strip('"').replace({"nan": ""})
    df["antibiotic_norm"] = df["Antibiotic"].str.lower().str.strip().replace(_ANTIBIOTIC_SYNONYMS)
    return df


def resolve_pair(phenotypes: pd.Series) -> str:
    vals = set(p for p in phenotypes if p in ("Resistant", "Susceptible",
                                               "Intermediate", "Nonsusceptible"))
    has_r = "Resistant" in vals
    has_s = "Susceptible" in vals
    if has_r and has_s:
        return "conflicting"
    if has_r:
        return "resistant"
    if has_s:
        return "susceptible"
    if vals:  # only Intermediate/Nonsusceptible present
        return "intermediate"
    return "unlabeled"


def build() -> None:
    df = load_raw()
    df = df[df["antibiotic_norm"].isin(PANEL)].copy()

    rows = []
    group_cols = ["Genome ID", "Genome Name", "antibiotic_norm"]
    for (genome_id, genome_name, ab), g in df.groupby(group_cols, sort=False):
        final = resolve_pair(g["Resistant Phenotype"])
        if final == "unlabeled":
            continue
        pubmeds = sorted({p for p in g["PubMed"] if p})
        methods = sorted({m for m in g["Laboratory Typing Method"] if m})
        standards = sorted({s for s in g["Testing Standard"] if s})
        rows.append({
            "genome_id": genome_id,
            "genome_name": genome_name,
            "antibiotic": ab,
            "drug_class": PANEL[ab],
            "final_label": final,
            "n_lab_records": len(g),
            "pubmed_ids": ";".join(pubmeds),
            "laboratory_typing_method": ";".join(methods),
            "testing_standard": ";".join(standards),
        })

    labels = pd.DataFrame(rows).sort_values(["antibiotic", "genome_id"])
    out_labels = os.path.join(_HERE, "labels_long.csv")
    labels.to_csv(out_labels, index=False)

    summary_rows = []
    for ab in PANEL:
        sub = labels[labels["antibiotic"] == ab]
        vc = sub["final_label"].value_counts()
        r, s = vc.get("resistant", 0), vc.get("susceptible", 0)
        bal = min(r, s) / max(r, s) if max(r, s) else 0.0
        summary_rows.append({
            "antibiotic": ab, "drug_class": PANEL[ab],
            "n_genomes_labeled": len(sub),
            "resistant": r, "susceptible": s,
            "intermediate": vc.get("intermediate", 0),
            "conflicting": vc.get("conflicting", 0),
            "R:S_balance": round(bal, 3),
        })
    summary = pd.DataFrame(summary_rows)
    out_summary = os.path.join(_HERE, "panel_summary.csv")
    summary.to_csv(out_summary, index=False)

    print(f"wrote {out_labels}  ({len(labels)} genome-antibiotic rows, "
          f"{labels['genome_id'].nunique()} unique genomes)")
    print(f"wrote {out_summary}")
    print()
    print(summary.to_string(index=False))


if __name__ == "__main__":
    build()
