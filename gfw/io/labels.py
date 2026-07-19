"""Load the organizer-pinned **laboratory-measured** AMR labels and join them to
the downloaded genomes + their MLST group / QC metadata.

The brief is explicit: use laboratory-measured test results, NOT computational
phenotype fields (which may be model-generated). We therefore keep only rows
whose ``Evidence`` is a laboratory method, and only Resistant/Susceptible calls.
"""

from __future__ import annotations

import os
import re
from functools import lru_cache
from typing import Dict, Optional

import pandas as pd

from ..config import LABELS_TSV, MANIFEST_CSV, METADATA_CSV


# --------------------------------------------------------------------------- #
# Normalization helpers.
# --------------------------------------------------------------------------- #
# Spelling variants in the BV-BRC TSV that denote the same drug.
_DRUG_ALIASES = {
    "ceftarolin": "ceftaroline",
    "phosphomycin": "fosfomycin",
    "co_trimoxazole": "trimethoprim_sulfamethoxazole",
    "cotrimoxazole": "trimethoprim_sulfamethoxazole",
}


def normalize_drug(name: str) -> str:
    """Canonical, filesystem-safe antibiotic id.

    ``"Trimethoprim/sulfamethoxazole"`` -> ``"trimethoprim_sulfamethoxazole"``.
    """
    s = (name or "").strip().lower()
    s = s.replace("+", "/")
    s = re.sub(r"[\s/]+", "_", s)
    s = re.sub(r"[^a-z0-9_]", "", s)
    s = s.strip("_")
    return _DRUG_ALIASES.get(s, s)


def display_drug(name: str) -> str:
    return (name or "").strip().title().replace("Sulfamethoxazole", "sulfamethoxazole")


def normalize_mlst_group(raw: Optional[str]) -> Optional[str]:
    """Collapse inconsistent MLST labels to a stable sequence-type group.

    ``MLST.saureus.22`` and ``MLST.Staphylococcus_aureus.22`` both -> ``ST22``.
    Non-numeric / novel STs keep their raw suffix; missing -> ``None``.
    """
    if raw is None or (isinstance(raw, float) and pd.isna(raw)):
        return None
    tok = str(raw).strip().split(".")[-1]
    if tok == "" or tok.lower() in ("nan", "none", "-", "novel", "nd"):
        return None
    return f"ST{tok}" if tok.isdigit() else f"ST_{tok}"


# --------------------------------------------------------------------------- #
# Loaders (cached — the TSV is ~7 MB).
# --------------------------------------------------------------------------- #
@lru_cache(maxsize=1)
def load_manifest() -> pd.DataFrame:
    man = pd.read_csv(MANIFEST_CSV, dtype=str)
    man["genome_id"] = man["genome_id"].astype(str)
    for c in ("genome_length", "n_contigs", "cluster_size"):
        if c in man.columns:
            man[c] = pd.to_numeric(man[c], errors="coerce")
    return man


@lru_cache(maxsize=1)
def load_metadata() -> pd.DataFrame:
    meta = pd.read_csv(METADATA_CSV, dtype=str)
    meta["genome_id"] = meta["genome_id"].astype(str)
    meta["mlst_group"] = meta["mlst"].map(normalize_mlst_group)
    for c in ("checkm_completeness", "checkm_contamination", "genome_length"):
        if c in meta.columns:
            meta[c] = pd.to_numeric(meta[c], errors="coerce")
    return meta


@lru_cache(maxsize=1)
def load_lab_labels() -> pd.DataFrame:
    """Long table of laboratory R/S calls for downloaded genomes.

    Columns: genome_id, drug (canonical), drug_display, resistant (1/0),
             mlst_group, hc10, checkm_completeness, checkm_contamination,
             genome_quality.
    """
    lab = pd.read_csv(LABELS_TSV, sep="\t", dtype=str)
    lab = lab[lab["Evidence"].str.contains("Laboratory", na=False)].copy()
    lab["genome_id"] = lab["Genome ID"].astype(str)
    lab["drug"] = lab["Antibiotic"].map(normalize_drug)
    lab["drug_display"] = lab["Antibiotic"].map(display_drug)
    phe = lab["Resistant Phenotype"].str.strip().str.lower()
    lab = lab[phe.isin(["resistant", "susceptible"])].copy()
    lab["resistant"] = (lab["Resistant Phenotype"].str.strip().str.lower()
                        == "resistant").astype(int)

    man = load_manifest()
    present = set(man.loc[man["status"] == "ok", "genome_id"])
    lab = lab[lab["genome_id"].isin(present)].copy()

    meta = load_metadata()[
        ["genome_id", "mlst_group", "cgmlst_hc10", "checkm_completeness",
         "checkm_contamination", "genome_quality"]
    ]
    lab = lab.merge(meta, on="genome_id", how="left")
    lab = lab.rename(columns={"cgmlst_hc10": "hc10"})

    # Collapse duplicate genome-drug measurements to a single label:
    # if any laboratory test called Resistant, treat the pair as Resistant
    # (resistance is the safety-relevant, harder-to-miss outcome).
    lab = (lab.sort_values("resistant", ascending=False)
              .drop_duplicates(subset=["genome_id", "drug"], keep="first"))

    keep = ["genome_id", "drug", "drug_display", "resistant", "mlst_group",
            "hc10", "checkm_completeness", "checkm_contamination", "genome_quality"]
    return lab[keep].reset_index(drop=True)
