"""Loads the zero-shot AMR knowledge base and matches AMR hits to panel drugs.

The knowledge base lives in ``data/antibiotics.json`` (the drug panel + molecular
targets) and ``data/gene_drug_rules.json`` (which determinants confer resistance
to which drug). Keeping it as data — not code — means a domain expert can extend
the panel or rules without touching the engine.
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass
from functools import lru_cache
from typing import Dict, List, Optional

from .annotate import AmrHit

_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


@dataclass
class DrugMatch:
    """A single reason a determinant implicates a drug (one hit may yield several)."""

    drug_id: str
    weight: float          # base rule-based confidence for "likely to fail"
    note: str
    rule_kind: str         # gene_rule | subclass_map | class_map
    hit: AmrHit
    component: Optional[str] = None  # for combination drugs (e.g. trimethoprim/sulfonamide)


@lru_cache(maxsize=1)
def load_panel() -> Dict[str, object]:
    with open(os.path.join(_DATA_DIR, "antibiotics.json")) as fh:
        return json.load(fh)


@lru_cache(maxsize=1)
def load_rules() -> Dict[str, object]:
    with open(os.path.join(_DATA_DIR, "gene_drug_rules.json")) as fh:
        return json.load(fh)


def panel_drugs() -> List[Dict[str, object]]:
    return list(load_panel()["panel"])


def drug_by_id(drug_id: str) -> Optional[Dict[str, object]]:
    for d in load_panel()["panel"]:
        if d["id"] == drug_id:
            return d
    return None


def supported_species() -> str:
    return str(load_panel()["species"])


def species_supported(species: Optional[str]) -> bool:
    """True if ``species`` is (loosely) the supported species. None = unknown -> True.

    We accept ``None`` as "caller did not specify" and do not block on it; the
    predictor still records that species was unconfirmed.
    """
    if not species:
        return True
    s = species.strip().lower()
    return any(tok in s for tok in load_panel()["species_tokens"])


def _quality_scale(hit: AmrHit) -> float:
    """Scale a rule weight by alignment quality when available (identity/coverage).

    Missing values -> 1.0 (no penalty). A strong acquired-gene hit (~100/100)
    keeps full weight; partial/low-identity hits are down-weighted, nudging weak
    evidence toward a no-call.
    """
    ident = hit.identity if hit.identity is not None else 100.0
    cov = hit.coverage if hit.coverage is not None else 100.0
    # Map identity 90..100 -> 0.7..1.0, coverage 60..100 -> 0.7..1.0 (clamped).
    id_factor = min(1.0, max(0.5, (ident - 80.0) / 20.0)) if ident < 100 else 1.0
    cov_factor = min(1.0, max(0.5, (cov - 50.0) / 50.0)) if cov < 100 else 1.0
    return round(min(1.0, id_factor * cov_factor), 4)


def match_hit(hit: AmrHit) -> List[DrugMatch]:
    """Return every drug this single AMR hit implicates, with precedence.

    Precedence: an explicit ``gene_rule`` wins; otherwise fall back to the
    determinant's ``Subclass`` map, then its broad ``Class`` map. Gene symbol,
    subclass and class are all provided by the annotation tool (AMRFinderPlus /
    cAMRah) — this is the honest zero-shot signal.
    """
    rules = load_rules()
    gene = hit.gene or ""
    scale = _quality_scale(hit)

    # 1) Gene-symbol-specific rules (highest precision).
    gene_matches: List[DrugMatch] = []
    for rule in rules["gene_rules"]:
        if _gene_rule_matches(rule, gene):
            for drug_id in rule["drugs"]:
                gene_matches.append(
                    DrugMatch(
                        drug_id=drug_id,
                        weight=round(rule["weight"] * scale, 4),
                        note=rule.get("note", ""),
                        rule_kind="gene_rule",
                        hit=hit,
                        component=rule.get("component"),
                    )
                )
    if gene_matches:
        return gene_matches

    # 2) Subclass map (specific, curated by the annotation tool per allele).
    sub = (hit.subclass or "").upper()
    if sub and sub in rules["subclass_map"]:
        w = rules["default_subclass_weight"] * scale
        return [
            DrugMatch(drug_id=d, weight=round(w, 4),
                      note=f"determinant subclass {sub}", rule_kind="subclass_map", hit=hit)
            for d in rules["subclass_map"][sub]
        ]

    # 3) Class map (broad fallback).
    cls = (hit.drug_class or "").upper()
    if cls and cls in rules["class_map"]:
        w = rules["default_class_weight"] * scale
        return [
            DrugMatch(drug_id=d, weight=round(w, 4),
                      note=f"determinant class {cls}", rule_kind="class_map", hit=hit)
            for d in rules["class_map"][cls]
        ]

    return []


def _gene_rule_matches(rule: Dict[str, object], gene: str) -> bool:
    kind = rule["match"]
    pattern = str(rule["pattern"])
    if kind == "prefix":
        return gene.lower().startswith(pattern.lower())
    if kind == "exact":
        return gene.lower() == pattern.lower()
    if kind == "regex":
        return re.search(pattern, gene, flags=re.IGNORECASE) is not None
    return False


def target_present(drug: Dict[str, object], hits: List[AmrHit],
                   species_ok: bool) -> str:
    """Deterministic molecular-target gate (brief, Module 02).

    Returns ``"present"`` / ``"absent"`` / ``"unknown"``. The gate exists so the
    system never reports *likely to work* purely from the absence of resistance
    markers — a drug can only work if its target exists.

    v0 behaviour: the panel's targets are **essential** genes (gyrase, PBPs,
    ribosome, DHFR/DHPS) that are present in every viable cell of the supported
    species. So for in-scope genomes we return ``"present"``. If species is not
    the supported one, the essentiality assumption no longer holds -> ``"unknown"``.
    The interface is built to accept real target-detection evidence later.
    """
    if drug.get("target_essential", False):
        return "present" if species_ok else "unknown"
    return "unknown"
