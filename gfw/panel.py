"""Data-driven antibiotic panel.

The panel is **not** hand-written: every antibiotic that appears in the
laboratory-labelled TSV is included, then sorted into a modelling tier by how
much balanced lab evidence exists for it (`config.PanelConfig`). This keeps the
scope statement honest (brief requirement) and reproducible.

* **Tier A** — both classes ≥ ``tier_a_min_per_class``  → train + calibrate.
* **Tier B** — 20 ≤ minority class < 100                → train but bias to no-call.
* **Tier C** — minority class < 20 (or too few groups)  → structural no-call.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional

import pandas as pd

from .config import PANEL, PanelConfig, PANEL_JSON, ensure_dirs
from .io.labels import load_lab_labels
from .targets.specs import spec_for


@dataclass
class DrugEntry:
    drug: str
    drug_display: str
    tier: str                 # "A" | "B" | "C"
    n_resistant: int
    n_susceptible: int
    n_total: int
    n_groups: int             # distinct MLST groups with a label for this drug
    drug_class: str
    mechanism: str
    target_kind: str          # protein | membrane | cell_wall | unknown
    target_genes: List[str] = field(default_factory=list)
    modelable: bool = False   # Tier A or B -> a model is trained
    note: str = ""

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


def _tier(n_r: int, n_s: int, n_groups: int, cfg: PanelConfig) -> str:
    minority = min(n_r, n_s)
    if n_groups < cfg.min_groups:
        return "C"
    if minority >= cfg.tier_a_min_per_class:
        return "A"
    if minority >= cfg.tier_b_min_per_class:
        return "B"
    return "C"


def build_panel(cfg: PanelConfig = PANEL,
                labels: Optional[pd.DataFrame] = None) -> List[DrugEntry]:
    lab = labels if labels is not None else load_lab_labels()
    entries: List[DrugEntry] = []
    for drug, sub in lab.groupby("drug"):
        n_r = int((sub["resistant"] == 1).sum())
        n_s = int((sub["resistant"] == 0).sum())
        n_groups = int(sub["mlst_group"].dropna().nunique())
        tier = _tier(n_r, n_s, n_groups, cfg)
        spec = spec_for(drug)
        entries.append(DrugEntry(
            drug=drug,
            drug_display=str(sub["drug_display"].iloc[0]),
            tier=tier,
            n_resistant=n_r, n_susceptible=n_s, n_total=n_r + n_s,
            n_groups=n_groups,
            drug_class=spec.drug_class if spec else "unknown",
            mechanism=spec.mechanism if spec else "",
            target_kind=spec.target_kind if spec else "unknown",
            target_genes=list(spec.target_genes) if spec else [],
            modelable=tier in ("A", "B"),
            note=spec.note if spec else "no curated target spec for this drug",
        ))
    # Most-evidence first, Tier A before B before C.
    order = {"A": 0, "B": 1, "C": 2}
    entries.sort(key=lambda e: (order[e.tier], -e.n_total))
    return entries


def modelable_drugs(panel: Optional[List[DrugEntry]] = None) -> List[str]:
    panel = panel or build_panel()
    return [e.drug for e in panel if e.modelable]


def save_panel(path: str = PANEL_JSON, cfg: PanelConfig = PANEL) -> List[DrugEntry]:
    ensure_dirs()
    panel = build_panel(cfg)
    with open(path, "w") as fh:
        json.dump([e.as_dict() for e in panel], fh, indent=2)
    return panel


def load_panel(path: str = PANEL_JSON) -> List[DrugEntry]:
    with open(path) as fh:
        return [DrugEntry(**d) for d in json.load(fh)]
