"""M5 — The Decision Report.

Assembles per-drug M4 decisions into one antibiotic-response report: JSON for the
API/UI and a human-readable Markdown rendering. Every report carries the
mandatory "confirm with standard laboratory testing" notice and, when the AMR
features are the M1 placeholder, a loud synthetic-data banner.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from . import SAFETY_NOTICE, SPECIES, __version__
from .decide import DrugDecision, FAIL, WORK, NO_CALL


@dataclass
class SampleReport:
    genome_id: str
    species: str
    scope_ok: bool
    qc: Dict[str, object]
    decisions: List[DrugDecision]
    n_proteins: int = 0
    features_synthetic: bool = False
    m1_metadata: Dict[str, object] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    version: str = __version__
    safety_notice: str = SAFETY_NOTICE

    def counts(self) -> Dict[str, int]:
        c = {FAIL: 0, WORK: 0, NO_CALL: 0}
        for d in self.decisions:
            c[d.call] = c.get(d.call, 0) + 1
        return c

    def as_dict(self) -> Dict[str, object]:
        return {
            "version": self.version,
            "genome_id": self.genome_id,
            "species": self.species,
            "scope_ok": self.scope_ok,
            "qc": self.qc,
            "n_proteins_predicted": self.n_proteins,
            "features_synthetic": self.features_synthetic,
            "m1": dict(self.m1_metadata),
            "summary": self.counts(),
            "decisions": [d.as_dict() for d in self.decisions],
            "warnings": list(self.warnings),
            "safety_notice": self.safety_notice,
        }

    def to_markdown(self) -> str:
        L: List[str] = []
        L.append(f"# Genome Firewall — antibiotic-response report")
        L.append(f"**Genome:** `{self.genome_id}`  •  **Species:** *{self.species}*  "
                 f"•  **ORFs:** {self.n_proteins}")
        if self.features_synthetic:
            L.append("> ⚠️ **SYNTHETIC M1 FEATURES (placeholder).** Predictions below are "
                     "illustrative of the pipeline only — not real performance. Swap in "
                     "teammates' AMRFinderPlus output to produce real calls.")
        if self.m1_metadata:
            available = self.m1_metadata.get("feature_row_available", False)
            count = self.m1_metadata.get("nonzero_feature_count", 0)
            L.append(f"**M1 AMRFinder feature row:** {'available' if available else 'unavailable'} "
                     f"({count} nonzero feature(s)); M3 probabilities use only this branch.")
        c = self.counts()
        L.append(f"\n**Summary:** {c[FAIL]} likely to fail • {c[WORK]} likely to work • "
                 f"{c[NO_CALL]} no-call")
        L.append("\n| Antibiotic | Call | Confidence | Evidence | Target | Supporting |")
        L.append("|---|---|---|---|---|---|")
        icon = {FAIL: "🔴", WORK: "🟢", NO_CALL: "⚪"}
        for d in self.decisions:
            support = ", ".join(sorted(set(d.supporting_determinants))) or "—"
            if d.call == WORK and d.target_evidence:
                support = "target: " + ", ".join(e.get("gene", "?") for e in d.target_evidence)
            L.append(f"| {d.drug_display} | {icon.get(d.call,'')} {d.call} | "
                     f"{d.confidence:.0%} | {d.evidence_category} | {d.target_status} | {support} |")
        L.append(f"\n---\n> ⚕️ **{self.safety_notice}**")
        return "\n".join(L)
