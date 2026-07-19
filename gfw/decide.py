"""M4 — the only M1/M3--M2 join point.

This module is the sole location that combines deterministic branch outputs:

* **p_R**  — calibrated P(resistant) from the M3 model,
* **target gate** — M2's proof that the drug's molecular target is present,
* **determinant evidence** — did M1 report a known resistance gene for this drug.

M2 evidence is not a feature column, training input, calibration input, or input
to M3's logistic-regression probability. It is used here only to gate a candidate
*likely to work* call; a *likely to fail* call remains M3-driven.

The function returns one of **likely to fail / likely to work / no-call** with a
confidence, an evidence category, and (for a no-call) an explicit reason.

The rule that encodes the challenge's core requirement (`# design notes.md`):
a **likely to work** call is only allowed when the drug's target is *provably
present* — never from the mere absence of resistance markers.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from . import SAFETY_NOTICE
from .config import DECISION, DecisionConfig
from .targets.specs import spec_for

# Calls.
FAIL = "likely to fail"
WORK = "likely to work"
NO_CALL = "no-call"

# Evidence categories from the brief.
EV_KNOWN = "known_resistance_determinant"   # (i)
EV_STATISTICAL = "statistical_association"    # (ii)
EV_NONE = "no_known_resistance_signal"        # (iii)


@dataclass
class DrugDecision:
    drug: str
    drug_display: str
    drug_class: str
    tier: str
    call: str
    confidence: float
    evidence_category: str
    target_status: str
    p_resistant: Optional[float] = None
    supporting_determinants: List[str] = field(default_factory=list)
    target_evidence: List[Dict[str, object]] = field(default_factory=list)
    rationale: str = ""
    no_call_reason: Optional[str] = None

    def as_dict(self) -> Dict[str, object]:
        d = self.__dict__.copy()
        if self.p_resistant is not None:
            d["p_resistant"] = round(self.p_resistant, 4)
        d["confidence"] = round(self.confidence, 4)
        return d


def _determinant_hits(drug: str, present_genes: List[str]) -> List[str]:
    """Which present M1 genes are *known* determinants for this drug (evidence i)."""
    spec = spec_for(drug)
    if not spec:
        return []
    hits = []
    for g in present_genes:
        for pat in spec.determinant_patterns:
            if re.search(pat, g, flags=re.IGNORECASE):
                hits.append(g)
                break
    return hits


def decide_drug(entry,
                p_resistant: Optional[float],
                target_status: str,
                present_genes: List[str],
                target_evidence: Optional[List[Dict[str, object]]] = None,
                scope_ok: bool = True,
                qc_ok: bool = True,
                cfg: DecisionConfig = DECISION) -> DrugDecision:
    """Apply M4 to completed M3 probability and M2 evidence for one drug.

    Callers calculate ``p_resistant`` before calling M4. This function never
    mutates or re-vectorizes M1 data, ensuring M2 cannot alter M3.
    """
    determinants = _determinant_hits(entry.drug, present_genes)
    base = dict(drug=entry.drug, drug_display=entry.drug_display,
                drug_class=entry.drug_class, tier=entry.tier,
                target_status=target_status, p_resistant=p_resistant,
                supporting_determinants=determinants,
                target_evidence=target_evidence or [])

    # --- Hard gates -> no-call ------------------------------------------------
    if not scope_ok:
        return DrugDecision(**base, call=NO_CALL, confidence=0.0,
                            evidence_category=EV_NONE,
                            rationale="Genome is outside the supported species scope.",
                            no_call_reason="out_of_scope_species")
    if not qc_ok:
        return DrugDecision(**base, call=NO_CALL, confidence=0.0,
                            evidence_category=EV_NONE,
                            rationale="Assembly failed QC; input is out-of-distribution.",
                            no_call_reason="failed_assembly_qc")
    if not entry.modelable or p_resistant is None:
        return DrugDecision(**base, call=NO_CALL, confidence=0.0,
                            evidence_category=EV_NONE,
                            rationale=(f"Insufficient laboratory evidence to model "
                                       f"{entry.drug_display} (tier {entry.tier}); "
                                       "no prediction is made."),
                            no_call_reason="drug_not_modelable")

    p = float(p_resistant)
    hi, lo = cfg.p_fail_hi, cfg.p_work_lo
    if entry.tier == "B":                       # widen the no-call band for low-power drugs
        hi = min(0.9, hi + cfg.tier_b_nocall_widen)
        lo = max(0.1, lo - cfg.tier_b_nocall_widen)

    # --- likely to fail: M3 drives this; M2 is report context only -----------
    if p >= hi:
        ev = EV_KNOWN if determinants else EV_STATISTICAL
        if determinants:
            why = (f"Model p(resistant)={p:.2f}; known determinant(s) detected "
                   f"({', '.join(sorted(set(determinants)))}).")
        else:
            why = (f"Model p(resistant)={p:.2f} from a statistical association; "
                   "no catalogued determinant was detected (treat as weaker evidence).")
        return DrugDecision(**base, call=FAIL, confidence=p,
                            evidence_category=ev, rationale=why)

    # --- likely to work: M2 may only downgrade this candidate to no-call -----
    if p <= lo:
        if cfg.require_target_for_work and target_status != "present":
            reason = ("drug_target_absent" if target_status == "absent"
                      else "target_gate_not_applicable")
            note = ("target not detected in the genome" if target_status == "absent"
                    else f"target presence cannot be established ({target_status})")
            return DrugDecision(**base, call=NO_CALL, confidence=1.0 - p,
                                evidence_category=EV_NONE,
                                rationale=(f"Low resistance probability (p={p:.2f}) but "
                                           f"{note}; a 'works' call is withheld so we never "
                                           "assert success from absence of markers alone."),
                                no_call_reason=reason)
        cited = ", ".join(e.get("gene", "?") for e in (target_evidence or [])) or "target"
        return DrugDecision(**base, call=WORK, confidence=1.0 - p,
                            evidence_category=EV_NONE,
                            rationale=(f"No known resistance signal (p={p:.2f}) AND the drug's "
                                       f"molecular target is present in the genome ({cited}) — "
                                       "so the drug has something to act on."))

    # --- honest no-call band --------------------------------------------------
    return DrugDecision(**base, call=NO_CALL, confidence=max(p, 1.0 - p),
                        evidence_category=(EV_KNOWN if determinants else EV_NONE),
                        rationale=(f"Model probability p(resistant)={p:.2f} is in the "
                                   f"uncertain band ({lo:.2f}–{hi:.2f}); returning no-call "
                                   "rather than forcing a yes/no answer."),
                        no_call_reason="uncertain_probability_band")
