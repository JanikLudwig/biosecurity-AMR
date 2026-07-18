"""Module 02 — The Predictor: AMR features -> per-antibiotic prediction.

Zero-shot rule engine. For each drug in the panel it decides one of:

    likely to fail  (RESISTANT)   — a known resistance determinant was detected
    likely to work  (SUSCEPTIBLE) — no known determinant, and the target is present
    no-call         (NO_CALL)     — evidence weak/conflicting, screening shallow,
                                     target absent/unknown, or sample out of scope

Key properties required by the brief:
  * a deterministic **target-presence gate** (never call "works" from absence alone);
  * an explicit **no-call** rather than forced yes/no;
  * an **evidence category** (known determinant vs statistical vs none);
  * confidence that is honest about being rule-based (v0 does no training, so it
    never claims category (ii) "statistical association").

There is no train/test split and therefore no leakage risk — a genuine advantage
of the zero-shot baseline. (A dedup utility is still provided for fair evaluation
of a genome collection; see ``genome_firewall.dedup``.)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from . import SAFETY_NOTICE, __version__
from .annotate import AmrHit, AnnotationResult
from . import knowledge as kb

# Prediction call labels (internal enum + human phrasing from the brief).
RESISTANT = "likely to fail"
SUSCEPTIBLE = "likely to work"
NO_CALL = "no-call"

# Evidence categories from the brief.
EV_KNOWN = "known_resistance_determinant"       # (i)
EV_STATISTICAL = "statistical_association"       # (ii) — not produced by v0
EV_NONE = "no_known_resistance_signal"           # (iii)


@dataclass
class PredictConfig:
    """Tunable thresholds for the zero-shot decision rule."""

    resistant_threshold: float = 0.60      # p(fail) at/above -> call RESISTANT
    resistant_nocall_floor: float = 0.40   # weak determinant band -> NO_CALL
    susceptible_base: float = 0.85         # confidence ceiling for "no marker found"
    susceptible_threshold: float = 0.50    # below -> screening too shallow -> NO_CALL


@dataclass
class DrugPrediction:
    drug_id: str
    drug_name: str
    drug_class: str
    call: str                              # RESISTANT / SUSCEPTIBLE / NO_CALL label
    confidence: float                      # rule-based heuristic in [0,1]
    evidence_category: str
    target_status: str                     # present / absent / unknown
    supporting_markers: List[Dict[str, object]] = field(default_factory=list)
    rationale: str = ""
    no_call_reason: Optional[str] = None
    confidence_basis: str = "rule_based_heuristic_uncalibrated"

    def as_dict(self) -> Dict[str, object]:
        return {
            "drug_id": self.drug_id,
            "drug_name": self.drug_name,
            "drug_class": self.drug_class,
            "call": self.call,
            "confidence": round(self.confidence, 3),
            "evidence_category": self.evidence_category,
            "target_status": self.target_status,
            "supporting_markers": self.supporting_markers,
            "rationale": self.rationale,
            "no_call_reason": self.no_call_reason,
            "confidence_basis": self.confidence_basis,
        }


@dataclass
class SamplePrediction:
    species: str
    species_supported: bool
    annotation_backend: str
    screening_completeness: float
    predictions: List[DrugPrediction]
    qc: Optional[Dict[str, object]] = None
    warnings: List[str] = field(default_factory=list)
    version: str = __version__
    safety_notice: str = SAFETY_NOTICE

    @property
    def no_call_rate(self) -> float:
        if not self.predictions:
            return 0.0
        n = sum(1 for p in self.predictions if p.call == NO_CALL)
        return round(n / len(self.predictions), 3)

    def as_dict(self) -> Dict[str, object]:
        return {
            "version": self.version,
            "species": self.species,
            "species_supported": self.species_supported,
            "annotation_backend": self.annotation_backend,
            "screening_completeness": self.screening_completeness,
            "no_call_rate": self.no_call_rate,
            "qc": self.qc,
            "warnings": list(self.warnings),
            "predictions": [p.as_dict() for p in self.predictions],
            "safety_notice": self.safety_notice,
        }


def _noisy_or(weights: List[float]) -> float:
    """Combine independent evidence weights: p = 1 - prod(1 - w_i)."""
    prod = 1.0
    for w in weights:
        prod *= (1.0 - max(0.0, min(1.0, w)))
    return 1.0 - prod


def predict_sample(annotation: AnnotationResult,
                   species: Optional[str] = None,
                   qc: Optional[Dict[str, object]] = None,
                   config: Optional[PredictConfig] = None) -> SamplePrediction:
    """Run the zero-shot predictor over the whole drug panel for one genome."""
    cfg = config or PredictConfig()
    species_ok = kb.species_supported(species)
    qc_failed = bool(qc) and not qc.get("passed", True)

    warnings = list(annotation.warnings)
    if not species_ok:
        warnings.append(
            f"Species '{species}' is outside v0 scope ({kb.supported_species()}); "
            "all drugs returned as no-call."
        )
    if qc_failed:
        warnings.append(
            "Assembly failed QC (" + ", ".join(qc.get("flags", [])) + "); "
            "input is out-of-distribution, all drugs returned as no-call."
        )

    # Pre-compute, per drug, every determinant that implicates it.
    matches_by_drug: Dict[str, list] = {}
    for hit in annotation.hits:
        for m in kb.match_hit(hit):
            matches_by_drug.setdefault(m.drug_id, []).append(m)

    predictions: List[DrugPrediction] = []
    for drug in kb.panel_drugs():
        predictions.append(
            _predict_drug(drug, matches_by_drug.get(drug["id"], []),
                          annotation, species_ok, qc_failed, cfg)
        )

    return SamplePrediction(
        species=species or "unspecified",
        species_supported=species_ok,
        annotation_backend=annotation.backend,
        screening_completeness=annotation.screening_completeness,
        predictions=predictions,
        qc=qc,
        warnings=warnings,
    )


def _marker_dict(m) -> Dict[str, object]:
    return {
        "gene": m.hit.gene,
        "type": "point_mutation" if m.hit.is_point_mutation else "acquired_gene",
        "drug_class": m.hit.drug_class,
        "subclass": m.hit.subclass,
        "identity": m.hit.identity,
        "coverage": m.hit.coverage,
        "source": m.hit.source,
        "rule_kind": m.rule_kind,
        "weight": m.weight,
        "note": m.note,
        "component": m.component,
    }


def _predict_drug(drug, matches, annotation, species_ok, qc_failed,
                  cfg: PredictConfig) -> DrugPrediction:
    drug_id, name, dclass = drug["id"], drug["name"], drug["drug_class"]

    # Sample-level gates first: out-of-scope or failed QC -> honest no-call.
    if qc_failed or not species_ok:
        reason = ("assembly failed QC" if qc_failed
                  else f"species outside supported scope ({kb.supported_species()})")
        return DrugPrediction(
            drug_id=drug_id, drug_name=name, drug_class=dclass, call=NO_CALL,
            confidence=0.0, evidence_category=EV_NONE, target_status="unknown",
            rationale=f"No-call: {reason}. Prediction withheld to avoid false confidence.",
            no_call_reason=reason,
        )

    target_status = kb.target_present(drug, annotation.hits, species_ok)

    # Case A: a resistance determinant was detected -> evidence category (i).
    if matches:
        weights = [m.weight for m in matches]
        p_fail = _noisy_or(weights)
        markers = [_marker_dict(m) for m in matches]
        gene_list = ", ".join(sorted({m.hit.gene for m in matches}))

        if p_fail >= cfg.resistant_threshold:
            return DrugPrediction(
                drug_id=drug_id, drug_name=name, drug_class=dclass, call=RESISTANT,
                confidence=p_fail, evidence_category=EV_KNOWN,
                target_status=target_status, supporting_markers=markers,
                rationale=(f"Known resistance determinant(s) detected ({gene_list}) "
                           f"consistent with {name} resistance."),
            )
        # A determinant exists but the combined evidence is weak/low-level.
        return DrugPrediction(
            drug_id=drug_id, drug_name=name, drug_class=dclass, call=NO_CALL,
            confidence=p_fail, evidence_category=EV_KNOWN,
            target_status=target_status, supporting_markers=markers,
            rationale=(f"A determinant was detected ({gene_list}) but the combined "
                       f"evidence is weak/low-level (p={p_fail:.2f}); confirm by testing."),
            no_call_reason="weak_or_low_level_resistance_evidence",
        )

    # Case B: no determinant found -> target gate decides work vs no-call.
    if target_status == "absent":
        return DrugPrediction(
            drug_id=drug_id, drug_name=name, drug_class=dclass, call=NO_CALL,
            confidence=0.0, evidence_category=EV_NONE, target_status=target_status,
            rationale=(f"Molecular target of {name} not detected; cannot assert the "
                       "drug will work. Withheld."),
            no_call_reason="drug_target_absent",
        )
    if target_status == "unknown":
        return DrugPrediction(
            drug_id=drug_id, drug_name=name, drug_class=dclass, call=NO_CALL,
            confidence=0.0, evidence_category=EV_NONE, target_status=target_status,
            rationale=(f"No resistance determinant found, but presence of the {name} "
                       "target could not be confirmed; withheld."),
            no_call_reason="drug_target_unconfirmed",
        )

    # Target present + no marker -> "likely to work", category (iii). Confidence
    # is bounded by how completely we screened (single tool vs cAMRah consensus):
    # absence of evidence is weaker than evidence of absence.
    conf = cfg.susceptible_base * annotation.screening_completeness
    if conf < cfg.susceptible_threshold:
        return DrugPrediction(
            drug_id=drug_id, drug_name=name, drug_class=dclass, call=NO_CALL,
            confidence=conf, evidence_category=EV_NONE, target_status=target_status,
            rationale=("No known resistance determinant found, but screening breadth "
                       f"is too limited (completeness={annotation.screening_completeness:.2f}) "
                       "to call susceptible; withheld."),
            no_call_reason="screening_too_shallow_to_confirm_susceptible",
        )
    return DrugPrediction(
        drug_id=drug_id, drug_name=name, drug_class=dclass, call=SUSCEPTIBLE,
        confidence=conf, evidence_category=EV_NONE, target_status=target_status,
        rationale=("No known resistance determinant detected and drug target present. "
                   "Absence of a marker is weaker evidence than a positive finding — "
                   "confidence is bounded accordingly."),
    )
