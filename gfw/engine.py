"""Orchestrator for the deterministic, two-branch decision path.

The branches deliberately stay separate until **M4**:

* ``predict_m1`` looks up teammate-produced AMRFinderPlus features and runs M3
  logistic-regression models. Its probabilities depend on M1 features only.
* ``predict_m2`` runs deterministic target detection from the assembly and
  returns per-drug target evidence. It is never a model feature or calibrator input.
* :meth:`predict_genome` is the sole M4 join point, followed by M5 reporting.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import pandas as pd

from . import SPECIES
from .decide import decide_drug
from .engine_support import scope_and_qc
from .io.fasta import read_fasta, genome_path
from .m1_adapter import load_features
from .panel import DrugEntry, load_panel, build_panel
from .predict import DrugModel, load_models
from .report import SampleReport
from .targets.detector import TargetReference, detect, drug_target_status, load_references


@dataclass
class M1Prediction:
    """Auditable output of Branch A (M1 feature lookup followed by M3)."""

    genome_id: str
    features: Dict[str, int]
    feature_row_available: bool
    feature_count: int
    model_feature_counts: Dict[str, int] = field(default_factory=dict)
    p_resistant: Dict[str, float] = field(default_factory=dict)


@dataclass
class M2Prediction:
    """Auditable output of Branch B (target detection only)."""

    detection: object
    targets: Dict[str, object]


class Engine:
    def __init__(self,
                 models: Optional[Dict[str, DrugModel]] = None,
                 panel: Optional[List[DrugEntry]] = None,
                 references: Optional[Dict[str, TargetReference]] = None,
                 features: Optional[pd.DataFrame] = None,
                 features_synthetic: bool = False):
        self.models = models if models is not None else load_models()
        self.panel = panel if panel is not None else _safe_panel()
        self.references = references if references is not None else load_references()
        if features is None:
            try:
                features, features_synthetic = load_features()
            except Exception:
                features, features_synthetic = None, False
        self.features = features
        self.features_synthetic = features_synthetic

    # -- Branch A: M1 AMRFinder features -> M3 probabilities ----------------
    def _feature_dict(self, genome_id: str) -> Dict[str, int]:
        if self.features is None or genome_id not in self.features.index:
            return {}
        row = self.features.loc[genome_id]
        return {g: int(v) for g, v in row.items() if int(v) != 0}

    def predict_m1(self, genome_id: str) -> M1Prediction:
        """Run Branch A only: M1 feature row -> M3 calibrated probabilities.

        M2 output is intentionally neither accepted nor consulted here: model
        vectorization uses each saved model's fixed M1 feature order.
        """
        feature_row_available = self.features is not None and genome_id in self.features.index
        feat = self._feature_dict(genome_id)
        probabilities: Dict[str, float] = {}
        feature_counts: Dict[str, int] = {}
        for entry in self.panel:
            model = self.models.get(entry.drug)
            if model and entry.modelable:
                feature_counts[entry.drug] = len(model.features)
                probabilities[entry.drug] = model.predict_proba(feat)
        return M1Prediction(
            genome_id=genome_id, features=feat,
            feature_row_available=feature_row_available,
            feature_count=len(feat), model_feature_counts=feature_counts,
            p_resistant=probabilities)

    # -- Branch B: FASTA -> M2 target evidence ------------------------------
    def predict_m2(self, genome_id: str, contigs) -> M2Prediction:
        """Run Branch B only; its output cannot affect any M3 probability."""
        detection = detect(genome_id, contigs, self.references)
        return M2Prediction(
            detection=detection,
            targets={entry.drug: drug_target_status(entry.drug, detection)
                     for entry in self.panel})

    # -- M4 join, then M5 report --------------------------------------------
    def predict_genome(self, genome_id: str,
                       fasta_path: Optional[str] = None) -> SampleReport:
        path = fasta_path or genome_path(genome_id)
        contigs = read_fasta(path)
        scope_ok, qc, species = scope_and_qc(genome_id, contigs)

        m1 = self.predict_m1(genome_id)
        m2 = self.predict_m2(genome_id, contigs)
        present_genes = list(m1.features)

        decisions = []
        for entry in self.panel:
            # M4 is the only M1/M3--M2 join: probabilities were fully computed
            # in Branch A before Branch B evidence is considered.
            dt = m2.targets[entry.drug]
            decisions.append(decide_drug(
                entry, p_resistant=m1.p_resistant.get(entry.drug), target_status=dt.target_status,
                present_genes=present_genes, target_evidence=dt.detected,
                scope_ok=scope_ok, qc_ok=qc["passed"]))

        return SampleReport(
            genome_id=genome_id, species=species or SPECIES, scope_ok=scope_ok,
            qc=qc, decisions=decisions, n_proteins=m2.detection.n_proteins,
            features_synthetic=self.features_synthetic,
            m1_metadata={"feature_row_available": m1.feature_row_available,
                         "nonzero_feature_count": m1.feature_count,
                         "model_feature_counts": m1.model_feature_counts,
                         "source": "AMRFinderPlus feature matrix (M1)"},
            warnings=m2.detection.warnings + ([] if m1.features else
                     ["no M1 features for this genome; resistance model saw an all-zero vector"]))


def _safe_panel() -> List[DrugEntry]:
    try:
        return load_panel()
    except Exception:
        return build_panel()
