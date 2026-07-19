"""Orchestrator — run one genome through the whole pipeline.

Wires the deterministic modules together for a single genome:
QC/scope gate → **M2** target detection → **M1** feature lookup → **M3**
per-drug probability → **M4** decision → **M5** report. Load-once, predict-many
via the :class:`Engine` class so references and models are read a single time.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Dict, List, Optional

import pandas as pd

from . import SPECIES
from .config import QC
from .decide import decide_drug
from .engine_support import scope_and_qc
from .io.fasta import read_fasta, genome_path, assembly_stats, qc_assembly
from .m1_adapter import load_features
from .panel import DrugEntry, load_panel, build_panel
from .predict import DrugModel, load_models
from .report import SampleReport
from .targets.detector import TargetReference, detect, drug_target_status, load_references


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

    # -- feature lookup (M1) ------------------------------------------------
    def _feature_dict(self, genome_id: str) -> Dict[str, int]:
        if self.features is None or genome_id not in self.features.index:
            return {}
        row = self.features.loc[genome_id]
        return {g: int(v) for g, v in row.items() if int(v) != 0}

    # -- full pipeline for one genome --------------------------------------
    def predict_genome(self, genome_id: str,
                       fasta_path: Optional[str] = None) -> SampleReport:
        path = fasta_path or genome_path(genome_id)
        contigs = read_fasta(path)
        scope_ok, qc, species = scope_and_qc(genome_id, contigs)

        detection = detect(genome_id, contigs, self.references)
        feat = self._feature_dict(genome_id)
        present_genes = list(feat.keys())

        decisions = []
        for entry in self.panel:
            dt = drug_target_status(entry.drug, detection)
            model = self.models.get(entry.drug)
            p_r = model.predict_proba(feat) if (model and entry.modelable) else None
            decisions.append(decide_drug(
                entry, p_resistant=p_r, target_status=dt.target_status,
                present_genes=present_genes, target_evidence=dt.detected,
                scope_ok=scope_ok, qc_ok=qc["passed"]))

        return SampleReport(
            genome_id=genome_id, species=species or SPECIES, scope_ok=scope_ok,
            qc=qc, decisions=decisions, n_proteins=detection.n_proteins,
            features_synthetic=self.features_synthetic,
            warnings=detection.warnings + ([] if feat else
                     ["no M1 features for this genome; resistance model saw an all-zero vector"]))


def _safe_panel() -> List[DrugEntry]:
    try:
        return load_panel()
    except Exception:
        return build_panel()
