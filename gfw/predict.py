"""M3 — The Predictor.

One **regularized logistic regression per antibiotic** (the brief's recommended
baseline) over the M1 AMR feature matrix, wrapped in a probability
**calibrator** fitted on a held-out calibration split. Per-drug L2 coefficients
are the model's own explanation of which genes drove a call.

Training uses only the ``train`` partition; calibration only the ``calibration``
partition; neither ever sees the grouped hidden ``test`` set (see
:mod:`gfw.split`). Artifacts are saved one file per drug under ``models/``.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import joblib
import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression

from .config import MODELS_DIR


@dataclass
class DrugModel:
    drug: str
    tier: str
    features: List[str]                 # column order the model expects
    model: object                       # calibrated estimator (predict_proba)
    base_coef: Dict[str, float]         # gene -> L2 coefficient (explanation)
    n_train: int
    n_calib: int
    calibration_method: str
    synthetic_features: bool = False

    def vectorize(self, feat: Dict[str, int]) -> np.ndarray:
        return np.array([[int(feat.get(f, 0)) for f in self.features]], dtype=np.int8)

    def predict_proba(self, feat: Dict[str, int]) -> float:
        """Calibrated P(resistant) for a single genome's feature dict."""
        return float(self.model.predict_proba(self.vectorize(feat))[0, 1])

    def top_genes(self, feat: Dict[str, int], k: int = 4) -> List[Tuple[str, float]]:
        """Present features with the largest positive resistance coefficients."""
        present = [(g, self.base_coef.get(g, 0.0)) for g in self.features
                   if feat.get(g, 0) and self.base_coef.get(g, 0.0) > 0]
        present.sort(key=lambda t: -t[1])
        return present[:k]


def _fit_one(drug: str, tier: str, feats: List[str],
             Xtr: np.ndarray, ytr: np.ndarray,
             Xca: np.ndarray, yca: np.ndarray,
             synthetic: bool) -> Optional[DrugModel]:
    if len(np.unique(ytr)) < 2:
        return None
    base = LogisticRegression(penalty="l2", C=1.0, class_weight="balanced",
                              max_iter=2000, solver="liblinear")
    base.fit(Xtr, ytr)

    # Calibrate on the dedicated calibration split (prefit base estimator).
    min_class = int(min((yca == 0).sum(), (yca == 1).sum())) if len(yca) else 0
    if min_class >= 25 and len(np.unique(yca)) == 2:
        method = "isotonic"
    elif len(np.unique(yca)) == 2 and min_class >= 5:
        method = "sigmoid"
    else:
        method = "none"
    if method == "none":
        model = base  # fall back to uncalibrated base probabilities
    else:
        try:                                     # sklearn ≥ 1.6 API
            from sklearn.frozen import FrozenEstimator
            cal = CalibratedClassifierCV(FrozenEstimator(base), method=method)
        except Exception:                        # older sklearn
            cal = CalibratedClassifierCV(base, method=method, cv="prefit")
        cal.fit(Xca, yca)
        model = cal

    coef = dict(zip(feats, base.coef_.ravel().tolist()))
    return DrugModel(drug=drug, tier=tier, features=feats, model=model,
                     base_coef=coef, n_train=len(ytr), n_calib=len(yca),
                     calibration_method=method, synthetic_features=synthetic)


def train_models(features: pd.DataFrame,
                 labels: pd.DataFrame,
                 split: pd.DataFrame,
                 panel,
                 synthetic: bool = False) -> Dict[str, DrugModel]:
    """Fit one calibrated model per modelable (Tier A/B) drug."""
    part = split.set_index("genome_id")["partition"]
    feat_cols = list(features.columns)
    models: Dict[str, DrugModel] = {}
    for entry in panel:
        if not entry.modelable:
            continue
        sub = labels[labels["drug"] == entry.drug][["genome_id", "resistant"]]
        sub = sub[sub["genome_id"].isin(features.index)]
        sub = sub.assign(partition=sub["genome_id"].map(part))
        tr = sub[sub["partition"] == "train"]
        ca = sub[sub["partition"] == "calibration"]
        if tr.empty:
            continue
        Xtr = features.loc[tr["genome_id"]].to_numpy()
        ytr = tr["resistant"].to_numpy()
        Xca = features.loc[ca["genome_id"]].to_numpy() if not ca.empty else np.empty((0, len(feat_cols)))
        yca = ca["resistant"].to_numpy() if not ca.empty else np.empty((0,))
        m = _fit_one(entry.drug, entry.tier, feat_cols, Xtr, ytr, Xca, yca, synthetic)
        if m is not None:
            models[entry.drug] = m
    return models


# --------------------------------------------------------------------------- #
# Persistence.
# --------------------------------------------------------------------------- #
def save_models(models: Dict[str, DrugModel], models_dir: str = MODELS_DIR) -> None:
    os.makedirs(models_dir, exist_ok=True)
    manifest = {}
    for drug, m in models.items():
        joblib.dump(m, os.path.join(models_dir, f"{drug}.joblib"))
        manifest[drug] = {"tier": m.tier, "n_train": m.n_train, "n_calib": m.n_calib,
                          "calibration": m.calibration_method,
                          "synthetic_features": m.synthetic_features}
    with open(os.path.join(models_dir, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)


def load_models(models_dir: str = MODELS_DIR) -> Dict[str, DrugModel]:
    manifest_path = os.path.join(models_dir, "manifest.json")
    if not os.path.exists(manifest_path):
        return {}
    with open(manifest_path) as fh:
        manifest = json.load(fh)
    out: Dict[str, DrugModel] = {}
    for drug in manifest:
        p = os.path.join(models_dir, f"{drug}.joblib")
        if os.path.exists(p):
            out[drug] = joblib.load(p)
    return out
