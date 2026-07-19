"""Evaluation on the grouped hidden test set (the brief's success criteria).

Per modelable drug, on the ``test`` partition only (clonal groups unseen in
training/calibration):

* balanced accuracy, recall for resistant and susceptible **reported separately**,
* F1, AUROC, PR-AUC (PR-AUC matters under class imbalance),
* **Brier score** + a reliability curve for confidence quality,
* **no-call rate** from the decision bands, and accuracy on the *called* subset.

Nothing here fits parameters — it only scores saved models on held-out groups.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from sklearn.metrics import (average_precision_score, balanced_accuracy_score,
                             brier_score_loss, f1_score, recall_score, roc_auc_score)

from .config import DECISION, DecisionConfig
from .predict import DrugModel


def _nan(x):
    return None if x is None or (isinstance(x, float) and np.isnan(x)) else float(x)


@dataclass
class DrugMetrics:
    drug: str
    tier: str
    n_test: int
    n_resistant: int
    n_susceptible: int
    auroc: Optional[float] = None
    pr_auc: Optional[float] = None
    balanced_accuracy: Optional[float] = None
    recall_resistant: Optional[float] = None
    recall_susceptible: Optional[float] = None
    f1: Optional[float] = None
    brier: Optional[float] = None
    no_call_rate: Optional[float] = None
    accuracy_on_called: Optional[float] = None
    reliability: Dict[str, List[float]] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, object]:
        d = self.__dict__.copy()
        return d


def _reliability(y: np.ndarray, p: np.ndarray, n_bins: int = 10) -> Dict[str, List[float]]:
    bins = np.linspace(0, 1, n_bins + 1)
    idx = np.clip(np.digitize(p, bins) - 1, 0, n_bins - 1)
    mean_pred, obs_freq, count = [], [], []
    for b in range(n_bins):
        m = idx == b
        if m.sum() == 0:
            continue
        mean_pred.append(float(p[m].mean()))
        obs_freq.append(float(y[m].mean()))
        count.append(int(m.sum()))
    return {"mean_predicted": mean_pred, "observed_frequency": obs_freq, "count": count}


def evaluate_drug(model: DrugModel, drug: str, tier: str,
                  Xte: np.ndarray, yte: np.ndarray,
                  cfg: DecisionConfig = DECISION) -> DrugMetrics:
    p = model.model.predict_proba(Xte)[:, 1]
    n_r, n_s = int((yte == 1).sum()), int((yte == 0).sum())
    m = DrugMetrics(drug=drug, tier=tier, n_test=len(yte),
                    n_resistant=n_r, n_susceptible=n_s)
    both = n_r > 0 and n_s > 0
    yhat = (p >= 0.5).astype(int)
    if both:
        m.auroc = _nan(roc_auc_score(yte, p))
        m.pr_auc = _nan(average_precision_score(yte, p))
        m.balanced_accuracy = _nan(balanced_accuracy_score(yte, yhat))
        m.recall_resistant = _nan(recall_score(yte, yhat, pos_label=1, zero_division=0))
        m.recall_susceptible = _nan(recall_score(yte, yhat, pos_label=0, zero_division=0))
        m.f1 = _nan(f1_score(yte, yhat, pos_label=1, zero_division=0))
    m.brier = _nan(brier_score_loss(yte, p)) if len(yte) else None
    m.reliability = _reliability(yte, p)

    # Decision bands -> no-call rate + accuracy on the called subset.
    hi, lo = cfg.p_fail_hi, cfg.p_work_lo
    if tier == "B":
        hi = min(0.9, hi + cfg.tier_b_nocall_widen)
        lo = max(0.1, lo - cfg.tier_b_nocall_widen)
    called = (p >= hi) | (p <= lo)
    m.no_call_rate = _nan(1.0 - called.mean()) if len(p) else None
    if called.sum() > 0:
        pred_called = (p[called] >= hi).astype(int)
        m.accuracy_on_called = _nan((pred_called == yte[called]).mean())
    return m


def evaluate_all(features: pd.DataFrame, labels: pd.DataFrame, split: pd.DataFrame,
                 models: Dict[str, DrugModel]) -> List[DrugMetrics]:
    part = split.set_index("genome_id")["partition"]
    out: List[DrugMetrics] = []
    for drug, model in models.items():
        sub = labels[labels["drug"] == drug][["genome_id", "resistant"]]
        sub = sub[sub["genome_id"].isin(features.index)]
        sub = sub.assign(partition=sub["genome_id"].map(part))
        te = sub[sub["partition"] == "test"]
        if te.empty:
            continue
        Xte = features.loc[te["genome_id"], model.features].to_numpy()
        yte = te["resistant"].to_numpy()
        out.append(evaluate_drug(model, drug, model.tier, Xte, yte))
    return out
