#!/usr/bin/env python3
"""Score saved models on the grouped hidden test set; write metrics + plots.

Outputs (under ``reports/``):
  * metrics.json / metrics.csv  — per-drug metrics table (consumed by the UI)
  * reliability.png             — calibration curves for the Tier-A drugs
  * performance.png             — AUROC / PR-AUC / balanced-accuracy bars
"""

from __future__ import annotations

import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gfw.config import REPORTS_DIR, ensure_dirs
from gfw.evaluate import evaluate_all
from gfw.io.labels import load_lab_labels
from gfw.m1_adapter import load_features
from gfw.predict import load_models
from gfw.split import load_split, make_split


def _plot_reliability(metrics, path, synthetic):
    tier_a = [m for m in metrics if m.tier == "A" and m.reliability.get("mean_predicted")]
    if not tier_a:
        return
    n = len(tier_a)
    cols = 4
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(3.2 * cols, 3.0 * rows), squeeze=False)
    for ax, m in zip(axes.ravel(), tier_a):
        ax.plot([0, 1], [0, 1], "--", color="#9aa5b1", lw=1)
        mp = m.reliability["mean_predicted"]; of = m.reliability["observed_frequency"]
        ax.plot(mp, of, "o-", color="#2b6cb0", lw=1.6, ms=4)
        ax.set_title(f"{m.drug}\nAUROC={m.auroc:.2f} Brier={m.brier:.2f}", fontsize=8)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xlabel("predicted p(R)", fontsize=7); ax.set_ylabel("observed", fontsize=7)
        ax.tick_params(labelsize=6)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    tag = "  [SYNTHETIC PLACEHOLDER — illustrative only]" if synthetic else ""
    fig.suptitle("Reliability (grouped hidden test)" + tag, fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(path, dpi=130)
    plt.close(fig)


def _plot_performance(df, path, synthetic):
    d = df[df["auroc"].notna()].sort_values("auroc", ascending=True)
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=(8, max(3, 0.4 * len(d))))
    y = np.arange(len(d))
    ax.barh(y - 0.2, d["auroc"], height=0.2, label="AUROC", color="#2b6cb0")
    ax.barh(y, d["pr_auc"], height=0.2, label="PR-AUC", color="#38a169")
    ax.barh(y + 0.2, d["balanced_accuracy"], height=0.2, label="Balanced acc.", color="#dd6b20")
    ax.set_yticks(y); ax.set_yticklabels(d["drug"], fontsize=8)
    ax.set_xlim(0, 1); ax.legend(fontsize=8, loc="lower right")
    tag = "  [SYNTHETIC PLACEHOLDER]" if synthetic else ""
    ax.set_title("Per-drug performance on grouped hidden test" + tag, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main() -> int:
    ensure_dirs()
    features, synthetic = load_features()
    labels = load_lab_labels()
    try:
        split = load_split()
    except Exception:
        split = make_split()
    models = load_models()
    if not models:
        print("No models found — run scripts/train.py first.")
        return 1

    metrics = evaluate_all(features, labels, split, models)
    rows = [m.as_dict() for m in metrics]
    for r in rows:
        r.pop("reliability", None)
    df = pd.DataFrame(rows)

    with open(os.path.join(REPORTS_DIR, "metrics.json"), "w") as fh:
        json.dump({"synthetic_features": bool(synthetic),
                   "metrics": [m.as_dict() for m in metrics]}, fh, indent=2)
    df.to_csv(os.path.join(REPORTS_DIR, "metrics.csv"), index=False)
    _plot_reliability(metrics, os.path.join(REPORTS_DIR, "reliability.png"), synthetic)
    _plot_performance(df, os.path.join(REPORTS_DIR, "performance.png"), synthetic)

    if synthetic:
        print("⚠️  SYNTHETIC placeholder features — metrics below are ILLUSTRATIVE ONLY.\n")
    show = ["drug", "tier", "n_test", "auroc", "pr_auc", "balanced_accuracy",
            "recall_resistant", "recall_susceptible", "brier", "no_call_rate"]
    with pd.option_context("display.width", 200, "display.max_columns", None):
        print(df[show].round(3).to_string(index=False))
    print(f"\nWrote metrics + plots to {REPORTS_DIR}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
