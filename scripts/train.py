#!/usr/bin/env python3
"""Train the M3 predictor: one calibrated logistic regression per modelable drug.

Materialises the reproducible artifacts (panel, split) then fits + calibrates a
model per Tier-A/B drug on the grouped split and saves them under ``models/``.
Consumes whatever ``data/artifacts/features.parquet`` currently holds — the
teammates' real AMRFinderPlus matrix, or the synthetic placeholder.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gfw.io.labels import load_lab_labels
from gfw.m1_adapter import load_features
from gfw.panel import save_panel
from gfw.predict import save_models, train_models
from gfw.split import save_split, summarize


def main() -> int:
    panel = save_panel()
    split = save_split()
    print(summarize(split), "\n")

    features, synthetic = load_features()
    if synthetic:
        print("⚠️  Using SYNTHETIC placeholder M1 features — metrics are illustrative "
              "only until real AMRFinderPlus output replaces data/artifacts/features.parquet\n")
    labels = load_lab_labels()

    models = train_models(features, labels, split, panel, synthetic=synthetic)
    save_models(models)

    print(f"Trained {len(models)} per-drug models:")
    print(f"  {'drug':32s} {'tier':4s} {'n_train':>8s} {'n_calib':>8s} {'calib':>9s}")
    for drug, m in models.items():
        print(f"  {drug:32s} {m.tier:4s} {m.n_train:8d} {m.n_calib:8d} {m.calibration_method:>9s}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
