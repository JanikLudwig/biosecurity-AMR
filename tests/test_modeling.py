from pathlib import Path

import numpy as np
import pandas as pd

from genome_firewall.modeling.baseline import (
    fit_feature_novelty_reference,
    learn_no_call_thresholds,
    maximum_binary_jaccard,
    train_drug_model,
)


def test_binary_feature_novelty_reference() -> None:
    values = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]], dtype="uint8")
    reference, threshold = fit_feature_novelty_reference(values, quantile=0.0)
    assert len(reference) == 3
    assert 0 <= threshold <= 1
    assert maximum_binary_jaccard(np.array([1, 0, 0]), reference) == 1.0


def test_learns_separated_no_call_thresholds() -> None:
    thresholds = learn_no_call_thresholds(
        np.array([0, 0, 0, 1, 1, 1]),
        np.array([0.05, 0.10, 0.20, 0.75, 0.85, 0.95]),
        max_resistant_fraction_in_susceptible_calls=0.1,
        max_susceptible_fraction_in_resistant_calls=0.1,
        minimum_calls=2,
    )
    assert thresholds is not None
    assert thresholds["lower"] == 0.20
    assert thresholds["upper"] == 0.75


def test_keeps_one_safe_call_boundary() -> None:
    thresholds = learn_no_call_thresholds(
        np.array([0, 0, 0, 1, 0, 1]),
        np.array([0.05, 0.10, 0.20, 0.60, 0.80, 0.90]),
        max_resistant_fraction_in_susceptible_calls=0.0,
        max_susceptible_fraction_in_resistant_calls=0.0,
        minimum_calls=3,
    )
    assert thresholds is not None
    assert thresholds["lower"] == 0.20
    assert thresholds["upper"] is None


def test_uncalibrated_model_returns_only_no_calls(tmp_path: Path) -> None:
    rows = []
    for index in range(14):
        split = "train" if index < 10 else "test"
        rows.append(
            {
                "genome_id": f"g{index}",
                "cluster_id": f"c{index}",
                "split": split,
                "antibiotic": "drug",
                "target": index % 2,
                "label": "Resistant" if index % 2 else "Susceptible",
                "gene::x": index % 2,
            }
        )
    table = pd.DataFrame(rows)
    metadata, predictions = train_drug_model(
        table,
        ["gene::x"],
        antibiotic="drug",
        config={
            "regularization_c": 1.0,
            "seed": 42,
            "min_train_per_class": 5,
            "min_calibration_per_class": 2,
            "no_call_lower": 0.35,
            "no_call_upper": 0.65,
        },
        output_directory=tmp_path,
    )
    assert metadata["status"] == "trained"
    assert metadata["calibration_status"] == "unavailable"
    assert set(predictions["model_signal"]) == {"no_call"}
    assert (tmp_path / "drug.joblib").is_file()
