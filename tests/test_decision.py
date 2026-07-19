from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from genome_firewall.decision import (
    build_feature_vector,
    predict_antibiotic,
    relevant_resistance_evidence,
)


def test_feature_vector_reports_unseen_amr_elements() -> None:
    evidence = pd.DataFrame({"feature_key": ["gene::known", "gene::new"]})
    vector, unknown = build_feature_vector(evidence, ["gene::known"])
    assert vector.tolist() == [1]
    assert unknown == ["gene::new"]


def test_drug_specific_subclass_mapping() -> None:
    evidence = pd.DataFrame(
        {
            "amr_subclass": ["CLINDAMYCIN/ERYTHROMYCIN", "TETRACYCLINE"],
            "element_symbol": ["erm(C)", "tet(K)"],
        }
    )
    relevant = relevant_resistance_evidence(evidence, ["ERYTHROMYCIN"])
    assert relevant["element_symbol"].tolist() == ["erm(C)"]


def test_known_resistance_conflict_forces_no_call(tmp_path: Path) -> None:
    columns = ["gene::erm(C)"]
    estimator = LogisticRegression().fit(
        pd.DataFrame([[0], [0], [1], [1]], columns=columns), [0, 0, 1, 1]
    )
    model_path = tmp_path / "erythromycin.joblib"
    joblib.dump(
        {
            "estimator": estimator,
            "feature_columns": columns,
            "decision_thresholds": {"lower": 1.0, "upper": None},
            "novelty_reference": np.array([[0], [1]], dtype="uint8"),
            "novelty_min_jaccard": 0.0,
        },
        model_path,
    )
    evidence = pd.DataFrame(
        {
            "feature_key": ["gene::erm(C)"],
            "element_symbol": ["erm(C)"],
            "evidence_category": ["known_resistance_gene"],
            "amr_subclass": ["CLINDAMYCIN/ERYTHROMYCIN"],
            "coverage": [100.0],
            "identity": [100.0],
        }
    )
    decision = predict_antibiotic(
        antibiotic="erythromycin",
        evidence=evidence,
        model_path=model_path,
        drug={
            "resistance_terms": ["ERYTHROMYCIN"],
            "target_label": "23S rRNA",
        },
        target_status={"status": "present", "detected": ["23S rRNA"]},
        lineage_status={
            "status": "in_distribution",
            "maximum_training_ani": 1.0,
            "minimum_training_ani": 0.95,
            "nearest_training_genome": "g1",
        },
        qc_passed=True,
        qc_reasons=[],
    )
    # The synthetic threshold makes this a susceptible model signal, but known
    # drug-specific resistance evidence prevents a likely-to-work call.
    assert decision["call"] == "no_call"
    assert "model_susceptible_but_resistance_evidence_detected" in decision["reasons"]


def test_target_gate_does_not_hide_resistant_signal(tmp_path: Path) -> None:
    columns = ["gene::mecA"]
    estimator = LogisticRegression().fit(
        pd.DataFrame([[0], [0], [1], [1]], columns=columns), [0, 0, 1, 1]
    )
    decision = predict_antibiotic(
        antibiotic="cefoxitin",
        evidence=pd.DataFrame(
            {
                "feature_key": ["gene::mecA"],
                "element_symbol": ["mecA"],
                "evidence_category": ["known_resistance_gene"],
                "amr_subclass": ["METHICILLIN"],
                "coverage": [100.0],
                "identity": [100.0],
            }
        ),
        model_artifact={
            "estimator": estimator,
            "feature_columns": columns,
            "decision_thresholds": {"lower": None, "upper": 0.0},
            "novelty_reference": np.array([[0], [1]], dtype="uint8"),
            "novelty_min_jaccard": 0.0,
        },
        drug={"resistance_terms": ["METHICILLIN"], "target_label": "PBPs"},
        target_status={"status": "not_verified", "detected": []},
        lineage_status={
            "status": "in_distribution",
            "maximum_training_ani": 1.0,
            "minimum_training_ani": 0.95,
            "nearest_training_genome": "g1",
        },
        qc_passed=True,
        qc_reasons=[],
    )
    assert decision["call"] == "likely_to_fail"
    assert decision["evidence_category"] == "known_resistance_gene_or_mutation"
