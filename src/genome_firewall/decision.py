from __future__ import annotations

import tomllib
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd

from genome_firewall.modeling.baseline import maximum_binary_jaccard

DEFAULT_DRUG_REGISTRY = Path("configs/drug_registry.toml")


def load_drug_registry(path: Path = DEFAULT_DRUG_REGISTRY) -> dict[str, Any]:
    with path.open("rb") as handle:
        registry = tomllib.load(handle)
    if "registry" not in registry or "drugs" not in registry:
        raise ValueError("Drug registry must contain [registry] and [drugs]")
    return registry


def build_feature_vector(
    evidence: pd.DataFrame, feature_columns: list[str]
) -> tuple[np.ndarray, list[str]]:
    detected = set(evidence.get("feature_key", pd.Series(dtype=object)).astype(str))
    known = set(feature_columns)
    vector = np.array([1 if feature in detected else 0 for feature in feature_columns], dtype=float)
    return vector, sorted(detected.difference(known))


def relevant_resistance_evidence(
    evidence: pd.DataFrame, resistance_terms: list[str]
) -> pd.DataFrame:
    if evidence.empty:
        return evidence.copy()
    terms = {term.upper() for term in resistance_terms}

    def matches(value: object) -> bool:
        tokens = {token.strip().upper() for token in str(value).split("/")}
        return bool(tokens.intersection(terms))

    return evidence.loc[evidence["amr_subclass"].map(matches)].copy()


def _strong_evidence(evidence: pd.DataFrame) -> pd.DataFrame:
    if evidence.empty:
        return evidence.copy()
    coverage = pd.to_numeric(evidence["coverage"], errors="coerce").fillna(0)
    identity = pd.to_numeric(evidence["identity"], errors="coerce").fillna(0)
    point = evidence["evidence_category"].eq("known_resistance_mutation")
    return evidence.loc[point | ((coverage >= 90) & (identity >= 90))].copy()


def _model_signal(probability: float, thresholds: dict[str, Any] | None) -> str:
    if not thresholds:
        return "no_call"
    lower = thresholds.get("lower")
    upper = thresholds.get("upper")
    if lower is not None and probability <= float(lower):
        return "susceptible_signal"
    if upper is not None and probability >= float(upper):
        return "resistant_signal"
    return "no_call"


def predict_antibiotic(
    *,
    antibiotic: str,
    evidence: pd.DataFrame,
    model_path: Path | None = None,
    model_artifact: dict[str, Any] | None = None,
    drug: dict[str, Any],
    target_status: dict[str, Any],
    lineage_status: dict[str, Any],
    qc_passed: bool,
    qc_reasons: list[str],
) -> dict[str, Any]:
    if model_artifact is None:
        if model_path is None:
            raise ValueError("model_path or model_artifact is required")
        artifact = joblib.load(model_path)
    else:
        artifact = model_artifact
    features, unknown_features = build_feature_vector(evidence, artifact["feature_columns"])
    feature_frame = pd.DataFrame([features], columns=artifact["feature_columns"])
    probability = float(artifact["estimator"].predict_proba(feature_frame)[0, 1])
    thresholds = artifact.get("decision_thresholds")
    signal = _model_signal(probability, thresholds)
    relevant = relevant_resistance_evidence(evidence, drug["resistance_terms"])
    strong = _strong_evidence(relevant)

    reference = artifact.get("novelty_reference", np.empty((0, len(features))))
    similarity = maximum_binary_jaccard(features, reference)
    novelty_floor = float(artifact.get("novelty_min_jaccard", 0.0))
    in_distribution = not unknown_features and similarity >= novelty_floor

    reasons: list[str] = []
    call = "no_call"
    evidence_category = "no_known_resistance_signal"
    if not qc_passed:
        reasons.extend(f"qc:{reason}" for reason in qc_reasons)
    elif lineage_status["status"] != "in_distribution":
        reasons.append(f"lineage_{lineage_status['status']}")
    elif not in_distribution:
        reasons.append("outside_amr_feature_distribution")
        if unknown_features:
            reasons.append("unseen_amr_features_detected")
    elif signal == "susceptible_signal" and not relevant.empty:
        reasons.append("model_susceptible_but_resistance_evidence_detected")
        evidence_category = "conflicting_known_resistance_evidence"
    elif signal == "susceptible_signal" and target_status["status"] != "present":
        reasons.append("molecular_target_not_verified")
    elif signal == "susceptible_signal":
        call = "likely_to_work"
        evidence_category = "no_known_resistance_signal"
    elif signal == "resistant_signal" and not strong.empty:
        call = "likely_to_fail"
        evidence_category = "known_resistance_gene_or_mutation"
    elif signal == "resistant_signal":
        call = "likely_to_fail"
        evidence_category = "statistical_association_only"
    else:
        reasons.append("calibrated_probability_inside_no_call_region")
        if not relevant.empty:
            evidence_category = "known_resistance_evidence_but_model_uncertain"

    confidence = 1 - probability if call == "likely_to_work" else probability
    if call == "no_call":
        confidence = max(probability, 1 - probability)
    return {
        "antibiotic": antibiotic,
        "call": call,
        "confidence": confidence,
        "probability_resistant": probability,
        "model_signal": signal,
        "evidence_category": evidence_category,
        "supporting_elements": strong["element_symbol"].astype(str).tolist(),
        "all_relevant_elements": relevant["element_symbol"].astype(str).tolist(),
        "target_status": target_status["status"],
        "target_label": drug["target_label"],
        "targets_detected": target_status.get("detected", []),
        "feature_similarity": similarity,
        "feature_similarity_floor": novelty_floor,
        "unknown_features": unknown_features,
        "lineage_status": lineage_status["status"],
        "maximum_training_ani": lineage_status.get("maximum_training_ani"),
        "minimum_training_ani": lineage_status.get("minimum_training_ani"),
        "nearest_training_genome": lineage_status.get("nearest_training_genome"),
        "reasons": reasons,
    }
