from __future__ import annotations

import json
import hashlib
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.frozen import FrozenEstimator
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    f1_score,
    recall_score,
    roc_auc_score,
)

LABEL_MAP = {"Susceptible": 0, "Resistant": 1}


def feature_schema_sha256(feature_columns: list[str]) -> str:
    return hashlib.sha256("\n".join(feature_columns).encode("utf-8")).hexdigest()


def maximum_binary_jaccard(vector: np.ndarray, reference: np.ndarray) -> float:
    """Return maximum Jaccard similarity to binary training feature profiles."""
    if reference.size == 0:
        return 0.0
    vector_bool = vector.astype(bool)
    reference_bool = reference.astype(bool)
    intersection = np.logical_and(reference_bool, vector_bool).sum(axis=1)
    union = np.logical_or(reference_bool, vector_bool).sum(axis=1)
    similarity = np.divide(
        intersection,
        union,
        out=np.ones_like(intersection, dtype=float),
        where=union != 0,
    )
    return float(similarity.max())


def fit_feature_novelty_reference(
    values: np.ndarray, *, quantile: float
) -> tuple[np.ndarray, float]:
    """Store unique training profiles and a leave-one-profile-out similarity floor."""
    reference = np.unique(values.astype("uint8"), axis=0)
    if len(reference) <= 1:
        return reference, 1.0
    similarities: list[float] = []
    for index, vector in enumerate(reference):
        others = np.delete(reference, index, axis=0)
        similarities.append(maximum_binary_jaccard(vector, others))
    return reference, float(np.quantile(similarities, quantile))


def learn_no_call_thresholds(
    target: np.ndarray,
    probability_resistant: np.ndarray,
    *,
    max_resistant_fraction_in_susceptible_calls: float,
    max_susceptible_fraction_in_resistant_calls: float,
    minimum_calls: int,
) -> dict[str, float | int | None] | None:
    """Learn conservative call boundaries from the held-out calibration partition."""
    candidates = np.unique(probability_resistant)
    valid_lowers: list[tuple[float, int, float]] = []
    for threshold in candidates:
        called = probability_resistant <= threshold
        count = int(called.sum())
        error = float(target[called].mean()) if count else 1.0
        if count >= minimum_calls and error <= max_resistant_fraction_in_susceptible_calls:
            valid_lowers.append((float(threshold), count, error))

    valid_uppers: list[tuple[float, int, float]] = []
    for threshold in candidates:
        called = probability_resistant >= threshold
        count = int(called.sum())
        error = float((target[called] == 0).mean()) if count else 1.0
        if count >= minimum_calls and error <= max_susceptible_fraction_in_resistant_calls:
            valid_uppers.append((float(threshold), count, error))

    if not valid_lowers and not valid_uppers:
        return None

    # Optimize the two call regions jointly. Independently taking the widest valid
    # prefix and suffix can make them overlap on noisy calibration predictions.
    valid_pairs = [
        (lower, upper)
        for lower in valid_lowers
        for upper in valid_uppers
        if lower[0] < upper[0]
    ]
    if valid_pairs:
        lower, upper = max(
            valid_pairs,
            key=lambda pair: (
                pair[0][1] + pair[1][1],
                min(pair[0][1], pair[1][1]),
                pair[0][1],
            ),
        )
    else:
        best_lower = max(valid_lowers, key=lambda item: item[1]) if valid_lowers else None
        best_upper = max(valid_uppers, key=lambda item: item[1]) if valid_uppers else None
        if best_lower is not None and (
            best_upper is None or best_lower[1] >= best_upper[1]
        ):
            lower, upper = best_lower, None
        else:
            lower, upper = None, best_upper

    return {
        "lower": lower[0] if lower is not None else None,
        "upper": upper[0] if upper is not None else None,
        "susceptible_calls": lower[1] if lower is not None else 0,
        "susceptible_call_error": lower[2] if lower is not None else None,
        "resistant_calls": upper[1] if upper is not None else 0,
        "resistant_call_error": upper[2] if upper is not None else None,
    }


def load_modeling_table(
    features_path: Path,
    phenotypes_path: Path,
    splits_path: Path,
) -> tuple[pd.DataFrame, list[str]]:
    features = pd.read_parquet(features_path)
    phenotypes = pd.read_csv(phenotypes_path, dtype=object, keep_default_na=False)
    splits = pd.read_csv(splits_path, dtype=object, keep_default_na=False)
    feature_columns = sorted(column for column in features.columns if column != "genome_id")
    table = phenotypes.merge(features, on="genome_id", how="inner", validate="many_to_one")
    table = table.merge(
        splits[["genome_id", "cluster_id", "split"]],
        on="genome_id",
        how="inner",
        validate="many_to_one",
    )
    table["target"] = table["label"].map(LABEL_MAP)
    if table["target"].isna().any():
        raise ValueError("Unexpected non-binary phenotype reached model training")
    table[feature_columns] = table[feature_columns].fillna(0).astype("float64")
    return table, feature_columns


def _class_counts(values: pd.Series) -> dict[str, int]:
    counts = values.value_counts()
    return {
        "susceptible": int(counts.get(0, 0)),
        "resistant": int(counts.get(1, 0)),
    }


def _evaluate(y_true: np.ndarray, probability: np.ndarray) -> dict[str, float | None]:
    if len(y_true) == 0:
        return {}
    predicted = (probability >= 0.5).astype(int)
    both_classes = len(np.unique(y_true)) == 2
    return {
        "balanced_accuracy": float(balanced_accuracy_score(y_true, predicted))
        if both_classes
        else None,
        "resistant_recall": float(recall_score(y_true, predicted, pos_label=1, zero_division=0)),
        "susceptible_recall": float(recall_score(y_true, predicted, pos_label=0, zero_division=0)),
        "f1": float(f1_score(y_true, predicted, zero_division=0)),
        "auroc": float(roc_auc_score(y_true, probability)) if both_classes else None,
        "pr_auc": float(average_precision_score(y_true, probability)) if both_classes else None,
        "brier": float(brier_score_loss(y_true, probability)),
    }


def _called_metrics(predictions: pd.DataFrame) -> dict[str, float | int | None]:
    if predictions.empty:
        return {"called_count": 0, "coverage": None, "accuracy": None}
    called = predictions.loc[predictions["model_signal"].ne("no_call")].copy()
    metrics: dict[str, float | int | None] = {
        "called_count": len(called),
        "coverage": float(len(called) / len(predictions)),
        "accuracy": None,
        "susceptible_calls": int(called["model_signal"].eq("susceptible_signal").sum()),
        "resistant_calls": int(called["model_signal"].eq("resistant_signal").sum()),
        "resistant_fraction_in_susceptible_calls": None,
        "susceptible_fraction_in_resistant_calls": None,
    }
    if called.empty:
        return metrics
    predicted = called["model_signal"].map(
        {"susceptible_signal": 0, "resistant_signal": 1}
    )
    metrics["accuracy"] = float(predicted.eq(called["target"]).mean())
    susceptible_calls = called.loc[called["model_signal"].eq("susceptible_signal")]
    resistant_calls = called.loc[called["model_signal"].eq("resistant_signal")]
    if len(susceptible_calls):
        metrics["resistant_fraction_in_susceptible_calls"] = float(
            susceptible_calls["target"].mean()
        )
    if len(resistant_calls):
        metrics["susceptible_fraction_in_resistant_calls"] = float(
            resistant_calls["target"].eq(0).mean()
        )
    return metrics


def reliability_table(predictions: pd.DataFrame, *, bins: int = 10) -> pd.DataFrame:
    """Build plot-ready reliability points from untouched test predictions."""
    columns = [
        "antibiotic",
        "probability_bin",
        "bin_lower",
        "bin_upper",
        "samples",
        "mean_probability_resistant",
        "observed_resistant_fraction",
    ]
    if predictions.empty:
        return pd.DataFrame(columns=columns)
    frame = predictions.copy()
    edges = np.linspace(0.0, 1.0, bins + 1)
    frame["probability_bin"] = pd.cut(
        frame["probability_resistant"], edges, labels=False, include_lowest=True
    )
    frame = frame.dropna(subset=["probability_bin"])
    frame["probability_bin"] = frame["probability_bin"].astype(int)
    result = (
        frame.groupby(["antibiotic", "probability_bin"])
        .agg(
            samples=("target", "size"),
            mean_probability_resistant=("probability_resistant", "mean"),
            observed_resistant_fraction=("target", "mean"),
        )
        .reset_index()
    )
    result["bin_lower"] = result["probability_bin"].map(lambda value: edges[value])
    result["bin_upper"] = result["probability_bin"].map(lambda value: edges[value + 1])
    return result[columns]


def train_drug_model(
    table: pd.DataFrame,
    feature_columns: list[str],
    *,
    antibiotic: str,
    config: dict[str, Any],
    output_directory: Path,
) -> tuple[dict[str, Any], pd.DataFrame]:
    rows = table.loc[table["antibiotic"].eq(antibiotic)].copy()
    train = rows.loc[rows["split"].eq("train")]
    calibration = rows.loc[rows["split"].eq("calibration")]
    test = rows.loc[rows["split"].eq("test")]
    counts = {
        "train": _class_counts(train["target"]),
        "calibration": _class_counts(calibration["target"]),
        "test": _class_counts(test["target"]),
    }
    metadata: dict[str, Any] = {
        "antibiotic": antibiotic,
        "status": "skipped",
        "calibration_status": "unavailable",
        "class_counts": counts,
        "feature_count": len(feature_columns),
    }
    minimum_train = config["min_train_per_class"]
    if min(counts["train"].values()) < minimum_train:
        metadata["reason"] = f"training requires at least {minimum_train} samples per class"
        return metadata, pd.DataFrame()

    base = LogisticRegression(
        C=config["regularization_c"],
        class_weight="balanced",
        solver="liblinear",
        max_iter=1000,
        random_state=config["seed"],
    )
    base.fit(train[feature_columns], train["target"].astype(int))
    novelty_reference, novelty_min_jaccard = fit_feature_novelty_reference(
        train[feature_columns].to_numpy(dtype="uint8"),
        quantile=config.get("novelty_quantile", 0.05),
    )
    estimator: Any = base
    minimum_calibration = config["min_calibration_per_class"]
    if min(counts["calibration"].values()) >= minimum_calibration:
        estimator = CalibratedClassifierCV(
            estimator=FrozenEstimator(base), method="sigmoid"
        )
        estimator.fit(
            calibration[feature_columns], calibration["target"].astype(int)
        )
        metadata["calibration_status"] = "sigmoid"

    thresholds = None
    if metadata["calibration_status"] == "sigmoid":
        calibration_probability = estimator.predict_proba(calibration[feature_columns])[:, 1]
        thresholds = learn_no_call_thresholds(
            calibration["target"].to_numpy(dtype=int),
            calibration_probability,
            max_resistant_fraction_in_susceptible_calls=config.get(
                "max_resistant_fraction_in_susceptible_calls", 0.10
            ),
            max_susceptible_fraction_in_resistant_calls=config.get(
                "max_susceptible_fraction_in_resistant_calls", 0.10
            ),
            minimum_calls=config.get("min_threshold_calls", 5),
        )
        if thresholds is None:
            metadata["calibration_status"] = "thresholds_unavailable"
        elif thresholds["lower"] is None or thresholds["upper"] is None:
            metadata["calibration_status"] = "sigmoid_partial_thresholds"
    metadata["decision_thresholds"] = thresholds

    output_directory.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "schema_version": "genome-firewall-model-v1",
            "antibiotic": antibiotic,
            "estimator": estimator,
            "base_estimator": base,
            "feature_columns": feature_columns,
            "calibration_status": metadata["calibration_status"],
            "decision_thresholds": thresholds,
            "novelty_reference": novelty_reference,
            "novelty_min_jaccard": novelty_min_jaccard,
            "feature_schema_sha256": feature_schema_sha256(feature_columns),
        },
        output_directory / f"{antibiotic.replace('/', '_')}.joblib",
    )

    predictions = test[["genome_id", "cluster_id", "label", "target"]].copy()
    if len(test):
        predictions["raw_probability_resistant"] = base.predict_proba(
            test[feature_columns]
        )[:, 1]
        predictions["probability_resistant"] = estimator.predict_proba(
            test[feature_columns]
        )[:, 1]
    else:
        predictions["raw_probability_resistant"] = pd.Series(dtype=float)
        predictions["probability_resistant"] = pd.Series(dtype=float)

    if metadata["calibration_status"].startswith("sigmoid") and thresholds is not None:
        probability = predictions["probability_resistant"]
        conditions = []
        choices = []
        if thresholds["lower"] is not None:
            conditions.append(probability <= thresholds["lower"])
            choices.append("susceptible_signal")
        if thresholds["upper"] is not None:
            conditions.append(probability >= thresholds["upper"])
            choices.append("resistant_signal")
        predictions["model_signal"] = np.select(
            conditions,
            choices,
            default="no_call",
        )
    else:
        predictions["model_signal"] = "no_call"

    metadata["status"] = "trained"
    metadata["test_metrics"] = _evaluate(
        predictions["target"].to_numpy(dtype=int),
        predictions["probability_resistant"].to_numpy(dtype=float),
    )
    metadata["no_call_rate"] = (
        float(predictions["model_signal"].eq("no_call").mean()) if len(predictions) else None
    )
    metadata["called_test_metrics"] = _called_metrics(predictions)
    return metadata, predictions


def train_all_drugs(
    table: pd.DataFrame,
    feature_columns: list[str],
    *,
    antibiotics: list[str],
    config: dict[str, Any],
    output_directory: Path,
    bundle_metadata: dict[str, Any] | None = None,
) -> pd.DataFrame:
    output_directory.mkdir(parents=True, exist_ok=True)
    summaries: list[dict[str, Any]] = []
    prediction_frames: list[pd.DataFrame] = []
    for antibiotic in antibiotics:
        summary, predictions = train_drug_model(
            table,
            feature_columns,
            antibiotic=antibiotic,
            config=config,
            output_directory=output_directory,
        )
        summaries.append(summary)
        if not predictions.empty:
            predictions.insert(1, "antibiotic", antibiotic)
            prediction_frames.append(predictions)
        (output_directory / f"{antibiotic.replace('/', '_')}.metadata.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    predictions = (
        pd.concat(prediction_frames, ignore_index=True)
        if prediction_frames
        else pd.DataFrame()
    )
    predictions.to_csv(output_directory / "test-predictions.csv", index=False)
    reliability_table(predictions).to_csv(
        output_directory / "test-reliability.csv", index=False
    )
    summary_frame = pd.json_normalize(summaries)
    summary_frame.to_csv(output_directory / "model-summary.csv", index=False)
    trained = [summary["antibiotic"] for summary in summaries if summary["status"] == "trained"]
    bundle = {
        "schema_version": "genome-firewall-bundle-v1",
        "antibiotics": trained,
        "feature_count": len(feature_columns),
        "feature_columns": feature_columns,
        "feature_schema_sha256": feature_schema_sha256(feature_columns),
        "models": {
            antibiotic: {
                "path": f"{antibiotic.replace('/', '_')}.joblib",
                "metadata_path": f"{antibiotic.replace('/', '_')}.metadata.json",
                "evaluation": {
                    key: next(
                        summary[key]
                        for summary in summaries
                        if summary["antibiotic"] == antibiotic
                    )
                    for key in [
                        "calibration_status",
                        "class_counts",
                        "test_metrics",
                        "no_call_rate",
                        "called_test_metrics",
                    ]
                },
            }
            for antibiotic in trained
        },
        **(bundle_metadata or {}),
    }
    (output_directory / "bundle-manifest.json").write_text(
        json.dumps(bundle, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary_frame
