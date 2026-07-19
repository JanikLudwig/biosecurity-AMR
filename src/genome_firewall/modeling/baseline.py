from __future__ import annotations

import json
from pathlib import Path
from typing import Any
import logging

import joblib
import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.isotonic import IsotonicRegression
from sklearn.frozen import FrozenEstimator
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, GroupShuffleSplit
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    f1_score,
    recall_score,
    precision_score,
    roc_auc_score,
    confusion_matrix,
)
import math

class RobustSigmoidCalibrator:
    def __init__(self):
        self.calibrator = LogisticRegression(solver="liblinear", C=1.0)

    def fit(self, scores: np.ndarray, y: np.ndarray):
        X = scores.reshape(-1, 1)
        self.calibrator.fit(X, y)
        return self

    def predict_proba(self, scores: np.ndarray) -> np.ndarray:
        X = scores.reshape(-1, 1)
        return self.calibrator.predict_proba(X)


def build_training_dataset(
    features_path: Path | str,
    labels_path: Path | str,
    antibiotic: str,
) -> tuple[pd.DataFrame, list[str]]:
    """Build a combined dataset for a single antibiotic."""
    features = pd.read_csv(features_path)
    labels = pd.read_csv(labels_path)

    # Filter for the specific antibiotic
    labels = labels[labels["antibiotic"] == antibiotic].copy()

    # Filter only R and S labels
    valid_labels = labels[labels["label"].isin(["R", "S"])].copy()
    valid_labels["target"] = valid_labels["label"].map({"R": 1, "S": 0})

    # Join features and labels
    table = valid_labels.merge(
        features,
        on="genome_id",
        how="inner",
        validate="one_to_one"
    )

    # Exclude genome_id and label columns from feature list
    non_feature_cols = {"genome_id", "target", "label", "antibiotic", "Evidence", "evidence"}
    feature_columns = sorted([col for col in features.columns if col not in non_feature_cols])

    # Missing features -> 0
    table[feature_columns] = table[feature_columns].fillna(0).astype(float)

    return table, feature_columns

def assign_development_split(
    table: pd.DataFrame,
    train_fraction: float = 0.70,
    calibration_fraction: float = 0.15,
    test_fraction: float = 0.15,
    seed: int = 42,
    groups_path: Path | str | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Assign development split using either random stratify or grouped split."""
    table = table.sort_values("genome_id").reset_index(drop=True)
    table = table.copy()
    table["split"] = "none"

    metadata = {}
    metadata["split_seed"] = seed

    if groups_path is None:
        train_idx, temp_idx = train_test_split(
            table.index,
            test_size=(1.0 - train_fraction),
            random_state=seed,
            stratify=table["target"]
        )

        rel_calib_fraction = calibration_fraction / (calibration_fraction + test_fraction)
        calib_idx, test_idx = train_test_split(
            temp_idx,
            test_size=(1.0 - rel_calib_fraction),
            random_state=seed,
            stratify=table.loc[temp_idx, "target"]
        )

        table.loc[train_idx, "split"] = "train"
        table.loc[calib_idx, "split"] = "calibration"
        table.loc[test_idx, "split"] = "test"

        metadata["evaluation_status"] = "development-only"
        metadata["split_method"] = "development-only random split"
        metadata["groups_used"] = False
        return table, metadata
    else:
        if not Path(groups_path).exists():
            raise FileNotFoundError(f"Groups file not found: {groups_path}")

        groups_df = pd.read_csv(groups_path)
        if "genome_id" not in groups_df.columns or "group_id" not in groups_df.columns:
            raise ValueError("groups_path CSV must contain 'genome_id' and 'group_id' columns")

        if groups_df["genome_id"].isnull().any():
            raise ValueError("Missing genome_id values in groups file")
        if groups_df["group_id"].isnull().any() or (groups_df["group_id"] == "").any():
            raise ValueError("Missing or empty group_id values in groups file")

        if groups_df["genome_id"].duplicated().any():
            raise ValueError("Duplicate genome_id values in groups file")

        group_map = dict(zip(groups_df["genome_id"], groups_df["group_id"]))

        missing_genomes = set(table["genome_id"]) - set(group_map.keys())
        if missing_genomes:
            raise ValueError(f"Genomes missing from groups file: {len(missing_genomes)}")

        table["group_id"] = table["genome_id"].map(group_map)

        extra_genomes = set(group_map.keys()) - set(table["genome_id"])
        if extra_genomes:
            logging.info(f"Ignored {len(extra_genomes)} extra genomes from groups file not in current dataset.")

        best_split = None
        best_diff = float("inf")

        for attempt in range(100):
            current_seed = seed + attempt
            gss1 = GroupShuffleSplit(n_splits=1, test_size=test_fraction, random_state=current_seed)
            try:
                train_cal_idx, test_idx = next(gss1.split(table, table["target"], groups=table["group_id"]))
            except Exception:
                continue

            train_cal_table = table.iloc[train_cal_idx].reset_index(drop=True)

            rel_calib_fraction = calibration_fraction / (1.0 - test_fraction)
            gss2 = GroupShuffleSplit(n_splits=1, test_size=rel_calib_fraction, random_state=current_seed)
            try:
                tr_idx_sub, cal_idx_sub = next(gss2.split(train_cal_table, train_cal_table["target"], groups=train_cal_table["group_id"]))
            except Exception:
                continue

            train_idx = train_cal_idx[tr_idx_sub]
            calib_idx = train_cal_idx[cal_idx_sub]

            if (len(np.unique(table.iloc[train_idx]["target"])) == 2 and
                len(np.unique(table.iloc[calib_idx]["target"])) == 2 and
                len(np.unique(table.iloc[test_idx]["target"])) == 2):

                actual_train = len(train_idx) / len(table)
                diff = abs(actual_train - train_fraction)

                if diff < best_diff:
                    best_diff = diff
                    best_split = (train_idx, calib_idx, test_idx)
                    metadata["split_attempt"] = attempt + 1

        if best_split is None:
            raise ValueError("Could not find a valid group split with both classes in all splits after 100 attempts.")

        train_idx, calib_idx, test_idx = best_split

        train_g = set(table.iloc[train_idx]["group_id"])
        calib_g = set(table.iloc[calib_idx]["group_id"])
        test_g = set(table.iloc[test_idx]["group_id"])

        if not train_g.isdisjoint(calib_g) or not train_g.isdisjoint(test_g) or not calib_g.isdisjoint(test_g):
            raise ValueError("Overlap detected in groups between splits!")

        train_gen = set(table.iloc[train_idx]["genome_id"])
        calib_gen = set(table.iloc[calib_idx]["genome_id"])
        test_gen = set(table.iloc[test_idx]["genome_id"])

        if not train_gen.isdisjoint(calib_gen) or not train_gen.isdisjoint(test_gen) or not calib_gen.isdisjoint(test_gen):
            raise ValueError("Overlap detected in genome_id between splits!")

        table.loc[train_idx, "split"] = "train"
        table.loc[calib_idx, "split"] = "calibration"
        table.loc[test_idx, "split"] = "test"

        metadata["evaluation_status"] = "grouped-evaluation"
        metadata["split_method"] = "grouped train/calibration/test split"
        metadata["groups_used"] = True
        metadata["groups_path_provided"] = Path(groups_path).name
        metadata["number_of_groups_total"] = len(set(table["group_id"]))
        metadata["number_of_groups_train"] = len(train_g)
        metadata["number_of_groups_calibration"] = len(calib_g)
        metadata["number_of_groups_test"] = len(test_g)
        metadata["group_overlap_check"] = "passed"
        metadata["genome_overlap_check"] = "passed"

        return table, metadata

def _safe_float(v):
    if v is None:
        return None
    val = float(v)
    if math.isnan(val) or math.isinf(val):
        return None
    return val

def _evaluate(y_true: np.ndarray, probability: np.ndarray) -> dict[str, Any]:
    if len(y_true) == 0:
        return {}
    predicted = (probability >= 0.5).astype(int)
    both_classes = len(np.unique(y_true)) == 2

    try:
        cm = confusion_matrix(y_true, predicted, labels=[0, 1]).tolist()
    except Exception:
        cm = None

    return {
        "balanced_accuracy": _safe_float(balanced_accuracy_score(y_true, predicted)) if both_classes else None,
        "precision": _safe_float(precision_score(y_true, predicted, zero_division=0)),
        "recall_resistant": _safe_float(recall_score(y_true, predicted, pos_label=1, zero_division=0)),
        "recall_susceptible": _safe_float(recall_score(y_true, predicted, pos_label=0, zero_division=0)),
        "f1": _safe_float(f1_score(y_true, predicted, zero_division=0)),
        "auroc": _safe_float(roc_auc_score(y_true, probability)) if both_classes else None,
        "pr_auc": _safe_float(average_precision_score(y_true, probability)) if both_classes else None,
        "brier_score": _safe_float(brier_score_loss(y_true, probability)),
        "confusion_matrix": cm
    }

def train_drug_model(
    table: pd.DataFrame,
    feature_columns: list[str],
    antibiotic: str,
    output_directory: Path,
    seed: int = 42,
    c_reg: float = 1.0,
    calibration_method: str = "sigmoid",
    split_metadata: dict[str, Any] | None = None,
) -> tuple[dict[str, Any], pd.DataFrame]:

    if calibration_method not in ["sigmoid", "isotonic", "none"]:
        raise ValueError(f"Invalid calibration_method: {calibration_method}")

    train = table[table["split"] == "train"]
    calibration = table[table["split"] == "calibration"]
    test = table[table["split"] == "test"]

    counts = {
        "train": {"susceptible": int((train["target"] == 0).sum()), "resistant": int((train["target"] == 1).sum())},
        "calibration": {"susceptible": int((calibration["target"] == 0).sum()), "resistant": int((calibration["target"] == 1).sum())},
        "test": {"susceptible": int((test["target"] == 0).sum()), "resistant": int((test["target"] == 1).sum())},
    }

    if min(counts["train"].values()) == 0 or min(counts["test"].values()) == 0:
        raise ValueError("Train or test split missing classes.")

    metadata: dict[str, Any] = {
        "antibiotic": antibiotic,
        "status": "skipped",
        "feature_count": len(feature_columns),
        "test_samples": len(test),
        "test_class_distribution": counts["test"],
        "base_estimator_training_samples": len(train),
        "base_estimator_class_distribution": counts["train"],
        "calibration_samples": len(calibration),
        "calibration_class_distribution": counts["calibration"],
        "hyperparameters": {"C": c_reg, "class_weight": "balanced", "solver": "liblinear"},
        "random_state": seed,
        "training_timestamp": pd.Timestamp.now("UTC").isoformat(),
        "software_versions": {"joblib": joblib.__version__},
        "scikit_learn_version": "unknown"
    }
    import sklearn
    metadata["scikit_learn_version"] = sklearn.__version__

    if split_metadata:
        metadata.update(split_metadata)
    else:
        metadata["split_method"] = "unknown"
        metadata["evaluation_status"] = "unknown"

    base = LogisticRegression(
        C=c_reg,
        class_weight="balanced",
        solver="liblinear",
        max_iter=2000,
        random_state=seed,
    )

    base.fit(train[feature_columns].to_numpy(), train["target"].to_numpy())

    calibrator = None
    applied_calibration = "none"

    if calibration_method != "none" and len(calibration) > 0 and min(counts["calibration"].values()) >= 2:
        scores_cal = base.decision_function(calibration[feature_columns].to_numpy())
        y_cal = calibration["target"].to_numpy()

        if calibration_method == "sigmoid":
            calibrator = RobustSigmoidCalibrator()
            calibrator.fit(scores_cal, y_cal)
            applied_calibration = "sigmoid"
        elif calibration_method == "isotonic":
            calibrator = IsotonicRegression(out_of_bounds="clip")
            calibrator.fit(scores_cal, y_cal)
            applied_calibration = "isotonic"

    metadata["calibration_method"] = applied_calibration

    output_directory.mkdir(parents=True, exist_ok=True)

    joblib.dump(
        {
            "antibiotic": antibiotic,
            "base_estimator": base,
            "calibrator": calibrator,
            "feature_columns": feature_columns,
            "calibration_method": metadata["calibration_method"],
            "hyperparameters": metadata["hyperparameters"],
            "random_state": metadata["random_state"],
            "split_method": metadata.get("split_method", "unknown"),
            "evaluation_status": metadata.get("evaluation_status", "unknown"),
        },
        output_directory / f"{antibiotic.replace('/', '_')}.joblib"
    )

    predictions = test[["genome_id", "label", "target"]].copy()
    if len(test):
        scores_test = base.decision_function(test[feature_columns].to_numpy())
        if applied_calibration == "sigmoid":
            probs = calibrator.predict_proba(scores_test)[:, 1]
        elif applied_calibration == "isotonic":
            probs = calibrator.predict(scores_test)
        else:
            probs = base.predict_proba(test[feature_columns].to_numpy())[:, 1]
        predictions["probability_resistant"] = probs
    else:
        predictions["probability_resistant"] = pd.Series(dtype=float)

    predictions["prediction"] = (predictions["probability_resistant"] >= 0.5).astype(int)
    predictions["antibiotic"] = antibiotic

    metadata["status"] = "trained"
    metrics = _evaluate(
        predictions["target"].to_numpy(dtype=int),
        predictions["probability_resistant"].to_numpy(dtype=float)
    )
    metrics["test_size"] = len(test)
    metrics["test_class_distribution"] = counts["test"]
    metrics["evaluation_status"] = metadata.get("evaluation_status", "unknown")
    metrics["split_method"] = metadata.get("split_method", "unknown")
    metrics["calibration_method"] = applied_calibration

    metadata["test_metrics"] = metrics

    return metadata, predictions

class AMRModel:
    def __init__(self, model_path: Path | str):
        data = joblib.load(model_path)
        self.antibiotic = data["antibiotic"]
        self.base_estimator = data["base_estimator"]
        self.calibrator = data.get("calibrator", None)
        self.feature_columns = data["feature_columns"]
        self.calibration_method = data.get("calibration_method", "none")
        self.split_method = data.get("split_method", "unknown")
        self.evaluation_status = data.get("evaluation_status", "unknown")
        self.random_state = data.get("random_state", None)

    def predict(self, features_df: pd.DataFrame) -> pd.DataFrame:
        df = features_df.copy()

        extra_features = [c for c in df.columns if c not in self.feature_columns and c != "genome_id"]
        if extra_features:
            logging.info(f"Ignoring extra features: {extra_features}")

        missing_features = [c for c in self.feature_columns if c not in df.columns]
        if missing_features:
            missing_df = pd.DataFrame(0.0, index=df.index, columns=missing_features)
            df = pd.concat([df, missing_df], axis=1)

        X = df[self.feature_columns].astype(float).to_numpy()

        scores = self.base_estimator.decision_function(X)

        if self.calibration_method == "sigmoid" and self.calibrator is not None:
            probabilities = self.calibrator.predict_proba(scores)[:, 1]
        elif self.calibration_method == "isotonic" and self.calibrator is not None:
            probabilities = self.calibrator.predict(scores)
        else:
            probabilities = self.base_estimator.predict_proba(X)[:, 1]

        predictions = (probabilities >= 0.5).astype(int)
        confidences = np.abs(probabilities - 0.5) * 2

        def get_evidence(row):
            return ",".join([f for f in self.feature_columns if row[f] > 0])

        results = pd.DataFrame({
            "genome_id": df.get("genome_id", np.arange(len(df))),
            "antibiotic": self.antibiotic,
            "probability_resistant": probabilities,
            "prediction": predictions,
            "confidence": confidences,
        })
        results["evidence_features"] = df[self.feature_columns].astype(float).apply(get_evidence, axis=1)

        return results
