import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from genome_firewall.modeling.baseline import (
    build_training_dataset,
    assign_development_split,
    train_drug_model,
    AMRModel,
)

@pytest.fixture
def synthetic_data(tmp_path):
    features_path = tmp_path / "features.csv"
    labels_path = tmp_path / "labels.csv"

    features_df = pd.DataFrame({
        "genome_id": [f"G{i}" for i in range(20)],
        "gene::A": [i % 2 for i in range(20)],
        "mutation::B": [(i+1) % 2 for i in range(20)],
        "gene::C": [1] * 10 + [0] * 10,
    })
    features_df.to_csv(features_path, index=False)

    labels_df = pd.DataFrame({
        "genome_id": [f"G{i}" for i in range(20)],
        "antibiotic": ["cefoxitin"] * 20,
        "label": ["R", "S"] * 10,
        "Evidence": ["Lab"] * 20,
    })
    labels_df.to_csv(labels_path, index=False)

    return features_path, labels_path

def test_feature_label_join_and_encoding(synthetic_data):
    features_path, labels_path = synthetic_data

    table, feature_columns = build_training_dataset(features_path, labels_path, "cefoxitin")

    assert "target" in table.columns
    assert len(table) == 20
    assert "target" not in feature_columns
    assert "label" not in feature_columns

def test_train_and_serialize_model_sigmoid(synthetic_data, tmp_path):
    features_path, labels_path = synthetic_data

    table, feature_columns = build_training_dataset(features_path, labels_path, "cefoxitin")
    table, split_metadata = assign_development_split(table, train_fraction=0.6, calibration_fraction=0.2, test_fraction=0.2)

    output_dir = tmp_path / "models"
    metadata, predictions = train_drug_model(
        table, feature_columns, antibiotic="cefoxitin", output_directory=output_dir, c_reg=1.0, calibration_method="sigmoid"
    )

    assert metadata["status"] == "trained"
    assert metadata["calibration_method"] == "sigmoid"
    assert (output_dir / "cefoxitin.joblib").exists()

    model = AMRModel(output_dir / "cefoxitin.joblib")
    assert model.antibiotic == "cefoxitin"
    assert model.calibration_method == "sigmoid"

    res = model.predict(table)
    assert "probability_resistant" in res.columns
    assert res["probability_resistant"].max() <= 1.0
    assert res["probability_resistant"].min() >= 0.0

def test_train_and_serialize_model_none(synthetic_data, tmp_path):
    features_path, labels_path = synthetic_data
    table, feature_columns = build_training_dataset(features_path, labels_path, "cefoxitin")
    table, split_metadata = assign_development_split(table, train_fraction=0.6, calibration_fraction=0.2, test_fraction=0.2)
    output_dir = tmp_path / "models"
    metadata, predictions = train_drug_model(
        table, feature_columns, antibiotic="cefoxitin", output_directory=output_dir, c_reg=1.0, calibration_method="none"
    )
    assert metadata["calibration_method"] == "none"

def test_train_and_serialize_model_isotonic(synthetic_data, tmp_path):
    features_path, labels_path = synthetic_data
    table, feature_columns = build_training_dataset(features_path, labels_path, "cefoxitin")
    table, split_metadata = assign_development_split(table, train_fraction=0.6, calibration_fraction=0.2, test_fraction=0.2)
    output_dir = tmp_path / "models"
    metadata, predictions = train_drug_model(
        table, feature_columns, antibiotic="cefoxitin", output_directory=output_dir, c_reg=1.0, calibration_method="isotonic"
    )
    assert metadata["calibration_method"] == "isotonic"

def test_group_split_success(tmp_path):
    df = pd.DataFrame({
        "genome_id": [f"G{i}" for i in range(20)],
        "target": [0, 1] * 10
    })
    groups_file = tmp_path / "groups.csv"
    # Create 5 groups of 4 genomes each. Ensure both classes are in groups to allow splits.
    group_ids = [f"group{i//4}" for i in range(20)]
    groups_file.write_text("genome_id,group_id\n" + "\n".join([f"G{i},{g}" for i, g in zip(range(20), group_ids)]))

    table, split_metadata = assign_development_split(df, groups_path=groups_file)
    assert split_metadata["groups_used"] is True

    train_ids = set(table[table["split"]=="train"]["genome_id"])
    calib_ids = set(table[table["split"]=="calibration"]["genome_id"])
    test_ids = set(table[table["split"]=="test"]["genome_id"])

    assert train_ids.isdisjoint(calib_ids)
    assert train_ids.isdisjoint(test_ids)
    assert calib_ids.isdisjoint(test_ids)

    train_groups = set(table[table["split"]=="train"]["group_id"])
    calib_groups = set(table[table["split"]=="calibration"]["group_id"])
    test_groups = set(table[table["split"]=="test"]["group_id"])

    assert train_groups.isdisjoint(calib_groups)
    assert train_groups.isdisjoint(test_groups)
    assert calib_groups.isdisjoint(test_groups)

def test_group_split_missing_file():
    df = pd.DataFrame({"genome_id": ["G1"], "target": [0]})
    with pytest.raises(FileNotFoundError):
        assign_development_split(df, groups_path="does_not_exist.csv")

def test_group_split_missing_cols(tmp_path):
    df = pd.DataFrame({"genome_id": ["G1"], "target": [0]})
    groups_file = tmp_path / "groups.csv"
    groups_file.write_text("id,group\n1,1")
    with pytest.raises(ValueError, match="must contain"):
        assign_development_split(df, groups_path=groups_file)

def test_group_split_missing_mapping(tmp_path):
    df = pd.DataFrame({"genome_id": ["G1", "G2"], "target": [0, 1]})
    groups_file = tmp_path / "groups.csv"
    groups_file.write_text("genome_id,group_id\nG1,group1\n")
    with pytest.raises(ValueError, match="Genomes missing from groups file"):
        assign_development_split(df, groups_path=groups_file)

def test_group_split_duplicate_genomes(tmp_path):
    df = pd.DataFrame({"genome_id": ["G1"], "target": [0]})
    groups_file = tmp_path / "groups.csv"
    groups_file.write_text("genome_id,group_id\nG1,group1\nG1,group2\n")
    with pytest.raises(ValueError, match="Duplicate genome_id"):
        assign_development_split(df, groups_path=groups_file)

from unittest.mock import patch

def test_group_split_overlap_validation(tmp_path):
    df = pd.DataFrame({
        "genome_id": [f"G{i}" for i in range(20)],
        "target": [0, 1] * 10
    })
    groups_file = tmp_path / "groups.csv"
    group_ids = [f"group{i//4}" for i in range(20)]
    groups_file.write_text("genome_id,group_id\n" + "\n".join([f"G{i},{g}" for i, g in zip(range(20), group_ids)]))

    # Mock GroupShuffleSplit to return overlapping splits
    with patch('genome_firewall.modeling.baseline.GroupShuffleSplit.split') as mock_split:
        # Force the split to return the exact same indices for train and test
        mock_split.side_effect = lambda *args, **kwargs: iter([(np.array([0,1,2,3,4,5,6,7]), np.array([0,1,2,3,4,5,6,7]))])
        with pytest.raises(ValueError, match="Overlap detected in groups between splits!"):
            assign_development_split(df, groups_path=groups_file)
