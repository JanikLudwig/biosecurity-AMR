import pytest
from pathlib import Path
import pandas as pd
import json
import csv
from unittest.mock import patch, MagicMock

from genome_firewall.prediction import (
    validate_input_fasta,
    parse_amrfinder_output,
    align_features,
    predict_antibiotics,
    load_models
)
from genome_firewall.modeling.baseline import AMRModel

def test_validate_fasta_missing(tmp_path):
    missing = tmp_path / "missing.fna"
    with pytest.raises(FileNotFoundError):
        validate_input_fasta(missing)

def test_validate_fasta_empty(tmp_path):
    empty = tmp_path / "empty.fna"
    empty.touch()
    with pytest.raises(ValueError, match="FASTA file is empty"):
        validate_input_fasta(empty)

def test_validate_fasta_invalid_extension(tmp_path):
    invalid = tmp_path / "file.txt"
    invalid.write_text(">seq\nACGT")
    with pytest.raises(ValueError, match="Invalid FASTA extension"):
        validate_input_fasta(invalid)

def test_validate_fasta_no_header(tmp_path):
    f = tmp_path / "nohead.fna"
    f.write_text("ACGT\n")
    with pytest.raises(ValueError, match="No FASTA header found"):
        validate_input_fasta(f)

def test_validate_fasta_invalid_chars(tmp_path):
    f = tmp_path / "bad.fna"
    f.write_text(">seq\nACGT123\n")
    with pytest.raises(ValueError, match="Sequence contains invalid nucleotide characters"):
        validate_input_fasta(f)

def test_validate_fasta_valid(tmp_path):
    f = tmp_path / "good.fna"
    f.write_text(">contig1\nACGTN-M\n>contig2\nTGCA")
    # Should not raise
    validate_input_fasta(f)

def test_parse_amrfinder_output(tmp_path):
    tsv = tmp_path / "amr.tsv"
    # Create mock TSV matching AMRFinderPlus
    # Gene symbol, Element name, Element type, Element subtype, Method
    content = "Gene symbol\tElement type\tElement subtype\tMethod\tClass\tSubclass\n"
    content += "mecA\tAMR\t\t\tBETA-LACTAM\t\n"
    content += "gyrA_S84L\tPOINT\tPOINT\tPOINT\tQUINOLONE\t\n"
    tsv.write_text(content)

    df, features = parse_amrfinder_output(tsv, "genome1")
    assert "genome_id" in df.columns
    assert df["genome_id"].iloc[0] == "genome1"

    expected_features = ["gene::mecA", "mutation::gyrA::S84L"]
    for ef in expected_features:
        assert ef in df.columns
        assert df[ef].iloc[0] == 1

    assert set(features) == set(expected_features)

def test_parse_amrfinder_empty_hits(tmp_path):
    tsv = tmp_path / "amr.tsv"
    tsv.write_text("Gene symbol\tElement type\tElement subtype\tMethod\tClass\tSubclass\n")
    df, features = parse_amrfinder_output(tsv, "genome2")
    assert df.columns.tolist() == ["genome_id"]
    assert features == []

def test_align_features():
    df = pd.DataFrame({"genome_id": ["g1"], "gene::mecA": [1], "unknown::feature": [1]})

    # Mock model
    mock_model = MagicMock(spec=AMRModel)
    mock_model.feature_columns = ["gene::mecA", "gene::blaZ", "mutation::gyrA::S84L"]

    models = {"test": mock_model}

    df_aligned, rec_count, unk_count = align_features(df, models)
    assert rec_count == 1
    assert unk_count == 1

def test_predict_antibiotics():
    df = pd.DataFrame({"genome_id": ["g1"], "gene::mecA": [1]})

    mock_model1 = MagicMock(spec=AMRModel)
    mock_model1.calibration_method = "sigmoid"
    mock_model1.evaluation_status = "development-only"
    res1 = pd.DataFrame([{
        "genome_id": "g1",
        "antibiotic": "cefoxitin",
        "probability_resistant": 0.95,
        "prediction": 1,
        "confidence": 0.90,
        "evidence_features": "gene::mecA"
    }])
    mock_model1.predict.return_value = res1

    mock_model2 = MagicMock(spec=AMRModel)
    mock_model2.calibration_method = "sigmoid"
    mock_model2.evaluation_status = "development-only"
    res2 = pd.DataFrame([{
        "genome_id": "g1",
        "antibiotic": "ciprofloxacin",
        "probability_resistant": 0.10,
        "prediction": 0,
        "confidence": 0.80,
        "evidence_features": ""
    }])
    mock_model2.predict.return_value = res2

    models = {"cefoxitin": mock_model1, "ciprofloxacin": mock_model2}
    preds = predict_antibiotics(df, models)

    assert len(preds) == 2
    assert preds[0]["prediction"] == "R"
    assert preds[0]["evidence_features"] == ["gene::mecA"]

    assert preds[1]["prediction"] == "S"
    assert preds[1]["evidence_features"] == []

@patch("genome_firewall.prediction.run_single_genome")
@patch("genome_firewall.prediction.load_models")
def test_predict_genome_e2e(mock_load_models, mock_run_single, tmp_path):
    from genome_firewall.prediction import predict_genome

    fasta = tmp_path / "test.fna"
    fasta.write_text(">seq1\nACGT\n")

    # Mock AMRFinderPlus run to create the TSV
    def fake_run(*args, **kwargs):
        tsv_path = kwargs['output_dir'] / f"{kwargs['raw_id']}.tsv"
        tsv_path.parent.mkdir(parents=True, exist_ok=True)
        tsv_path.write_text("Gene symbol\tElement type\tElement subtype\tMethod\tClass\tSubclass\n"
                            "mecA\tAMR\t\t\tBETA-LACTAM\t\n")
        return {'status': 'success'}

    mock_run_single.side_effect = fake_run

    # Mock models
    mock_model1 = MagicMock(spec=AMRModel)
    mock_model1.calibration_method = "sigmoid"
    mock_model1.evaluation_status = "development-only"
    mock_model1.feature_columns = ["gene::mecA"]
    res1 = pd.DataFrame([{
        "genome_id": "test",
        "antibiotic": "cefoxitin",
        "probability_resistant": 0.95,
        "prediction": 1,
        "confidence": 0.90,
        "evidence_features": "gene::mecA"
    }])
    mock_model1.predict.return_value = res1
    mock_load_models.return_value = {"cefoxitin": mock_model1}

    out_dir = tmp_path / "out"

    exit_code = predict_genome(
        fasta_path=fasta,
        model_dir=tmp_path / "models", # doesn't matter, mocked
        output_dir=out_dir,
        organism="Staphylococcus_aureus",
        docker_image="fake",
        backend="docker",
        genome_id="test"
    )

    assert exit_code == 0

    # Check output files
    run_dir = list((out_dir / "test").glob("*"))[0]
    assert (run_dir / "prediction.json").exists()
    assert (run_dir / "prediction.csv").exists()

    with open(run_dir / "prediction.json") as f:
        data = json.load(f)
        assert data["genome_id"] == "test"
        assert len(data["predictions"]) == 1
        assert data["predictions"][0]["prediction"] == "R"

    df = pd.read_csv(run_dir / "prediction.csv")
    assert len(df) == 1
    assert df["prediction"].iloc[0] == "R"
