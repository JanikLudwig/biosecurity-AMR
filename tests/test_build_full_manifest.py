import pytest
from pathlib import Path
import pandas as pd
import json
import subprocess
import sys

from scripts.build_full_manifest import main, is_computational, find_fasta

def test_cli_arguments(tmp_path):
    labels = tmp_path / "labels.tsv"
    labels.write_text("Dummy\n") # Avoid EmptyDataError
    fasta_dir = tmp_path / "genomes"
    fasta_dir.mkdir()
    
    # Missing required args
    with pytest.raises(SystemExit):
        main([])
        
    # Valid arguments but empty file (will exit 1 due to missing columns)
    assert main(["--labels", str(labels), "--fasta-dir", str(fasta_dir)]) == 1

def test_missing_input_files(tmp_path):
    # Fasta dir exists, labels missing
    fasta_dir = tmp_path / "genomes"
    fasta_dir.mkdir()
    assert main(["--labels", str(tmp_path / "missing.tsv"), "--fasta-dir", str(fasta_dir)]) == 1
    
    # Labels exist, fasta dir missing
    labels = tmp_path / "labels.tsv"
    labels.write_text("Dummy\n")
    assert main(["--labels", str(labels), "--fasta-dir", str(tmp_path / "missing_dir")]) == 1

def test_rs_normalization(tmp_path):
    labels = tmp_path / "labels.tsv"
    fasta_dir = tmp_path / "genomes"
    fasta_dir.mkdir()
    
    # Create valid mock data
    df = pd.DataFrame({
        "Genome ID": ["1.1", "1.2", "1.3", "1.4"],
        "Antibiotic": ["cefoxitin"] * 4,
        "Resistant Phenotype": ["Resistant", "Susceptible", "RESISTANT ", " Intermediate"],
        "Evidence": ["Lab", "Lab", "Lab", "Lab"],
        "Computational Method": ["", "", "", ""]
    })
    df.to_csv(labels, sep="\t", index=False)
    
    # Create fastas
    for gid in ["1.1", "1.2", "1.3"]:
        f = fasta_dir / f"{gid}.fna"
        f.write_text(">seq\nATGC")
        
    out_manifest = tmp_path / "manifest.csv"
    out_labels = tmp_path / "labels.csv.gz"
    
    assert main(["--labels", str(labels), "--fasta-dir", str(fasta_dir), 
                 "--manifest-out", str(out_manifest), "--labels-out", str(out_labels)]) == 0
                 
    res_labels = pd.read_csv(out_labels, dtype=str)
    assert len(res_labels) == 3
    assert set(res_labels["label"]) == {"R", "S"}

def test_conflict_exclusion(tmp_path):
    labels = tmp_path / "labels.tsv"
    fasta_dir = tmp_path / "genomes"
    fasta_dir.mkdir()
    
    df = pd.DataFrame({
        "Genome ID": ["1.1", "1.1", "1.2", "1.2", "1.3"],
        "Antibiotic": ["cefoxitin", "cefoxitin", "ciprofloxacin", "ciprofloxacin", "cefoxitin"],
        "Resistant Phenotype": ["R", "S", "R", "R", "S"], # 1.1 conflict, 1.2 valid, 1.3 valid
        "Evidence": ["Lab"] * 5,
        "Computational Method": [""] * 5
    })
    df.to_csv(labels, sep="\t", index=False)
    
    for gid in ["1.1", "1.2", "1.3"]:
        f = fasta_dir / f"{gid}.fna"
        f.write_text(">seq\nATGC")
        
    out_labels = tmp_path / "labels.csv.gz"
    main(["--labels", str(labels), "--fasta-dir", str(fasta_dir), 
          "--manifest-out", str(tmp_path / "m.csv"), "--labels-out", str(out_labels)])
          
    res_labels = pd.read_csv(out_labels, dtype=str)
    # 1.1 should be excluded
    assert "1.1" not in res_labels["genome_id"].values
    assert "1.2" in res_labels["genome_id"].values
    assert "1.3" in res_labels["genome_id"].values
    assert len(res_labels) == 2

def test_computational_prediction_exclusion():
    # Test the is_computational logic directly
    assert is_computational({"Computational Method": "AMRFinder"}) == True
    assert is_computational({"Evidence": "Predicted"}) == True
    assert is_computational({"Laboratory Typing Method": "in silico analysis"}) == True
    
    assert is_computational({
        "Computational Method": "",
        "Evidence": "Laboratory Method",
        "Laboratory Typing Method": "MIC"
    }) == False
    
    assert is_computational({"Computational Method": pd.NA}) == False

def test_exact_genome_id_matching(tmp_path):
    fasta_dir = tmp_path / "genomes"
    fasta_dir.mkdir()
    
    # 1.1 is correct, 1.10 shouldn't match 1.1
    (fasta_dir / "1.10.fna").write_text(">seq\nATGC")
    (fasta_dir / "1.1.fna").write_text(">seq\nATGC")
    
    # Test function directly
    assert find_fasta("1.1", fasta_dir).name == "1.1.fna"
    assert find_fasta("1.10", fasta_dir).name == "1.10.fna"
    assert find_fasta("1.2", fasta_dir) is None

def test_windows_paths_with_spaces(tmp_path):
    space_dir = tmp_path / "my genomes folder"
    space_dir.mkdir()
    
    labels = tmp_path / "my labels file.tsv"
    
    df = pd.DataFrame({
        "Genome ID": ["1.1"],
        "Antibiotic": ["cefoxitin"],
        "Resistant Phenotype": ["R"]
    })
    df.to_csv(labels, sep="\t", index=False)
    
    (space_dir / "1.1.fna").write_text(">seq\nATGC")
    
    out_manifest = tmp_path / "out manifest.csv"
    out_labels = tmp_path / "out labels.csv.gz"
    
    assert main([
        "--labels", str(labels), 
        "--fasta-dir", str(space_dir), 
        "--manifest-out", str(out_manifest), 
        "--labels-out", str(out_labels)
    ]) == 0
    
    assert out_manifest.exists()
    assert out_labels.exists()

def test_no_hardcoded_paths():
    script_path = Path("scripts/build_full_manifest.py")
    content = script_path.read_text(encoding="utf-8")
    
    assert "C:\\Users\\" not in content
    assert "C:/Users/" not in content
    assert "benja" not in content
