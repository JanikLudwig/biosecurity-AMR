import pytest
import pandas as pd
from pathlib import Path
from genome_firewall.import_local_dataset import validate_fasta, run_import

class DummyArgs:
    def __init__(self, tsv, fasta_dir, manifest, labels, reports_dir, seed=42):
        self.tsv = tsv
        self.fasta_dir = fasta_dir
        self.manifest = manifest
        self.labels = labels
        self.reports_dir = reports_dir
        self.seed = seed

def test_full_pipeline_synthetic(tmp_path):
    tsv_path = tmp_path / "labels.tsv"
    fasta_dir = tmp_path / "genomes"
    fasta_dir.mkdir()
    
    # 1. TSV Schema, Normalization, Conflicts, Exclude Computational
    tsv_content = (
        "Genome ID\tGenome Name\tAntibiotic\tResistant Phenotype\tEvidence\n"
        "1.1\tGen1\tcefoxitin\tResistant\tLaboratory Method\n"  # Valid R
        "1.1\tGen1\tciprofloxacin\tsusceptible\tLaboratory Method\n"  # Valid S
        "2.2\tGen2\tcefoxitin\tIntermediate\tLaboratory Method\n"  # Excluded I
        "3.3\tGen3\tcefoxitin\tResistant\tLaboratory Method\n"  # Conflict R
        "3.3\tGen3\tcefoxitin\tSusceptible\tLaboratory Method\n"  # Conflict S -> both dropped
        "4.4\tGen4\tcefoxitin\tResistant\tComputational Prediction\n"  # Excluded Evidence
        "5.5\tGen5\terythromycin\tNonSusceptible\tLaboratory Method\n" # Excluded NS
        "6.1\tGen6\tcefoxitin\tResistant\tLaboratory Method\n"  # valid
        "6.10\tGen6.10\tcefoxitin\tResistant\tLaboratory Method\n"  # valid, to test mapping 6.1 vs 6.10
    )
    # create more genomes so we can test the 100 limit if needed, but for now just basic tests
    # let's create 105 valid genomes
    extra = ""
    for i in range(100, 205):
        extra += f"{i}.1\tGen{i}\tcefoxitin\tResistant\tLaboratory Method\n"
    
    tsv_path.write_text(tsv_content + extra, encoding="utf-8")
    
    # 2. FASTAs
    def write_fasta(gid, content=">h\n" + "A"*1000):
        (fasta_dir / f"{gid}.fna").write_text(content)
        
    write_fasta("1.1")
    write_fasta("2.2") # no R/S labels
    write_fasta("3.3") # conflict
    write_fasta("4.4") # computational
    write_fasta("6.1")
    write_fasta("6.10")
    
    # Duplicate FASTA test
    (fasta_dir / "6.1.fasta").write_text(">h\n" + "A"*1000) 
    
    for i in range(100, 205):
        write_fasta(f"{i}.1")
        
    # invalid FASTA
    (fasta_dir / "101.1.fna").write_text("invalid")
    
    manifest_out = tmp_path / "out/manifest.csv"
    labels_out = tmp_path / "out/labels.csv.gz"
    reports_dir = tmp_path / "out/reports"
    
    args = DummyArgs(str(tsv_path), str(fasta_dir), str(manifest_out), str(labels_out), str(reports_dir))
    
    code = run_import(args)
    assert code == 0
    
    # Manifest compatibility
    man_df = pd.read_csv(manifest_out, dtype=str)
    assert "genome_id" in man_df.columns
    assert "fasta_path" in man_df.columns
    
    # Exactly 100 selected
    assert len(man_df) == 100
    assert len(man_df["genome_id"].unique()) == 100
    
    # Check 6.1 vs 6.10 exact mapping: 6.1 had duplicate FASTA (6.1.fna and 6.1.fasta) so it should be skipped
    assert "6.1" not in man_df["genome_id"].values
    assert "6.10" in man_df["genome_id"].values
    
    # Check invalid fasta
    assert "101.1" not in man_df["genome_id"].values
    
    # Check determinism: run again with same seed, result should match
    manifest_out2 = tmp_path / "out/manifest2.csv"
    args2 = DummyArgs(str(tsv_path), str(fasta_dir), str(manifest_out2), str(labels_out), str(reports_dir))
    run_import(args2)
    man_df2 = pd.read_csv(manifest_out2, dtype=str)
    assert man_df["genome_id"].tolist() == man_df2["genome_id"].tolist()

def test_validate_fasta_valid(tmp_path):
    f = tmp_path / "valid.fna"
    f.write_text(">header\nATGC\n")
    res = validate_fasta(f)
    assert not res["valid"]  # too short
    
    f.write_text(">header\n" + "A"*1500 + "\n")
    res = validate_fasta(f)
    assert res["valid"]
    assert res["size"] > 0
    assert res["seq_len"] == 1500

def test_validate_fasta_invalid_empty(tmp_path):
    f = tmp_path / "empty.fna"
    f.write_text("")
    res = validate_fasta(f)
    assert not res["valid"]
    assert res["reason"] == "empty"

def test_validate_fasta_invalid_no_header(tmp_path):
    f = tmp_path / "noheader.fna"
    f.write_text("ATGC\nATGC")
    res = validate_fasta(f)
    assert not res["valid"]
    assert res["reason"] == "no_header"

def test_windows_paths(tmp_path):
    path_with_space = tmp_path / "my folder"
    path_with_space.mkdir()
    f = path_with_space / "test.fna"
    f.write_text(">header\n" + "A"*1000)
    res = validate_fasta(f)
    assert res["valid"]
