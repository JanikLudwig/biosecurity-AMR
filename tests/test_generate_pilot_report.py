import pytest
import pandas as pd
import json
from pathlib import Path
from unittest.mock import patch

from scripts.generate_pilot_report import main

def test_generate_pilot_report_no_absolute_paths(tmp_path):
    # Setup mock data for the report script
    
    # 1. Labels
    labels_file = tmp_path / "labels.csv"
    pd.DataFrame({
        "genome_id": ["id1", "id1", "id2"],
        "antibiotic": ["cefoxitin", "ciprofloxacin", "cefoxitin"],
        "label": ["R", "S", "R"]
    }).to_csv(labels_file, index=False)
    
    # 2. Manifest
    manifest_file = tmp_path / "manifest.csv"
    pd.DataFrame({
        "genome_id": ["id1", "id2"]
    }).to_csv(manifest_file, index=False)
    
    # 3. Run report
    run_report_file = tmp_path / "runs.csv"
    pd.DataFrame({
        "original_genome_id": ["id1", "id2"],
        "status": ["success", "failed"],
        "runtime_seconds": [1.5, None]
    }).to_csv(run_report_file, index=False)
    
    # 4. Feature Summary
    feature_summary_file = tmp_path / "feature_summary.json"
    with open(feature_summary_file, "w") as f:
        json.dump({
            "number_of_gene_features": 10,
            "number_of_mutation_features": 2,
            "number_of_features": 12,
            "number_of_genomes": 2,
            "number_of_zero_hit_genomes": 0
        }, f)
        
    # 5. Features
    features_file = tmp_path / "features.csv"
    pd.DataFrame({
        "genome_id": ["id1", "id2"],
        "gene::blaZ": [1, 0],
        "mutation::gyrA::S84L": [1, 0]
    }).to_csv(features_file, index=False)
    
    # Output paths
    out_json = tmp_path / "summary.json"
    out_md = tmp_path / "summary.md"
    
    # Execute
    args = [
        "--labels", str(labels_file),
        "--manifest", str(manifest_file),
        "--run-report", str(run_report_file),
        "--feature-summary", str(feature_summary_file),
        "--features", str(features_file),
        "--output-json", str(out_json),
        "--output-markdown", str(out_md)
    ]
    
    # Make sure it runs without exceptions
    main(args)
    
    # Verify outputs were created and contain expected relative info
    assert out_json.exists()
    assert out_md.exists()
    
    with open(out_json, "r", encoding="utf-8") as f:
        data = json.load(f)
        
    assert data["pilot_genome_count"] == 2
    assert data["amrfinder_success"] == 1
    assert data["amrfinder_failed"] == 1
    assert data["genomes_with_mutations"] == 1
    assert data["top_mutations"] == {"mutation::gyrA::S84L": 1}
    
    # Check proper path sanitization
    assert data["data_sources"]["labels"] == "labels.csv"
    assert data["data_sources"]["manifest"] == "manifest.csv"
    
    # Ensure no local absolute paths in the output
    md_content = out_md.read_text(encoding="utf-8")
    
    # The markdown file contains the safe filename, but not the full absolute path.
    assert "labels.csv" in md_content
    assert str(tmp_path) not in md_content
    assert "C:\\Users\\" not in md_content
    
    assert str(tmp_path) not in json.dumps(data)
    assert "C:\\Users\\" not in json.dumps(data)
