import pytest
import os
import json
import subprocess
from pathlib import Path
from unittest.mock import patch, MagicMock

import genome_firewall.run_amrfinder as run_amrfinder

def test_workers_must_be_positive_main(monkeypatch, tmp_path):
    # Mock sys.exit to raise an exception instead of exiting the test
    with patch("sys.exit") as mock_exit:
        mock_exit.side_effect = SystemExit
        with patch("sys.argv", ["python", "--workers", "0", "--threads", "2"]):
            with pytest.raises(SystemExit):
                run_amrfinder.main()

def test_warnings_cpu_count_and_memory(monkeypatch, tmp_path, caplog):
    # Mock sys.exit and discover_fastas to just exit cleanly without doing real work
    with patch("sys.exit") as mock_exit:
        mock_exit.side_effect = SystemExit
        with patch("genome_firewall.run_amrfinder.discover_fastas", return_value=[]):
            with patch("os.cpu_count", return_value=4):
                with patch("sys.argv", ["python", "--workers", "5", "--threads", "2"]):
                    with pytest.raises(SystemExit):
                        run_amrfinder.main()
                        
    # Check if warnings were logged
    assert any("greater than os.cpu_count()" in record.message for record in caplog.records)
    assert any("Docker on Windows may run out of memory" in record.message for record in caplog.records)

@patch('genome_firewall.run_amrfinder.subprocess.run')
def test_parallel_execution_success(mock_run, tmp_path):
    # We will mock the behavior to succeed.
    mock_run.return_value = MagicMock(returncode=0, stdout="amrfinder 3.2.1\nDatabase: 2026-05-15.1\n", stderr="")
    
    # We'll create some fake fasta paths
    genomes = [
        ("id1", tmp_path / "id1.fasta"),
        ("id2", tmp_path / "id2.fasta"),
        ("id3", tmp_path / "id3.fasta")
    ]
    
    # We must patch validate_existing_output to pretend the output file was created correctly
    with patch('genome_firewall.run_amrfinder.validate_existing_output') as mock_validate:
        # Initially, output does not exist (skip is False). 
        # Then, when checked after subprocess, it returns True.
        mock_validate.side_effect = [False, True, False, True, False, True]
        
        with patch('pathlib.Path.replace') as mock_replace:
            results = []
            for raw_id, fasta in genomes:
                res = run_amrfinder.run_single_genome(
                    raw_id=raw_id,
                    fasta_path=fasta,
                    output_dir=tmp_path / "out",
                    log_dir=tmp_path / "logs",
                    backend="docker",
                    image="img",
                    organism="org",
                    threads=1,
                    plus=False,
                    force=False
                )
                results.append(res)
            
            assert len(results) == 3
            assert all(r['status'] == 'success' for r in results)
            
def test_deterministic_report_order():
    # results are appended in whatever order threads complete them.
    # main() should sort them by original_genome_id.
    
    genomes = [("id2", Path("")), ("id1", Path("")), ("id3", Path(""))]
    
    # mock discover_fastas to return our list
    # mock run_single_genome to return instantly
    with patch("genome_firewall.run_amrfinder.discover_fastas", return_value=genomes):
        with patch("genome_firewall.run_amrfinder.run_single_genome") as mock_run:
            mock_run.side_effect = lambda raw_id, **kwargs: {'original_genome_id': raw_id, 'status': 'success'}
            
            with patch("genome_firewall.run_amrfinder.write_run_report") as mock_write:
                with patch("sys.exit") as mock_exit:
                    mock_exit.side_effect = SystemExit
                    with patch("sys.argv", ["python", "--workers", "3", "--input-dir", "in", "--output-dir", "out"]):
                        with pytest.raises(SystemExit):
                            run_amrfinder.main()
                            
                # Check that write_run_report was called with SORTED results
                call_args = mock_write.call_args[0]
                results = call_args[1]
                assert [r['original_genome_id'] for r in results] == ["id1", "id2", "id3"]

@patch('genome_firewall.run_amrfinder.subprocess.run')
def test_fail_fast_behavior(mock_run, tmp_path):
    genomes = [("id1", Path("")), ("id2", Path("")), ("id3", Path(""))]
    
    with patch("genome_firewall.run_amrfinder.discover_fastas", return_value=genomes):
        with patch("genome_firewall.run_amrfinder.run_single_genome") as mock_run_single:
            # First one fails, the others shouldn't even be submitted if fail-fast works 
            # (or they get cancelled).
            def side_effect(raw_id, **kwargs):
                if raw_id == "id1":
                    return {'original_genome_id': raw_id, 'status': 'failed'}
                import time
                time.sleep(0.1) # Wait a bit to simulate work and let fail-fast cancel it
                return {'original_genome_id': raw_id, 'status': 'success'}
                
            mock_run_single.side_effect = side_effect
            
            with patch("genome_firewall.run_amrfinder.write_run_report"):
                with patch("sys.exit") as mock_exit:
                    mock_exit.side_effect = SystemExit
                    # Running sequentially (workers=1) makes fail-fast deterministic
                    with patch("sys.argv", ["python", "--workers", "1", "--fail-fast"]):
                        with pytest.raises(SystemExit):
                            run_amrfinder.main()
                            
            # Because it's sequential and fail-fast, only id1 should be processed
            assert mock_run_single.call_count == 1
            assert mock_run_single.call_args[1]['raw_id'] == 'id1'

def test_skip_existing_outputs(tmp_path):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / "id1.tsv").write_text("Header\nData\n")
    
    with patch("genome_firewall.run_amrfinder.subprocess.run") as mock_run:
        res = run_amrfinder.run_single_genome(
            raw_id="id1",
            fasta_path=tmp_path / "id1.fa",
            output_dir=out_dir,
            log_dir=tmp_path / "logs",
            backend="docker",
            image="img",
            organism="org",
            threads=1,
            plus=False,
            force=False
        )
        assert res['status'] == 'skipped'
        assert mock_run.call_count == 0

        # Now with force=True
        mock_run.return_value = MagicMock(returncode=1) # Just mock it to do something
        res = run_amrfinder.run_single_genome(
            raw_id="id1",
            fasta_path=tmp_path / "id1.fa",
            output_dir=out_dir,
            log_dir=tmp_path / "logs",
            backend="docker",
            image="img",
            organism="org",
            threads=1,
            plus=False,
            force=True
        )
        assert res['status'] == 'failed' # Because exit code 1
        assert mock_run.call_count == 1
