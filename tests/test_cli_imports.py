import pytest
import sys
from unittest.mock import patch
from genome_firewall import audit_bvbrc, prepare_cohort, download_fastas

def test_imports_do_not_execute_logic():
    # If importing them executed logic, it would fail because data is missing.
    # Since we can import them here, it passes.
    assert hasattr(audit_bvbrc, "main")
    assert hasattr(prepare_cohort, "main")
    assert hasattr(download_fastas, "main")

def test_audit_bvbrc_help(capsys):
    with pytest.raises(SystemExit) as excinfo:
        audit_bvbrc.main(["--help"])
    assert excinfo.value.code == 0
    out, err = capsys.readouterr()
    assert "Audit BV-BRC Metadata" in out

def test_prepare_cohort_help(capsys):
    with pytest.raises(SystemExit) as excinfo:
        prepare_cohort.main(["--help"])
    assert excinfo.value.code == 0
    out, err = capsys.readouterr()
    assert "Prepare cohort" in out

def test_download_fastas_help(capsys):
    with pytest.raises(SystemExit) as excinfo:
        # download_fastas in HackNation has a `def main():` without args parameter.
        # Wait, I didn't change download_fastas.py to accept args=None. Let's patch sys.argv.
        with patch.object(sys, 'argv', ['download_fastas.py', '--help']):
            download_fastas.main()
    assert excinfo.value.code == 0
    out, err = capsys.readouterr()
    assert "Download FASTA" in out

def test_audit_bvbrc_missing_data(caplog):
    # Execute without --help, should fail gracefully returning 1
    # We patch sys.argv to just be the script name (no args needed for audit)
    code = audit_bvbrc.main([])
    assert code == 1
    assert "Fehlende Rohdaten" in caplog.text

def test_prepare_cohort_missing_data(caplog):
    code = prepare_cohort.main([])
    assert code == 1
    assert "Fehlende Rohdaten" in caplog.text
