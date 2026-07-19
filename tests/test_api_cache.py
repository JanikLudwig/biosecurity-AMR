from __future__ import annotations

import json
from pathlib import Path

from genome_firewall.api.app import _load_cached_report


def test_cached_report_requires_the_same_fasta_digest(tmp_path: Path) -> None:
    directory = tmp_path / "analysis"
    directory.mkdir()
    report = {"qc": {"sha256": "expected"}, "decisions": [{"call": "no_call"}]}
    (directory / "analysis.report.json").write_text(json.dumps(report), encoding="utf-8")

    assert _load_cached_report(directory, "expected") == report
    assert _load_cached_report(directory, "different") is None
