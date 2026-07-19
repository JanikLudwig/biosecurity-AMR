import json
from pathlib import Path

from genome_firewall.annotation.batch import _is_current


def test_cached_annotation_requires_matching_versions(tmp_path: Path) -> None:
    output = tmp_path / "sample.tsv"
    output.write_text("header\n", encoding="utf-8")
    output.with_suffix(".provenance.json").write_text(
        json.dumps({"executable_version": "4.2.7", "database_version": "2026-05-15.1"}),
        encoding="utf-8",
    )
    assert _is_current(output, software="4.2.7", database="2026-05-15.1")
    assert not _is_current(output, software="4.2.7", database="different")
