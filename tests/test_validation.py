import json
from pathlib import Path

import pandas as pd

from genome_firewall.validation import evaluate_report_directory


def test_external_report_evaluation_counts_no_calls(tmp_path: Path) -> None:
    reports = tmp_path / "reports"
    reports.mkdir()
    (reports / "g1.report.json").write_text(
        json.dumps(
            {
                "genome_id": "g1",
                "decisions": [
                    {"antibiotic": "drug-a", "call": "likely_to_fail", "confidence": 0.9},
                    {"antibiotic": "drug-b", "call": "no_call", "confidence": 0.6},
                ],
            }
        ),
        encoding="utf-8",
    )
    phenotypes = tmp_path / "phenotypes.csv"
    pd.DataFrame(
        {
            "genome_id": ["g1", "g1"],
            "antibiotic": ["drug-a", "drug-b"],
            "label": ["Resistant", "Susceptible"],
        }
    ).to_csv(phenotypes, index=False)

    summary, matched = evaluate_report_directory(reports, phenotypes)
    by_drug = summary.set_index("antibiotic")
    assert by_drug.loc["drug-a", "called_accuracy"] == 1.0
    assert by_drug.loc["drug-b", "no_call_rate"] == 1.0
    assert len(matched) == 2
