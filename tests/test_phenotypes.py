from pathlib import Path

import pandas as pd

from genome_firewall.data.phenotypes import audit_source, summarize_labels


def test_summary_reports_both_classes() -> None:
    labels = pd.DataFrame(
        {
            "genome_id": ["a", "b", "c"],
            "antibiotic": ["drug", "drug", "drug"],
            "label": ["Resistant", "Susceptible", "Susceptible"],
        }
    )
    row = summarize_labels(labels).iloc[0]
    assert row["Resistant"] == 1
    assert row["Susceptible"] == 2


def test_audit_preserves_measurement_and_standard_provenance(tmp_path: Path) -> None:
    source = tmp_path / "ast.csv"
    pd.DataFrame(
        {
            "Taxon ID": ["1280", "1280"],
            "Genome ID": ["a", "b"],
            "Genome Name": ["Staphylococcus aureus A", "Staphylococcus aureus B"],
            "Antibiotic": ["Ciprofloxacin", "Ciprofloxacin"],
            "Resistant Phenotype": ["Susceptible", "Intermediate"],
            "Evidence": ["Laboratory Method", "Laboratory Method"],
            "Measurement": ["<=0.5", ""],
            "Testing Standard": ["EUCAST", ""],
            "Testing Standard Year": ["2025", ""],
            "Laboratory Typing Method": ["Broth dilution", "Disk diffusion"],
        }
    ).to_csv(source, index=False)

    audit = audit_source(
        source,
        species="Staphylococcus aureus",
        taxon_id=1280,
        evidence="Laboratory Method",
        antibiotics=["ciprofloxacin"],
    )

    row = audit.summary.iloc[0]
    assert row["observation_rows"] == 2
    assert row["measurement_rows"] == 1
    assert row["rows_with_standard"] == 1
    assert set(audit.standards["testing_standard"]) == {"EUCAST", "missing"}
