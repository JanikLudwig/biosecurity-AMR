from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


def evaluate_report_directory(
    report_directory: Path, phenotype_path: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Score frozen decision reports against external laboratory labels."""
    rows = []
    for path in sorted(report_directory.rglob("*.report.json")):
        report = json.loads(path.read_text(encoding="utf-8"))
        for decision in report["decisions"]:
            rows.append(
                {
                    "genome_id": report["genome_id"],
                    "antibiotic": decision["antibiotic"],
                    "call": decision["call"],
                    "confidence": decision["confidence"],
                }
            )
    if not rows:
        raise ValueError(f"No decision reports found under {report_directory}")
    reports = pd.DataFrame(rows)
    labels = pd.read_csv(phenotype_path, dtype=object, keep_default_na=False)
    required = {"genome_id", "antibiotic", "label"}
    missing = required.difference(labels.columns)
    if missing:
        raise ValueError(f"External phenotype file is missing columns: {sorted(missing)}")
    matched = reports.merge(
        labels[["genome_id", "antibiotic", "label"]],
        on=["genome_id", "antibiotic"],
        how="inner",
        validate="one_to_one",
    )
    if matched.empty:
        raise ValueError("No report decisions matched external phenotype labels")
    matched["called"] = matched["call"].ne("no_call")
    matched["correct"] = (
        matched["call"].eq("likely_to_fail") & matched["label"].eq("Resistant")
    ) | (matched["call"].eq("likely_to_work") & matched["label"].eq("Susceptible"))

    summaries = []
    for antibiotic, group in matched.groupby("antibiotic"):
        called = group.loc[group["called"]]
        resistant = group.loc[group["label"].eq("Resistant")]
        susceptible = group.loc[group["label"].eq("Susceptible")]
        summaries.append(
            {
                "antibiotic": antibiotic,
                "total": len(group),
                "called": len(called),
                "coverage": len(called) / len(group),
                "no_call_rate": group["call"].eq("no_call").mean(),
                "called_accuracy": called["correct"].mean() if len(called) else None,
                "resistant_recall_including_no_calls": (
                    resistant["call"].eq("likely_to_fail").mean()
                    if len(resistant)
                    else None
                ),
                "susceptible_recall_including_no_calls": (
                    susceptible["call"].eq("likely_to_work").mean()
                    if len(susceptible)
                    else None
                ),
            }
        )
    return pd.DataFrame(summaries).sort_values("antibiotic"), matched
