from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd

from genome_firewall.annotation.amrfinder import (
    database_version,
    executable_version,
    parse_output,
    resolve_executable,
    run_nucleotide,
)
from genome_firewall.data.qc import evaluate_quality, inspect_fasta
from genome_firewall.decision import load_drug_registry, predict_antibiotic
from genome_firewall.lineage import evaluate_lineage
from genome_firewall.modeling.bundle import ModelBundle
from genome_firewall.targets import analyze_drug_targets

CONFIRMATION_WARNING = (
    "Research prototype only. Every prediction must be confirmed with standard "
    "laboratory antimicrobial-susceptibility testing."
)


def predict_fasta(
    fasta: Path,
    *,
    config: dict[str, Any],
    model_directory: Path,
    output_directory: Path,
    registry_path: Path,
    lineage_artifact: Path | None = None,
    model_bundle: ModelBundle | None = None,
    threads: int = 2,
) -> dict[str, Any]:
    """Run the defensive FASTA-to-report path for one assembled S. aureus genome."""
    executable = resolve_executable(config["amrfinder"]["executable"])
    if executable is None or database_version(executable) is None:
        raise RuntimeError("AMRFinderPlus executable/database is not ready")
    registry = load_drug_registry(registry_path)
    expected_species = registry["registry"]["species"]
    if expected_species != config["dataset"]["species"]:
        raise ValueError("Drug registry species does not match experiment configuration")

    bundle = model_bundle or ModelBundle(model_directory)
    unsupported = sorted(set(bundle.antibiotics).difference(registry["drugs"]))
    if unsupported:
        raise ValueError(f"Model bundle drugs are absent from the registry: {unsupported}")
    drugs = {antibiotic: registry["drugs"][antibiotic] for antibiotic in bundle.antibiotics}

    genome_id = fasta.stem
    metrics = inspect_fasta(fasta)
    assembly_quality = dict(config["quality"])
    assembly_quality["allowed_genome_quality"] = []
    qc_reasons = evaluate_quality(metrics, {}, assembly_quality)
    qc_passed = not qc_reasons
    lineage_status = evaluate_lineage(
        fasta, lineage_artifact or model_directory / "lineage-reference.joblib"
    )

    output_directory.mkdir(parents=True, exist_ok=True)
    raw_amr = output_directory / f"{genome_id}.amrfinder.tsv"
    if qc_passed:
        run_nucleotide(
            executable,
            fasta,
            raw_amr,
            organism=config["amrfinder"]["organism"],
            threads=threads,
        )
        evidence = parse_output(raw_amr, genome_id=genome_id)
        target_analysis = analyze_drug_targets(
            fasta,
            amrfinder_executable=executable,
            drugs=drugs,
        )
        target_statuses = target_analysis["drugs"]
    else:
        evidence = pd.DataFrame(
            columns=[
                "genome_id",
                "feature_key",
                "element_symbol",
                "evidence_category",
                "amr_subclass",
                "coverage",
                "identity",
            ]
        )
        target_statuses = {
            antibiotic: {"status": "not_evaluated", "detected": []}
            for antibiotic in drugs
        }
        target_analysis = {
            "workflow": "M2",
            "status": "not_evaluated",
            "reason": "assembly_qc_failed",
            "drugs": target_statuses,
        }

    decisions = []
    for antibiotic, drug in drugs.items():
        decisions.append(
            predict_antibiotic(
                antibiotic=antibiotic,
                evidence=evidence,
                model_artifact=bundle.artifact(antibiotic),
                drug=drug,
                target_status=target_statuses[antibiotic],
                lineage_status=lineage_status,
                qc_passed=qc_passed,
                qc_reasons=qc_reasons,
            )
        )

    report = {
        "schema_version": "1.0",
        "generated_at": datetime.now(UTC).isoformat(),
        "genome_id": genome_id,
        "species_scope": expected_species,
        "warning": CONFIRMATION_WARNING,
        "defensive_use_only": True,
        "qc": {"passed": qc_passed, "reasons": qc_reasons, **metrics.as_dict()},
        "provenance": {
            "amrfinder_version": executable_version(executable),
            "amrfinder_database": database_version(executable),
            "drug_registry_version": registry["registry"]["version"],
            "model_directory": str(model_directory),
        },
        "evidence_sources": {
            "amrfinder_mapping": registry["registry"]["amrfinder_mapping_source"],
            "drug_mechanisms": registry["registry"]["mechanism_source"],
        },
        "workflows": {
            "M1": {
                "name": "AMRFinderPlus feature extraction",
                "feature_schema": bundle.manifest["feature_schema_sha256"],
                "recognized_features": sorted(
                    set(evidence["feature_key"]).intersection(bundle.feature_columns)
                ),
                "unknown_features": sorted(
                    set(evidence["feature_key"]).difference(bundle.feature_columns)
                ),
                "evidence": json.loads(evidence.to_json(orient="records")),
            },
            "M2": target_analysis,
        },
        "lineage": lineage_status,
        "decisions": decisions,
    }
    evidence.to_csv(output_directory / f"{genome_id}.evidence.csv", index=False)
    pd.DataFrame(decisions).to_csv(output_directory / f"{genome_id}.report.csv", index=False)
    (output_directory / f"{genome_id}.report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return report
