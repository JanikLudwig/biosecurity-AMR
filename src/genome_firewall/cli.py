from __future__ import annotations

import asyncio
import hashlib
from pathlib import Path

import typer
import pandas as pd
from rich.console import Console
from rich.table import Table

from genome_firewall.annotation.batch import annotate_batch, write_batch_outputs
from genome_firewall.annotation.amrfinder import (
    database_version,
    executable_version,
    parse_output,
    resolve_executable,
    run_nucleotide,
)
from genome_firewall.config import DEFAULT_CONFIG, load_config
from genome_firewall.decision import DEFAULT_DRUG_REGISTRY
from genome_firewall.data.bvbrc import download_and_qc
from genome_firewall.data.training_v1 import prepare_training_v1
from genome_firewall.data.phenotypes import (
    audit_source,
    load_and_clean,
    select_genomes,
    summarize_labels,
    write_phenotype_audit,
    write_selection,
)
from genome_firewall.splitting.homology import (
    assign_grouped_splits,
    cluster_fastas,
    phenotype_split_support,
    write_split_outputs,
)
from genome_firewall.validation import evaluate_report_directory
from genome_firewall.modeling.baseline import load_modeling_table, train_all_drugs
from genome_firewall.inference import predict_fasta
from genome_firewall.lineage import build_lineage_reference

app = typer.Typer(help="Defensive S. aureus antimicrobial-resistance research pipeline.")
data_app = typer.Typer(help="Inspect and prepare BV-BRC phenotype data.")
amrfinder_app = typer.Typer(help="Check and run AMRFinderPlus.")
model_app = typer.Typer(help="Train and evaluate calibrated per-antibiotic models.")
app.add_typer(data_app, name="data")
app.add_typer(amrfinder_app, name="amrfinder")
app.add_typer(model_app, name="model")
console = Console()


@app.command("api")
def serve_api(
    host: str = typer.Option("127.0.0.1", "--host"),
    port: int = typer.Option(8000, "--port", min=1, max=65535),
    reload: bool = typer.Option(False, "--reload"),
) -> None:
    """Serve the FastAPI backend and its small demonstration frontend."""
    import uvicorn

    uvicorn.run(
        "genome_firewall.api.app:app",
        host=host,
        port=port,
        reload=reload,
    )


@app.command("predict")
def predict(
    fasta: Path = typer.Argument(..., exists=True, readable=True),
    model_directory: Path = typer.Option(Path("artifacts/models"), "--model-directory"),
    output: Path = typer.Option(Path("artifacts/reports"), "--output"),
    registry: Path = typer.Option(DEFAULT_DRUG_REGISTRY, "--registry"),
    config_path: Path = typer.Option(DEFAULT_CONFIG, "--config"),
    threads: int = typer.Option(2, min=1),
) -> None:
    """Generate a defensive antibiotic-response report from one assembled FASTA."""
    report = predict_fasta(
        fasta,
        config=load_config(config_path),
        model_directory=model_directory,
        output_directory=output,
        registry_path=registry,
        threads=threads,
    )
    frame = pd.DataFrame(report["decisions"])
    console.print(
        frame[["antibiotic", "call", "confidence", "evidence_category", "target_status"]]
    )
    console.print(report["warning"])
    console.print(f"Report: {output / f'{fasta.stem}.report.json'}")


def _clean(config_path: Path):
    config = load_config(config_path)
    dataset = config["dataset"]
    result = load_and_clean(
        Path(dataset["source_csv"]),
        species=dataset["species"],
        taxon_id=dataset["taxon_id"],
        evidence=dataset["evidence"],
        antibiotics=dataset["antibiotics"],
    )
    return config, result


@data_app.command("audit")
def audit(
    config_path: Path = typer.Option(DEFAULT_CONFIG, "--config"),
    output: Path = typer.Option(Path("data/processed/phenotype-audit"), "--output"),
) -> None:
    """Write AST label, measurement, testing-standard, and method provenance reports."""
    config = load_config(config_path)
    dataset = config["dataset"]
    result = audit_source(
        Path(dataset["source_csv"]),
        species=dataset["species"],
        taxon_id=dataset["taxon_id"],
        evidence=dataset["evidence"],
        antibiotics=dataset["antibiotics"],
    )
    write_phenotype_audit(output, result)
    console.print(result.summary)
    console.print(f"Phenotype provenance audit: {output}")


@data_app.command("summarize")
def summarize(config_path: Path = typer.Option(DEFAULT_CONFIG, "--config")) -> None:
    """Print cleaned label counts for the configured species and drugs."""
    config, result = _clean(config_path)
    summary = summarize_labels(result.labels)
    table = Table(title=f"{config['dataset']['species']} laboratory AST labels")
    for column in ["antibiotic", "Resistant", "Susceptible", "total", "resistant_fraction"]:
        table.add_column(column)
    for row in summary.to_dict(orient="records"):
        table.add_row(
            str(row["antibiotic"]),
            str(row["Resistant"]),
            str(row["Susceptible"]),
            str(row["total"]),
            f"{row['resistant_fraction']:.1%}",
        )
    console.print(table)
    console.print(f"Usable genomes: {result.labels['genome_id'].nunique():,}")
    console.print(f"Excluded: {result.excluded_counts}")


@data_app.command("select")
def select(
    config_path: Path = typer.Option(DEFAULT_CONFIG, "--config"),
    output: Path = typer.Option(Path("data/manifests/initial"), "--output"),
) -> None:
    """Build a deterministic initial genome and phenotype manifest."""
    config, result = _clean(config_path)
    dataset = config["dataset"]
    manifest, labels = select_genomes(
        result.labels,
        max_genomes=dataset["max_genomes"],
        minimum_per_class=dataset["minimum_per_class"],
        seed=dataset["seed"],
    )
    source = Path(dataset["source_csv"])
    metadata = {
        "species": dataset["species"],
        "taxon_id": dataset["taxon_id"],
        "antibiotics": dataset["antibiotics"],
        "genome_count": len(manifest),
        "label_count": len(labels),
        "seed": dataset["seed"],
        "source": str(source),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "excluded_counts": result.excluded_counts,
    }
    write_selection(output, manifest, labels, metadata)
    console.print(f"Wrote {len(manifest):,} genomes and {len(labels):,} labels to {output}")


@data_app.command("prepare-training-v1")
def prepare_v1(
    source: Path = typer.Option(Path("data/training-data-v1"), "--source"),
    features: Path = typer.Option(
        Path("data/processed/training-v1/features.parquet"), "--features"
    ),
    phenotypes: Path = typer.Option(
        Path("data/processed/training-v1/phenotypes.csv"), "--phenotypes"
    ),
    genomes: Path = typer.Option(
        Path("data/processed/training-v1/genomes.csv"), "--genomes"
    ),
    manifest: Path = typer.Option(
        Path("data/processed/training-v1/dataset-manifest.json"), "--manifest"
    ),
) -> None:
    """Normalize the supplied 3k tables without rerunning AMRFinderPlus."""
    summary = prepare_training_v1(
        source,
        features_output=features,
        phenotypes_output=phenotypes,
        genomes_output=genomes,
        manifest_output=manifest,
    )
    console.print(
        f"Prepared {summary['genomes']:,} genomes, {summary['labels']:,} labels, "
        f"and {summary['features']:,} features"
    )
    console.print(f"Download manifest: {genomes}")


@data_app.command("download")
def download(
    config_path: Path = typer.Option(DEFAULT_CONFIG, "--config"),
    manifest: Path = typer.Option(Path("data/manifests/initial/genomes.csv"), "--manifest"),
    output: Path = typer.Option(Path("data/raw/genomes"), "--output"),
    qc_output: Path = typer.Option(Path("data/manifests/qc.csv"), "--qc-output"),
    limit: int | None = typer.Option(None, min=1, help="Download only N genomes for a smoke run."),
    sample_seed: int | None = typer.Option(
        None, help="Deterministically sample across the manifest instead of taking its first N rows."
    ),
) -> None:
    """Download selected BV-BRC FASTAs and validate assembly quality/provenance."""
    config = load_config(config_path)
    dataset = config["dataset"]
    qc = asyncio.run(
        download_and_qc(
            manifest,
            output,
            qc_output,
            species=dataset["species"],
            taxon_id=dataset["taxon_id"],
            quality=config["quality"],
            bvbrc=config["bvbrc"],
            limit=limit,
            sample_seed=sample_seed,
        )
    )
    passed = int(qc["passed_qc"].sum())
    console.print(f"Downloaded/checked {len(qc):,} genomes: {passed:,} passed QC")
    console.print(f"QC manifest: {qc_output}")


@data_app.command("cluster-split")
def cluster_split(
    fasta_directory: Path = typer.Option(Path("data/raw/genomes"), "--fasta-directory"),
    qc_manifest: Path = typer.Option(Path("data/manifests/qc-dev.csv"), "--qc-manifest"),
    output: Path = typer.Option(Path("data/processed/splits"), "--output"),
    config_path: Path = typer.Option(DEFAULT_CONFIG, "--config"),
    phenotypes: Path = typer.Option(
        Path("data/manifests/initial/phenotypes.csv"), "--phenotypes"
    ),
) -> None:
    """Cluster related genomes and assign leakage-safe dataset splits."""
    config = load_config(config_path)
    similarity = config["similarity"]
    qc = pd.read_csv(qc_manifest, dtype=object)
    passed_ids = set(
        qc.loc[qc["passed_qc"].astype(str).str.casefold().eq("true"), "genome_id"]
    )
    fasta_paths = sorted(
        path for path in fasta_directory.glob("*.fna") if path.stem in passed_ids
    )
    if not fasta_paths:
        console.print("[red]No QC-passing FASTA files were found.[/red]")
        raise typer.Exit(1)

    membership, edges = cluster_fastas(
        fasta_paths,
        ksize=similarity["ksize"],
        scaled=similarity["scaled"],
        ani_threshold=similarity["ani_threshold"],
    )
    splits = assign_grouped_splits(
        membership,
        train_fraction=similarity["train_fraction"],
        calibration_fraction=similarity["calibration_fraction"],
        test_fraction=similarity["test_fraction"],
        seed=similarity["seed"],
    )
    write_split_outputs(output, splits, edges, similarity)
    labels = pd.read_csv(phenotypes, dtype=object, keep_default_na=False)
    support = phenotype_split_support(splits, labels)
    support.to_csv(output / "phenotype-support.csv", index=False)
    summary = splits.groupby("split").agg(
        genomes=("genome_id", "size"), clusters=("cluster_id", "nunique")
    )
    console.print(summary)
    console.print(support)
    console.print(f"Near-duplicate edges at threshold: {len(edges):,}")
    console.print(f"Split outputs: {output}")


@data_app.command("split-support")
def split_support(
    splits: Path = typer.Option(
        Path("data/processed/splits/genome-splits.csv"), "--splits"
    ),
    phenotypes: Path = typer.Option(
        Path("data/manifests/initial/phenotypes.csv"), "--phenotypes"
    ),
    output: Path | None = typer.Option(None, "--output"),
) -> None:
    """Report per-drug phenotype class support for an existing grouped split."""
    split_frame = pd.read_csv(splits, dtype=object, keep_default_na=False)
    labels = pd.read_csv(phenotypes, dtype=object, keep_default_na=False)
    support = phenotype_split_support(split_frame, labels)
    destination = output or splits.parent / "phenotype-support.csv"
    destination.parent.mkdir(parents=True, exist_ok=True)
    support.to_csv(destination, index=False)
    console.print(support)
    console.print(f"Phenotype support: {destination}")


@amrfinder_app.command("doctor")
def amrfinder_doctor(
    config_path: Path = typer.Option(DEFAULT_CONFIG, "--config"),
) -> None:
    """Confirm that the configured AMRFinderPlus executable is usable."""
    config = load_config(config_path)
    executable = resolve_executable(config["amrfinder"]["executable"])
    if executable is None:
        console.print("[red]AMRFinderPlus was not found.[/red]")
        console.print("Run ./scripts/setup-amrfinder.sh after installing micromamba/mamba/conda.")
        raise typer.Exit(1)
    console.print(f"Executable: {executable}")
    console.print(f"Version: {executable_version(executable)}")
    installed_database = database_version(executable)
    if installed_database is None:
        console.print("[red]AMRFinderPlus database was not found.[/red]")
        console.print("Run ./scripts/setup-amrfinder.sh to download it.")
        raise typer.Exit(1)
    console.print(f"Database: {installed_database}")


@amrfinder_app.command("run")
def amrfinder_run(
    fasta: Path = typer.Argument(..., exists=True, readable=True),
    genome_id: str | None = typer.Option(None, help="Defaults to the FASTA filename stem."),
    output: Path | None = typer.Option(None, help="Raw AMRFinder TSV destination."),
    config_path: Path = typer.Option(DEFAULT_CONFIG, "--config"),
    threads: int = typer.Option(2, min=1),
) -> None:
    """Annotate one assembled nucleotide FASTA and write raw + parsed evidence."""
    config = load_config(config_path)
    executable = resolve_executable(config["amrfinder"]["executable"])
    if executable is None or database_version(executable) is None:
        console.print("[red]AMRFinderPlus executable/database is not ready.[/red]")
        console.print("Run ./scripts/setup-amrfinder.sh, then rerun this command.")
        raise typer.Exit(1)

    resolved_genome_id = genome_id or fasta.stem
    raw_output = output or Path("data/interim/amrfinder") / f"{resolved_genome_id}.tsv"
    run_nucleotide(
        executable,
        fasta,
        raw_output,
        organism=config["amrfinder"]["organism"],
        threads=threads,
    )
    evidence = parse_output(raw_output, genome_id=resolved_genome_id)
    parsed_output = raw_output.with_suffix(".evidence.csv")
    evidence.to_csv(parsed_output, index=False)
    console.print(f"Detected {len(evidence)} AMR elements")
    console.print(f"Raw output: {raw_output}")
    console.print(f"Parsed evidence: {parsed_output}")


@amrfinder_app.command("batch")
def amrfinder_batch(
    fasta_directory: Path = typer.Option(Path("data/raw/genomes"), "--fasta-directory"),
    output: Path = typer.Option(Path("data/processed/features"), "--output"),
    raw_output: Path = typer.Option(Path("data/interim/amrfinder"), "--raw-output"),
    qc_manifest: Path | None = typer.Option(
        None, "--qc-manifest", help="If provided, annotate only genomes that passed QC."
    ),
    limit: int | None = typer.Option(None, min=1),
    workers: int = typer.Option(2, min=1),
    threads_per_worker: int = typer.Option(1, min=1),
    config_path: Path = typer.Option(DEFAULT_CONFIG, "--config"),
) -> None:
    """Annotate a resumable genome batch and build AMR feature tables."""
    config = load_config(config_path)
    executable = resolve_executable(config["amrfinder"]["executable"])
    if executable is None or database_version(executable) is None:
        console.print("[red]AMRFinderPlus executable/database is not ready.[/red]")
        raise typer.Exit(1)

    fasta_paths = sorted(fasta_directory.glob("*.fna"))
    if qc_manifest is not None:
        qc = pd.read_csv(qc_manifest, dtype=object)
        passed_ids = set(qc.loc[qc["passed_qc"].astype(str).str.casefold().eq("true"), "genome_id"])
        fasta_paths = [path for path in fasta_paths if path.stem in passed_ids]
    if limit is not None:
        fasta_paths = fasta_paths[:limit]
    if not fasta_paths:
        console.print("[red]No eligible FASTA files were found.[/red]")
        raise typer.Exit(1)

    status, evidence, features = annotate_batch(
        fasta_paths,
        raw_output,
        executable=executable,
        organism=config["amrfinder"]["organism"],
        workers=workers,
        threads_per_worker=threads_per_worker,
    )
    write_batch_outputs(output, status, evidence, features)
    failures = int(status["status"].eq("failed").sum())
    console.print(
        f"Processed {len(status):,} genomes: {len(evidence):,} evidence rows, "
        f"{features.shape[1]:,} features, {failures:,} failures"
    )
    console.print(f"Batch outputs: {output}")


@model_app.command("train")
def model_train(
    features: Path = typer.Option(
        Path("data/processed/features/amr-features.parquet"), "--features"
    ),
    phenotypes: Path = typer.Option(
        Path("data/manifests/initial/phenotypes.csv"), "--phenotypes"
    ),
    splits: Path = typer.Option(
        Path("data/processed/splits/genome-splits.csv"), "--splits"
    ),
    output: Path = typer.Option(Path("artifacts/models"), "--output"),
    antibiotics: list[str] | None = typer.Option(
        None, "--antibiotic", help="Train only this antibiotic; repeat the option for several."
    ),
    evaluation_status: str = typer.Option(
        "grouped-development", "--evaluation-status"
    ),
    config_path: Path = typer.Option(DEFAULT_CONFIG, "--config"),
) -> None:
    """Train one baseline model per drug using fixed grouped splits."""
    config = load_config(config_path)
    table, feature_columns = load_modeling_table(features, phenotypes, splits)
    summary = train_all_drugs(
        table,
        feature_columns,
        antibiotics=antibiotics or config["dataset"]["antibiotics"],
        config=config["model"],
        output_directory=output,
        bundle_metadata={
            "species": config["dataset"]["species"],
            "evaluation_status": evaluation_status,
            "split_genomes": int(
                pd.read_csv(splits, dtype=object, usecols=["genome_id"])["genome_id"].nunique()
            ),
            "split_path": str(splits),
            "features_path": str(features),
            "phenotypes_path": str(phenotypes),
        },
    )
    console.print(summary[["antibiotic", "status", "calibration_status"]])
    console.print(f"Model artifacts: {output}")


@model_app.command("lineage")
def model_lineage(
    splits: Path = typer.Option(
        Path("data/processed/splits-500/genome-splits.csv"), "--splits"
    ),
    fasta_directory: Path = typer.Option(Path("data/raw/genomes"), "--fasta-directory"),
    output: Path = typer.Option(
        Path("artifacts/models/lineage-reference.joblib"), "--output"
    ),
    config_path: Path = typer.Option(DEFAULT_CONFIG, "--config"),
) -> None:
    """Build the calibrated sequence-lineage novelty reference used at inference."""
    config = load_config(config_path)
    artifact = build_lineage_reference(
        splits_path=splits,
        fasta_directory=fasta_directory,
        output_path=output,
        ksize=config["similarity"]["ksize"],
        scaled=config["similarity"]["scaled"],
        calibration_quantile=config["lineage"]["calibration_quantile"],
    )
    console.print(
        f"Lineage reference: {len(artifact['training_genome_ids']):,} training genomes, "
        f"minimum estimated ANI {artifact['minimum_training_ani']:.4%}"
    )
    console.print(f"Artifact: {output}")


@model_app.command("evaluate-reports")
def model_evaluate_reports(
    reports: Path = typer.Option(Path("artifacts/reports/external"), "--reports"),
    phenotypes: Path = typer.Option(..., "--phenotypes"),
    output: Path = typer.Option(Path("artifacts/external-validation"), "--output"),
) -> None:
    """Evaluate frozen decision reports against an external laboratory cohort."""
    summary, matched = evaluate_report_directory(reports, phenotypes)
    output.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output / "summary.csv", index=False)
    matched.to_csv(output / "matched-decisions.csv", index=False)
    console.print(summary)
    console.print(f"External validation outputs: {output}")
