#!/usr/bin/env python3
"""Explore how an AMRFinderPlus report becomes Genome Firewall model features.

This is an educational sandbox, not a predictor. It displays:

1. AMRFinderPlus's raw TSV rows.
2. Genome Firewall's normalized gene/mutation evidence.
3. A one-genome binary feature vector.

Examples:
    uv run python scripts/explore_amrfinder.py
    uv run python scripts/explore_amrfinder.py data/raw/genomes/1280.9342.fna --force
    uv run python scripts/explore_amrfinder.py --show-all
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from genome_firewall.annotation.amrfinder import (
    database_version,
    executable_version,
    parse_output,
    resolve_executable,
    run_nucleotide,
)
from genome_firewall.config import DEFAULT_CONFIG, load_config

DEFAULT_FASTA = Path("data/raw/genomes/1280.9342.fna")
DEFAULT_RAW_DIRECTORY = Path("data/interim/amrfinder")
DEFAULT_EXPLORER_DIRECTORY = Path("data/processed/amrfinder-explorer")


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run or reuse AMRFinderPlus and inspect how its rows become model features."
    )
    parser.add_argument(
        "fasta",
        nargs="?",
        type=Path,
        default=DEFAULT_FASTA,
        help=f"Assembled nucleotide FASTA (default: {DEFAULT_FASTA})",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="Genome Firewall TOML configuration.",
    )
    parser.add_argument(
        "--raw-directory",
        type=Path,
        default=DEFAULT_RAW_DIRECTORY,
        help="Directory containing raw AMRFinder TSV files.",
    )
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=DEFAULT_EXPLORER_DIRECTORY,
        help="Where to write the normalized evidence and sample feature vector.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rerun AMRFinderPlus even when a raw TSV already exists.",
    )
    parser.add_argument(
        "--show-all",
        action="store_true",
        help="Print every result row instead of the first 30.",
    )
    return parser.parse_args()


def rich_table(frame: pd.DataFrame, columns: list[str], *, limit: int | None) -> Table:
    visible = frame.loc[:, columns]
    if limit is not None:
        visible = visible.head(limit)
    table = Table(show_lines=False)
    for column in columns:
        table.add_column(column, overflow="fold", max_width=38)
    for row in visible.itertuples(index=False, name=None):
        table.add_row(*(str(value) for value in row))
    return table


def feature_vector(evidence: pd.DataFrame, *, genome_id: str) -> pd.DataFrame:
    """Build the one-row representation used as input to the baseline models."""
    if evidence.empty:
        return pd.DataFrame({"genome_id": [genome_id]})
    vector = (
        evidence.pivot_table(
            index="genome_id",
            columns="feature_key",
            values="feature_value",
            aggfunc="max",
            fill_value=0,
        )
        .astype("uint8")
        .reset_index()
    )
    vector.columns.name = None
    return vector


def main() -> int:
    args = arguments()
    console = Console()
    config = load_config(args.config)

    if not args.fasta.is_file():
        console.print(f"[red]FASTA does not exist: {args.fasta}[/red]")
        return 1
    executable = resolve_executable(config["amrfinder"]["executable"])
    if executable is None or database_version(executable) is None:
        console.print("[red]AMRFinderPlus executable/database is not ready.[/red]")
        console.print("Run ./scripts/setup-amrfinder.sh first.")
        return 1

    genome_id = args.fasta.stem
    raw_output = args.raw_directory / f"{genome_id}.tsv"
    if args.force or not raw_output.is_file():
        console.print(f"Running AMRFinderPlus on [cyan]{args.fasta}[/cyan] ...")
        run_nucleotide(
            executable,
            args.fasta,
            raw_output,
            organism=config["amrfinder"]["organism"],
            threads=2,
        )
    else:
        console.print(f"Reusing existing raw report: [cyan]{raw_output}[/cyan]")

    raw = pd.read_csv(raw_output, sep="\t", dtype=object, keep_default_na=False)
    evidence = parse_output(raw_output, genome_id=genome_id)
    vector = feature_vector(evidence, genome_id=genome_id)

    args.output_directory.mkdir(parents=True, exist_ok=True)
    evidence_path = args.output_directory / f"{genome_id}.evidence.csv"
    vector_path = args.output_directory / f"{genome_id}.feature-vector.csv"
    evidence.to_csv(evidence_path, index=False)
    vector.to_csv(vector_path, index=False)

    console.print(
        Panel.fit(
            f"Genome: [bold]{genome_id}[/bold]\n"
            f"AMRFinderPlus: {executable_version(executable)}\n"
            f"Database: {database_version(executable)}\n"
            f"Raw elements: {len(raw)}\n"
            f"Binary features present: {len(vector.columns) - 1}",
            title="Run provenance",
        )
    )

    limit = None if args.show_all else 30
    raw_columns = [
        "Element symbol",
        "Element name",
        "Type",
        "Subtype",
        "Class",
        "Subclass",
        "Method",
        "% Coverage of reference",
        "% Identity to reference",
    ]
    console.rule("1. Raw AMRFinderPlus output")
    console.print(rich_table(raw, raw_columns, limit=limit))
    if limit is not None and len(raw) > limit:
        console.print(f"[dim]Showing {limit} of {len(raw)} rows; pass --show-all for all rows.[/dim]")

    console.rule("2. Normalized Genome Firewall evidence")
    evidence_columns = [
        "feature_key",
        "element_symbol",
        "evidence_category",
        "amr_class",
        "amr_subclass",
        "method",
        "coverage",
        "identity",
        "feature_value",
    ]
    console.print(rich_table(evidence, evidence_columns, limit=limit))

    console.rule("3. Binary feature vector")
    present = [column for column in vector.columns if column != "genome_id"]
    feature_rows = pd.DataFrame({"feature": present, "value": 1})
    console.print(rich_table(feature_rows, ["feature", "value"], limit=limit))

    console.rule("How the transformation currently works")
    console.print(
        "[bold]Known mutation:[/bold] AMRFinder Subtype=POINT or Method starts with POINT "
        "becomes [cyan]mutation::<symbol>[/cyan].\n"
        "[bold]Known gene/element:[/bold] every other reported row becomes "
        "[cyan]gene::<symbol>[/cyan].\n"
        "[bold]Presence:[/bold] a detected element receives value 1. In a cohort matrix, "
        "elements absent from a genome receive value 0.\n\n"
        "[yellow]Important:[/yellow] this stage performs deterministic feature construction, "
        "not statistical feature selection. The cohort matrix uses the union of all detected "
        "feature keys. Regularized logistic regression later shrinks less useful coefficients. "
        "Coverage, identity, and AMRFinder method remain in the evidence table for explanation "
        "and future filtering, but are not currently numeric model inputs."
    )
    console.print(f"\nNormalized evidence written to: [cyan]{evidence_path}[/cyan]")
    console.print(f"Feature vector written to: [cyan]{vector_path}[/cyan]")
    console.print(
        "[bold red]Research prototype only:[/bold red] this script does not make an antibiotic "
        "treatment prediction."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

