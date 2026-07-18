"""Command-line interface for Genome Firewall v0.

Examples
--------
    # Precomputed AMRFinderPlus table (runs anywhere, no bioinformatics install):
    python -m genome_firewall.cli predict \
        --tsv genome_firewall/examples/ecoli_resistant_amrfinder.tsv \
        --tsv-source amrfinderplus --species "Escherichia coli"

    # From an assembly, letting the tool pick a backend (AMRFinderPlus/cAMRah):
    python -m genome_firewall.cli predict --fasta assembly.fasta --annotator auto

    # De-duplicate a genome collection for leakage-free evaluation:
    python -m genome_firewall.cli dedup a.fasta b.fasta c.fasta --threshold 0.9
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import List, Optional

from . import __version__
from .annotate import annotate
from .fasta import compute_qc
from .predict import predict_sample
from .report import render_markdown, render_text


def _cmd_predict(args) -> int:
    qc = None
    if args.fasta:
        qc = compute_qc(args.fasta).as_dict()

    annotation = annotate(
        fasta_path=args.fasta,
        backend=args.annotator,
        tsv_path=args.tsv,
        organism=args.organism,
        tsv_source=args.tsv_source,
    )
    sample = predict_sample(annotation, species=args.species, qc=qc)

    if args.format == "json":
        print(json.dumps(sample.as_dict(), indent=2))
    elif args.format == "md":
        print(render_markdown(sample))
    else:
        print(render_text(sample))
    return 0


def _cmd_dedup(args) -> int:
    from .dedup import cluster_genomes

    clusters = cluster_genomes(args.fasta, k=args.k, num_hashes=args.num_hashes,
                               threshold=args.threshold)
    print(f"{len(args.fasta)} genomes -> {len(clusters)} clusters "
          f"(Jaccard threshold {args.threshold})")
    for i, c in enumerate(clusters, 1):
        print(f"  cluster {i}: representative={c.representative} "
              f"({len(c.members)} member(s))")
        for m in c.members:
            if m != c.representative:
                print(f"      - {m}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="genome-firewall",
        description="Genome Firewall v0 — zero-shot, rule-based AMR decision support "
                    "(research prototype; confirm every result with lab testing).",
    )
    p.add_argument("--version", action="version", version=f"genome-firewall {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    pr = sub.add_parser("predict", help="Predict antibiotic response for one genome.")
    src = pr.add_argument_group("input (provide --fasta and/or --tsv)")
    src.add_argument("--fasta", help="Reconstructed assembly (FASTA, optionally .gz).")
    src.add_argument("--tsv", help="Precomputed AMR table (AMRFinderPlus/cAMRah/Abricate).")
    pr.add_argument("--annotator", default="auto",
                    choices=["auto", "amrfinderplus", "camrah", "tsv"],
                    help="Annotation backend (default: auto).")
    pr.add_argument("--tsv-source", default="amrfinderplus",
                    choices=["amrfinderplus", "camrah", "abricate", "resfinder", "tsv"],
                    help="Tool that produced --tsv (affects screening-completeness).")
    pr.add_argument("--species", default=None,
                    help="Bacterial species label (e.g. 'Escherichia coli').")
    pr.add_argument("--organism", default=None,
                    help="AMRFinderPlus --organism value (e.g. Escherichia).")
    pr.add_argument("--format", default="text", choices=["text", "md", "json"])
    pr.set_defaults(func=_cmd_predict)

    dd = sub.add_parser("dedup", help="Cluster genomes by sequence homology (evaluation).")
    dd.add_argument("fasta", nargs="+", help="Genome FASTA files to cluster.")
    dd.add_argument("--k", type=int, default=21, help="k-mer size (default 21).")
    dd.add_argument("--num-hashes", type=int, default=200, help="MinHash sketch size.")
    dd.add_argument("--threshold", type=float, default=0.9,
                    help="Jaccard threshold to merge near-identical genomes.")
    dd.set_defaults(func=_cmd_dedup)
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "predict" and not args.fasta and not args.tsv:
        print("error: predict needs --fasta and/or --tsv", file=sys.stderr)
        return 2
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
