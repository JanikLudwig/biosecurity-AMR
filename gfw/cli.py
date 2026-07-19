"""Command-line entry point.

    python -m gfw.cli predict --genome 1280.10000            # markdown report
    python -m gfw.cli predict --genome 1280.10000 --json     # JSON
    python -m gfw.cli predict --fasta path/to/assembly.fna --genome-id my_sample
"""

from __future__ import annotations

import argparse
import json
import sys

from .engine import Engine


def _cmd_predict(args) -> int:
    engine = Engine()
    gid = args.genome or args.genome_id or "sample"
    report = engine.predict_genome(gid, fasta_path=args.fasta)
    if args.json:
        print(json.dumps(report.as_dict(), indent=2))
    else:
        print(report.to_markdown())
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="gfw", description="Genome Firewall — S. aureus")
    sub = p.add_subparsers(dest="cmd", required=True)

    pp = sub.add_parser("predict", help="predict antibiotic response for one genome")
    pp.add_argument("--genome", help="genome_id present in genomes/")
    pp.add_argument("--fasta", help="path to an assembly FASTA (instead of --genome)")
    pp.add_argument("--genome-id", help="label to use when --fasta is given")
    pp.add_argument("--json", action="store_true", help="emit JSON instead of Markdown")
    pp.set_defaults(func=_cmd_predict)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
