#!/usr/bin/env python3
"""Precompute full pipeline reports for a subset of genomes (for the UI / demo).

Runs the real end-to-end Engine (M2 target detection included) over a sample of
genomes — by default from the grouped hidden ``test`` partition — and writes one
JSON report per genome plus an ``index.json`` summary. This doubles as the
plan's "verify end-to-end on a subset" step.

    python scripts/precompute_reports.py --n 150 --partition test
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gfw.config import REPORTS_DIR, ensure_dirs
from gfw.engine import Engine
from gfw.split import load_split, make_split


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=150)
    ap.add_argument("--partition", default="test",
                    choices=["test", "train", "calibration", "all"])
    ap.add_argument("--seed", type=int, default=1280)
    args = ap.parse_args()

    ensure_dirs()
    out_dir = os.path.join(REPORTS_DIR, "genomes")
    os.makedirs(out_dir, exist_ok=True)

    try:
        split = load_split()
    except Exception:
        split = make_split()
    reps = split[split["is_representative"]]
    if args.partition != "all":
        reps = reps[reps["partition"] == args.partition]
    reps = reps.sample(min(args.n, len(reps)), random_state=args.seed)

    engine = Engine()
    index = []
    t0 = time.time()
    for i, row in enumerate(reps.itertuples(), 1):
        gid = str(row.genome_id)
        try:
            report = engine.predict_genome(gid)
        except Exception as exc:  # noqa: BLE001
            print(f"  ! {gid}: {exc}")
            continue
        d = report.as_dict()
        d["mlst_group"] = row.mlst_group
        d["partition"] = row.partition
        with open(os.path.join(out_dir, f"{gid}.json"), "w") as fh:
            json.dump(d, fh)
        index.append({"genome_id": gid, "mlst_group": row.mlst_group,
                      "partition": row.partition, "n_proteins": d["n_proteins_predicted"],
                      "summary": d["summary"]})
        if i % 20 == 0:
            print(f"  {i}/{len(reps)} ({(time.time()-t0)/i:.2f}s/genome)")
    index.sort(key=lambda x: x["genome_id"])
    with open(os.path.join(REPORTS_DIR, "index.json"), "w") as fh:
        json.dump({"genomes": index, "count": len(index)}, fh, indent=2)
    print(f"\nWrote {len(index)} reports -> {out_dir}  in {time.time()-t0:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
