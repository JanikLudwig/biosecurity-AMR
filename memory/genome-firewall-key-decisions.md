---
name: genome-firewall-key-decisions
description: "Design decisions for Genome Firewall — panel tiers, MLST split, M2 stack, how to run"
metadata: 
  node_type: memory
  type: project
  originSessionId: 9318c317-d76e-47bb-b922-2fc433f264f8
  modified: 2026-07-19T00:48:33.654Z
---

Design choices baked into `gfw/` (see [[project-genome-firewall]]):

- **Panel is data-driven** from the lab TSV (`gfw/panel.py`): every antibiotic included, tiered
  by balanced lab evidence — Tier A (train+calibrate, both classes ≥100), Tier B (low-power,
  biased to no-call), Tier C (structural no-call "insufficient lab evidence").
- **Leakage-safe split** (`gfw/split.py`): dedup near-identical genomes by **cgMLST hc10**
  (nearly all singletons here), then hold out whole **MLST sequence-type** lineages for
  calibration/test (70/15/15, seeded). hc10 too fine to group; MLST is the grouping unit.
- **M2 stack:** pyrodigal (ORFs) + pyhmmer (phmmer vs curated S. aureus target proteins in
  `data/references/targets/`, fetched from UniProt via `scripts/fetch_references.py`).
  Present if identity ≥0.80, ref-coverage ≥0.60, E ≤1e-10. Membrane/cell-wall drugs
  (daptomycin, vancomycin) → target gate `not_applicable`.
- **M3:** one L2 LogisticRegression per drug + calibration (isotonic if calib minority ≥25 else
  sigmoid) on the calibration split.
- **Run order:** `scripts/fetch_references.py` → `make_placeholder_features.py` →
  `train.py` → `evaluate.py`; UI: `precompute_reports.py` then `uvicorn web.api:app`.
  Always use `.venv/bin/python` (pyrodigal/pyhmmer/sklearn live there).
