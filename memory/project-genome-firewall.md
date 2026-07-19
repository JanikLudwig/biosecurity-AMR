---
name: project-genome-firewall
description: "What the Genome Firewall build is, scope split with teammates, and the synthetic-features caveat"
metadata: 
  node_type: memory
  type: project
  originSessionId: 9318c317-d76e-47bb-b922-2fc433f264f8
  modified: 2026-07-19T00:48:08.116Z
---

Repo `/work3/janlud/biosecurity-AMR` builds **Genome Firewall** (Hack-Nation challenge 06):
predict per-antibiotic response from a bacterial genome — **likely to fail / likely to work /
no-call** — before lab results.

- **Species pivoted to *Staphylococcus aureus*** (taxon 1280); the old `genome_firewall/` v0
  (E. coli, zero-shot rules) is **superseded/ignored**. New system lives in `gfw/`.
- **Scope split:** teammates own **M1** (AMRFinderPlus feature extraction). My scope =
  **M2 Drug-Target Detector (the deliverable)**, M3 predictor, M4 decision, M5 report, web UI.
  I consume M1 via a fixed contract (`gfw/m1_adapter.py`); I do **not** run AMRFinderPlus.
- **M1 features are currently a SYNTHETIC placeholder** (`scripts/make_placeholder_features.py`,
  parquet tagged `__synthetic__=1`) so M3–M5/UI run end-to-end; its metrics are illustrative
  only. Make real by dropping AMRFinderPlus TSVs in `data/amrfinder/` → `fold_amrfinder_dir()`
  → re-run `scripts/train.py`. Nothing else changes.
- **Data (pre-defined, keep):** `genomes/` (4,523 .fna, 4,218 labelled), `bvbrc_data/` lab AMR
  labels (Laboratory-Method only), manifest (hc10 clusters), metadata (MLST + CheckM).
- Env: **no conda / no AMRFinderPlus / no bioinformatics binaries**; venv at `.venv`
  (`--system-site-packages`) with pyrodigal, pyhmmer, scikit-learn, joblib, pyarrow, fastapi.

See [[genome-firewall-key-decisions]] and [[user-deterministic-no-llm-agents]].
