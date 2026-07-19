# Genome Firewall

Genome Firewall is a defensive research prototype that predicts antimicrobial resistance
for **Staphylococcus aureus** from an assembled genome. It is decision support only: every
result must be confirmed by standard laboratory susceptibility testing.

The project is currently building the reproducible data and annotation foundation. It uses
laboratory-measured BV-BRC phenotypes, AMRFinderPlus genes/mutations, homology-aware data
splits, calibrated per-antibiotic models, deterministic target gates, and an explicit no-call.

## Quick start

```bash
uv sync
uv run genome-firewall data summarize
uv run genome-firewall data audit
uv run genome-firewall data select
uv run streamlit run sandbox/app.py
```

## Generate a decision report

The end-to-end path accepts an already assembled *S. aureus* FASTA, runs assembly QC,
AMRFinderPlus, molecular-target searches, feature-novelty checks, and the calibrated models:

```bash
uv run genome-firewall predict data/raw/genomes/1280.9342.fna
```

Reports are written under `artifacts/reports/` as JSON and CSV. A known resistance element,
statistical association, and absence of a known resistance signal are reported as distinct evidence
categories. Failed QC, an unverified target, an unseen AMR feature profile, conflicting evidence, or
a probability inside the calibrated uncertainty interval produces `no_call`.

Frozen reports from a separately deduplicated external cohort can be scored without retraining:

```bash
uv run genome-firewall model evaluate-reports \
  --reports artifacts/reports/external \
  --phenotypes data/external/phenotypes.csv
```

The current label acceptance and mixed-standard limitations are pinned in
`docs/phenotype-policy.md`.

## FastAPI: M1 + M2

The API loads the versioned model bundle once at startup. `M1` is the AMRFinderPlus
gene/mutation workflow; `M2` calls ORFs with Pyrodigal and verifies protein targets with PyHMMER,
using nucleotide search for the erythromycin 23S rRNA target.

```bash
uv run genome-firewall api
```

Open <http://127.0.0.1:8000> for the small upload UI, or call it directly:

```bash
curl http://127.0.0.1:8000/api/v1/health
curl http://127.0.0.1:8000/api/v1/models
curl -X POST http://127.0.0.1:8000/api/v1/analyses \
  -F fasta=@data/raw/genomes/1280.9342.fna
```

The JSON response keeps `workflows.M1` and `workflows.M2` separate and then combines them only in
the final target-gated `decisions` list. Generated API reports are stored under
`artifacts/api-reports/`. All model output is research-only and requires standard laboratory
confirmation.

## Rebuild the 3k model bundle

The checked-in training-v1 tables already contain AMRFinderPlus features, so these steps do not
rerun AMRFinderPlus over the training cohort:

```bash
# Normalize IDs, feature schema, laboratory labels, and the download manifest
uv run genome-firewall data prepare-training-v1

# Resumable: existing FASTAs are reused
uv run genome-firewall data download \
  --manifest data/processed/training-v1/genomes.csv \
  --output data/raw/genomes \
  --qc-output data/processed/training-v1/qc.csv

# Whole-genome sourmash grouping and leakage-safe partitions
uv run genome-firewall data cluster-split \
  --qc-manifest data/processed/training-v1/qc.csv \
  --phenotypes data/processed/training-v1/phenotypes.csv \
  --output data/processed/training-v1/splits

# Calibrated models for the supported three-drug panel
uv run genome-firewall model train \
  --features data/processed/training-v1/features.parquet \
  --phenotypes data/processed/training-v1/phenotypes.csv \
  --splits data/processed/training-v1/splits/genome-splits.csv \
  --output artifacts/models-v1 \
  --evaluation-status grouped-held-out \
  --antibiotic cefoxitin \
  --antibiotic ciprofloxacin \
  --antibiotic erythromycin

uv run genome-firewall model lineage \
  --splits data/processed/training-v1/splits/genome-splits.csv \
  --fasta-directory data/raw/genomes \
  --output artifacts/models-v1/lineage-reference.joblib
```

Until the full assembly download and regrouping complete, a locally generated bundle may be
marked `provisional-500-genome-grouped-development` in `/api/v1/models`; the API exposes that label
so provisional results cannot be mistaken for the final 3k evaluation.


The supplied source export is expected at `data/BVBRC_genome_amr.csv`. Frozen selection and
QC manifests are written under `data/manifests/` and should be committed; downloaded FASTA
files remain ignored because they are large and reproducible from those manifests.

## AMRFinderPlus without Docker

AMRFinderPlus is a compiled bioinformatics program, not a Python library. `uv` manages this
application's Python dependencies; a small Bioconda environment manages the AMRFinderPlus
binary and its native dependencies:

```bash
./scripts/setup-amrfinder.sh
uv run genome-firewall amrfinder doctor
```

The setup script supports `micromamba`, `mamba`, or `conda`. The executable is placed under
`.tools/amrfinder/`; Docker and Podman are not required. The AMRFinder database must also be
versioned and recorded for every dataset build. Database setup will be added alongside the
FASTA acquisition stage.

### Explore AMRFinder output and feature construction

```bash
# Reuse the existing smoke-genome annotation
uv run python scripts/explore_amrfinder.py

# Rerun AMRFinderPlus and show every output row
uv run python scripts/explore_amrfinder.py \
  data/raw/genomes/1280.9342.fna --force --show-all
```

The explorer prints the raw AMRFinder report, normalized gene/mutation evidence, and the
one-genome binary feature vector. It also writes small CSV artifacts under
`data/processed/amrfinder-explorer/` for manual inspection.

## Dataset commands

```bash
# Inspect the supplied BV-BRC laboratory phenotype export
uv run genome-firewall data summarize

# Preserve and inspect MIC, testing-standard, and laboratory-method provenance
uv run genome-firewall data audit

# Rebuild the deterministic 2,000-genome cohort
uv run genome-firewall data select

# Fast acquisition/QC smoke run
uv run genome-firewall data download --limit 5 --qc-output data/manifests/qc-smoke.csv

# Full acquisition/QC run (several GB; cached and resumable)
uv run genome-firewall data download
```

The downloader requests an explicit 10,000-record sequence limit from the BV-BRC API. This is
necessary because its default response limit can otherwise silently return only the first 25
contigs. Downloaded assembly length and contig count are checked against API metadata.

## Reproducible development run

```bash
# Sample broadly across the fixed 2,000-genome manifest
uv run genome-firewall data download \
  --limit 100 --sample-seed 42 --qc-output data/manifests/qc-dev-100.csv

# AMRFinder is resumable and records software/database provenance per genome
uv run genome-firewall amrfinder batch \
  --qc-manifest data/manifests/qc-dev-100.csv --workers 4

# Keep near-identical genomes in a single partition
uv run genome-firewall data cluster-split \
  --qc-manifest data/manifests/qc-dev-500.csv --output data/processed/splits-500

# Recheck phenotype support without recomputing sequence sketches
uv run genome-firewall data split-support \
  --splits data/processed/splits-500/genome-splits.csv

# Train regularized per-drug baselines and calibrate only on the calibration split
uv run genome-firewall model train \
  --splits data/processed/splits-500/genome-splits.csv

# Package train-genome sketches and calibrate the lineage novelty floor
uv run genome-firewall model lineage
```

The 100-genome development run is plumbing validation, not a benchmark. If a calibration
partition lacks the configured minimum support for either class, the model is recorded as
uncalibrated and every decision is forced to `no_call`.

## Initial scope

- Species: *Staphylococcus aureus* (NCBI taxon 1280)
- Maximum selected genomes: 2,000
- Antibiotics: erythromycin, ciprofloxacin, gentamicin, tetracycline, clindamycin
- Labels: only explicit `Resistant` and `Susceptible` laboratory results
- Excluded initially: blank, intermediate, nonsusceptible, and contradictory labels

This repository never designs, modifies, or recommends modifications to an organism.
