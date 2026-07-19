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
