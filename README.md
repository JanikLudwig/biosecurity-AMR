# GyraseX

GyraseX is a defensive antimicrobial-resistance (AMR) research prototype for **one already assembled, already identified _Staphylococcus aureus_ genome**. It combines known AMR markers, calibrated statistical models, molecular-target checks, and conservative safety gates to return one of three research-only outcomes per supported antibiotic:

- `likely_to_work` — a low resistance-probability signal cleared every required gate;
- `likely_to_fail` — a high resistance-probability signal cleared its calibrated boundary; or
- `no_call` — the evidence is uncertain, incomplete, out of distribution, or conflicting.

> **Not a clinical device.** A result is not a treatment recommendation, a patient-level efficacy prediction, or a replacement for antimicrobial-susceptibility testing (AST). Every output requires confirmation by standard laboratory AST.

This repository implements the challenge brief's genome-to-decision concept while remaining strictly defensive: it analyses existing assembled genomes and does not design, modify, or optimise organisms.

## What is and is not in scope

| In scope | Out of scope |
| --- | --- |
| A single assembled nucleotide FASTA (`.fna`, `.fa`, or `.fasta`) from an isolated _S. aureus_ genome | FASTQ/read data, metagenomes, raw clinical specimens, mixed samples, or assembly/reconstruction |
| AMR decision support for the bundle's listed antibiotics | Species identification, diagnosis, patient-specific treatment selection, or clinical deployment |
| Known resistance evidence, calibrated resistance probability, target verification, and explicit abstention | Claiming that absence of a marker proves susceptibility |

The program is configured for _S. aureus_ (NCBI taxon `1280`) but does **not** perform taxonomic identification. Supply a species-verified, single-isolate _S. aureus_ assembly; an assembly from another species can pass the generic FASTA/assembly checks but is outside the model's validated scope.

## Repository status and reproducibility

The source, frozen feature tables, model manifests, evaluation metadata, and example reports are versioned. A fresh checkout does **not** contain the large/generated files required for live inference:

| Available in the repository | Must be created or supplied locally |
| --- | --- |
| `data/training-data-v1/`: 3,356 feature rows, labels, feature dictionary, checksums, and AMRFinderPlus provenance | `data/processed/` normalised tables and grouped splits |
| `data/demo-data.zip`: three assembled demo FASTAs | Downloaded FASTAs in `data/raw/genomes/` |
| `artifacts/models-v1/`: bundle manifest, per-drug metadata, held-out predictions, reliability table, and score summary | Serialized `*.joblib` models and `lineage-reference.joblib` |
| `artifacts/reports-v1/` and `artifacts/api-reports/`: frozen example outputs | A local AMRFinderPlus installation and database |

The missing files are intentionally ignored by Git because they are generated, large, or environment-specific. Consequently, `genome-firewall api` and `genome-firewall predict` will not work from a pristine clone until a compatible bundle and lineage reference are rebuilt (or otherwise provided). The commands in [Rebuild a local bundle](#rebuild-a-local-bundle) produce them.

There are two artifact families in the tree:

- **`artifacts/models-v1/`** is the default API bundle definition: three antibiotics, 158 features, and a declared `provisional-500-genome-grouped-development` evaluation. Its serialized model files are absent from this checkout.
- **`artifacts/models/`** contains metadata for a separate five-antibiotic, 104-feature development run, but has no `bundle-manifest.json`; it is not a runnable bundle from this checkout. The CLI's historical default points here, so pass `--model-directory artifacts/models-v1` after rebuilding the v1 bundle.

## Architecture

The implementation expands the brief's three modules into five user-facing stages. `M4` is the decision module referenced in the web client; its backend implementation is [`src/genome_firewall/decision.py`](src/genome_firewall/decision.py).

```text
assembled, species-verified S. aureus FASTA
        |
        +--> assembly QC and input validation ------------------------+
        |                                                             |
        +--> M1: AMRFinderPlus --> binary AMR features --> M3 model --+--> M4 decision policy --> M5 report
        |                                                             |
        +--> M2: molecular-target verification ----------------------+
              (independent of the statistical model)
        |
        +--> feature-profile novelty + lineage novelty gates ---------+
```

| Module | Implementation | Purpose and output |
| --- | --- | --- |
| Input/QC | [`data/qc.py`](src/genome_firewall/data/qc.py), [`inference.py`](src/genome_firewall/inference.py) | Parses an assembled DNA FASTA, calculates length/contigs/N50/ambiguity/SHA-256, and stops downstream annotation when assembly QC fails. |
| M1 — resistance features | [`annotation/amrfinder.py`](src/genome_firewall/annotation/amrfinder.py) | Runs AMRFinderPlus on the nucleotide assembly, normalises detected genes and point mutations, and aligns them to the ordered model schema. It keeps known biological evidence auditable. |
| M2 — target verification | [`targets.py`](src/genome_firewall/targets.py) | Independently checks that an antibiotic's molecular target can be found. Pyrodigal + PyHMMER are used for protein targets; BLASTN is used for rRNA targets. Target evidence never changes the model probability. |
| M3 — predictor | [`modeling/baseline.py`](src/genome_firewall/modeling/baseline.py), [`modeling/bundle.py`](src/genome_firewall/modeling/bundle.py) | Loads one validated, calibrated logistic-regression model per antibiotic and estimates `P(resistant)`. The bundle validates the exact ordered feature schema and SHA-256 hash before inference. |
| M4 — decision policy | [`decision.py`](src/genome_firewall/decision.py) | Applies calibrated probability thresholds and safety gates to choose `likely_to_work`, `likely_to_fail`, or `no_call`. This is the join between the M1/M3 and M2 branches. |
| M5 — report/API | [`inference.py`](src/genome_firewall/inference.py), [`api/app.py`](src/genome_firewall/api/app.py) | Writes versioned JSON/CSV/TSV evidence and exposes the same analysis through FastAPI. |

Cross-cutting safeguards are implemented separately:

- [`splitting/homology.py`](src/genome_firewall/splitting/homology.py) groups near-identical assemblies before data splitting, reducing train/test leakage.
- [`lineage.py`](src/genome_firewall/lineage.py) checks an input genome against training-genome sketches.
- M3 also compares its binary AMR-feature profile to training profiles using Jaccard similarity.
- [`configs/drug_registry.toml`](configs/drug_registry.toml) pins the species, per-drug resistance terms, target types, and target-count requirements.

## M4: how the decision module works

M4 is deliberately conservative. It starts with the model's calibrated resistance probability and then may turn a directional signal into `no_call`; it does not turn a weak signal into a directional call.

| Order | Condition | Result |
| --- | --- | --- |
| 1 | Assembly QC fails | `no_call` with one or more QC reasons. |
| 2 | The lineage reference is unavailable or the genome is outside it | `no_call`. |
| 3 | M1 reports an unseen feature or the feature profile is below the training-similarity floor | `no_call`. |
| 4 | `P(resistant)` is at or below the drug's calibrated lower boundary, but drug-relevant AMR evidence is present | `no_call` for conflicting evidence. |
| 5 | The same low-probability signal lacks a verified M2 target | `no_call`. |
| 6 | The low-probability signal clears all gates | `likely_to_work`. |
| 7 | `P(resistant)` is at or above the calibrated upper boundary | `likely_to_fail`; the report distinguishes strong known resistance evidence from a statistical-association-only result. |
| 8 | The probability lies between valid boundaries, or a required boundary is unavailable | `no_call`. |

For a `likely_to_work` call, relevant AMR evidence is determined from the drug registry's AMRFinder subclass terms. For a `likely_to_fail` call, a supporting marker is classed as strong when it is a known point mutation or has at least 90% reference coverage and 90% identity. A strong marker does not bypass an unavailable/uncertain model boundary; uncertainty remains a `no_call`.

The `confidence` field is derived from the final direction (`1 - P(resistant)` for `likely_to_work`, otherwise `P(resistant)`; for `no_call`, the larger side). It is a report field, not a clinical confidence interval or a guarantee of treatment outcome.

## Molecular-target checks (M2)

The v1 service bundle supports the following three drugs. The source registry also contains development entries for gentamicin, tetracycline, and clindamycin, but those are not part of the v1 API bundle.

| Antibiotic | Required target evidence |
| --- | --- |
| Cefoxitin | Both native penicillin-binding proteins: `pbpA` and `pbpB` |
| Ciprofloxacin | All four: `gyrA`, `gyrB`, `grlA`, and `grlB` |
| Erythromycin | One 23S rRNA target |

Protein target calls use frozen _S. aureus_ reference sequences in [`src/genome_firewall/resources/targets/`](src/genome_firewall/resources/targets/) whose checksums are verified against [`manifest.json`](src/genome_firewall/resources/targets/manifest.json). A protein hit must meet identity ≥80%, reference coverage ≥60%, and E-value ≤1e-10. RNA targets are queried with BLASTN using the AMRFinderPlus database reference; the implementation requires ≥70% identity and ≥80% coverage.

`M2` requires the `blastn` executable alongside the configured AMRFinderPlus executable. The setup script verifies AMRFinderPlus and its database, but does not separately test `blastn`; verify it in your local environment before relying on an end-to-end M2 run.

## Data and model construction

### Data contract

The primary data workflow uses BV-BRC laboratory phenotype records and assemblies. The repository's pinned phenotype policy is [`docs/phenotype-policy.md`](docs/phenotype-policy.md): only explicit `Resistant` and `Susceptible` laboratory labels are accepted; intermediate, missing, non-susceptible, and contradictory genome–antibiotic pairs are excluded. Submitter-provided testing-standard and method provenance are retained for audit, but raw measurements are not silently reinterpreted across standards.

`data/training-data-v1/` is a separate frozen training input:

- 3,356 genomes;
- 158 binary AMRFinderPlus features: 58 genes and 100 mutations;
- AMRFinderPlus 4.2.7 with database `2026-05-15.1` recorded in `DATASET_INFO.txt`;
- a feature dictionary and SHA-256 checksums for all supplied tables.

Run `genome-firewall data prepare-training-v1` to convert these CSV/GZip inputs into the current Parquet/CSV contracts; it does not re-run AMRFinderPlus over the 3,356 genomes.

### Leakage-aware splitting and training

For a rebuilt cohort, assemblies are sketched with sourmash (k=31, scaled=1000). The code estimates ANI from k-mer Jaccard similarity, joins genomes at estimated ANI ≥99.92% by single linkage, and assigns whole groups—not individual genomes—to fixed 70% train / 15% calibration / 15% test partitions. The supplied v1 metadata describes a 495-genome, 103-group development cohort; its intended split is 347 train / 74 calibration / 74 test genomes.

One regularised, class-balanced `liblinear` logistic-regression model is fit per antibiotic (`C=1.0`, maximum 1,000 iterations, seed 42). The base model sees only training groups. If calibration contains at least five samples of each class, a sigmoid calibrator is fit only on calibration groups. Candidate lower and upper thresholds are then chosen from calibration probabilities only when each candidate region contains at least five calls and its wrong-class fraction is no more than 10% under the prototype policy. A missing safe boundary disables that direction rather than weakening the rule.

The model artifact additionally stores:

- the exact ordered feature columns and their schema hash;
- a unique training-profile matrix plus a 5th-percentile leave-one-profile-out Jaccard floor; and
- calibration status and no-call thresholds.

The lineage artifact stores training sourmash sketches and sets its minimum accepted estimated ANI to the 5th percentile of calibration genomes' nearest-training ANI. The untouched test partition produces raw/calibrated predictions, metric summaries, and ten-bin reliability data.

## Current model scores

The following are the score metadata committed in [`artifacts/models-v1/model-summary.csv`](artifacts/models-v1/model-summary.csv). They describe the **provisional 495-genome grouped-development evaluation**, not an external validation or clinical performance claim. The underlying serialized v1 models are absent from this checkout; these numbers are preserved evaluation metadata.

| Drug | Test labels (S / R) | Balanced accuracy | Resistant recall | Susceptible recall | AUROC | PR-AUC | Brier ↓ | Model-signal coverage | No-call rate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Cefoxitin | 15 / 21 | 79.0% | 71.4% | 86.7% | 0.852 | 0.884 | 0.176 | 25.0% | 75.0% |
| Ciprofloxacin | 48 / 11 | 92.3% | 90.9% | 93.8% | 0.987 | 0.956 | 0.047 | 96.6% | 3.4% |
| Erythromycin | 44 / 21 | 88.2% | 81.0% | 95.5% | 0.960 | 0.900 | 0.087 | 61.5% | 38.5% |

Balanced accuracy averages the two class recalls; PR-AUC is useful when classes are imbalanced; Brier score measures probability quality, where lower is better. Model-signal coverage and called accuracy are computed before the full M2/QC/novelty decision gates, so they are not end-to-end clinical performance measures.

| Drug | Calibrated work boundary `P(resistant)` | Calibrated fail boundary `P(resistant)` | Called-test accuracy |
| --- | ---: | ---: | ---: |
| Cefoxitin | ≤23.7% | Unavailable | 77.8% (9 calls) |
| Ciprofloxacin | ≤46.7% | ≥56.3% | 94.7% (57 calls) |
| Erythromycin | ≤25.0% | ≥89.4% | 100.0% (40 calls) |

These boundaries were selected on calibration data, not the test set. In particular, cefoxitin has no calibrated `likely_to_fail` boundary and its high no-call rate is expected behaviour, not a missing value to be filled with 0.5.

## Input requirements and _S. aureus_ limitations

For an API upload, the file must be ≤25 MiB, non-empty, use a `.fna`, `.fa`, or `.fasta` suffix, and be a well-formed nucleotide FASTA. The parser permits IUPAC DNA symbols `ACGTNRYKMSWBDHV`; invalid symbols, a missing header, or an empty record are rejected.

An end-to-end run applies these assembly thresholds:

| Check | Accepted value |
| --- | ---: |
| Total assembly length | 2,400,000–3,100,000 bp |
| Contigs | ≤300 |
| N50 | ≥10,000 bp |
| Ambiguous-base fraction | ≤1% |

When BV-BRC metadata are available during cohort download, the pipeline also checks source length/contig-count agreement, allows only `Good` genome quality, and rejects CheckM completeness <95% or contamination >5%. For a user-supplied local FASTA, only the sequence-derived checks above are available; the pipeline does not calculate CheckM values itself.

Passing these generic checks does not establish species identity, clonality, purity, clinical relevance, or membership in the model's validated population. The lineage and feature-novelty gates are therefore required to make a directional call; failure to establish any of those conditions becomes `no_call`.

## Run locally

### Prerequisites

- Python 3.12 or later and [`uv`](https://docs.astral.sh/uv/);
- one of `micromamba`, `mamba`, or `conda` to install AMRFinderPlus 4.2.7 from [`environments/amrfinder.yml`](environments/amrfinder.yml);
- `blastn` available next to the AMRFinderPlus binary for the M2 RNA check; and
- optionally, Bun for the separate React/TanStack demonstration frontend in `bio_sentinel_ui/`.

Install the Python environment and AMRFinderPlus:

```bash
uv sync
./scripts/setup-amrfinder.sh
uv run genome-firewall amrfinder doctor
```

The setup script installs under `.tools/amrfinder/`, downloads the AMRFinder database when needed, and does not require Docker or Podman.

### Inspect the included data without inference

These commands work without downloaded assemblies or serialized models:

```bash
uv run genome-firewall data summarize
uv run genome-firewall data audit
uv run genome-firewall data prepare-training-v1
```

They inspect the supplied BV-BRC export, write an AST provenance audit under `data/processed/phenotype-audit/`, and normalise the frozen v1 tables under `data/processed/training-v1/`.

### Rebuild a local bundle

The following is the complete local-rebuild path. It downloads assemblies from BV-BRC, can take substantial time and storage, and creates ignored files. It also creates a new local evaluation; do not present it as the committed provisional result unless its inputs and split provenance match.

```bash
# 1. Materialise the frozen feature/label contract.
uv run genome-firewall data prepare-training-v1

# 2. Download the matching assemblies, resume safely, and record QC.
uv run genome-firewall data download \
  --manifest data/processed/training-v1/genomes.csv \
  --output data/raw/genomes \
  --qc-output data/processed/training-v1/qc.csv

# 3. Keep near-identical genomes in one data partition.
uv run genome-firewall data cluster-split \
  --fasta-directory data/raw/genomes \
  --qc-manifest data/processed/training-v1/qc.csv \
  --phenotypes data/processed/training-v1/phenotypes.csv \
  --output data/processed/training-v1/splits

# 4. Train the v1 service panel from the supplied 158 AMR features.
uv run genome-firewall model train \
  --features data/processed/training-v1/features.parquet \
  --phenotypes data/processed/training-v1/phenotypes.csv \
  --splits data/processed/training-v1/splits/genome-splits.csv \
  --output artifacts/models-v1 \
  --antibiotic cefoxitin \
  --antibiotic ciprofloxacin \
  --antibiotic erythromycin

# 5. Build the sequence-lineage novelty reference required at inference.
uv run genome-firewall model lineage \
  --splits data/processed/training-v1/splits/genome-splits.csv \
  --fasta-directory data/raw/genomes \
  --output artifacts/models-v1/lineage-reference.joblib
```

`model train` writes the `*.joblib` files, `bundle-manifest.json`, per-drug metadata, test predictions, reliability bins, and `model-summary.csv`. The bundle loader rejects a model with a different antibiotic, feature order, or schema hash.

### Analyse one demo FASTA

The repository includes three assembled demo FASTAs in `data/demo-data.zip`. After rebuilding the compatible bundle, extract one and run:

```bash
unzip -q data/demo-data.zip -d /tmp/genome-firewall-demo

uv run genome-firewall predict \
  /tmp/genome-firewall-demo/demo-data/1280.51926.fna \
  --model-directory artifacts/models-v1 \
  --output artifacts/reports/local-demo
```

The command writes four auditable outputs for the genome ID: the raw AMRFinderPlus TSV, parsed evidence CSV, decision CSV, and JSON report. The JSON contains QC, tool/database provenance, M1/M2 evidence, lineage status, probabilities, target status, reasons, and the mandatory laboratory-confirmation warning.

### Serve the FastAPI demo

With a rebuilt `artifacts/models-v1/` directory and a ready AMRFinderPlus installation:

```bash
uv run genome-firewall api
```

Open <http://127.0.0.1:8000> for the small bundled upload page. The API provides:

| Endpoint | Purpose |
| --- | --- |
| `GET /api/v1/health` | Readiness, AMRFinderPlus/database versions, species, and loaded drugs |
| `GET /api/v1/models` | Bundle metadata, thresholds, evaluation metadata, and reliability bins |
| `POST /api/v1/analyses` | Upload one FASTA as multipart field `fasta` |
| `GET /api/v1/analyses/{analysis_id}` | Retrieve a completed cached report |
| `GET /api/v1/analyses/{analysis_id}/raw/m1` | Download raw AMRFinderPlus output |
| `GET /api/v1/analyses/{analysis_id}/raw/m2` | Retrieve structured target-detection evidence |

For example:

```bash
curl -X POST http://127.0.0.1:8000/api/v1/analyses \
  -F fasta=@/tmp/genome-firewall-demo/demo-data/1280.51926.fna
```

The API caches reports by the first 20 hexadecimal characters of the uploaded FASTA's SHA-256 digest and stores them under `artifacts/api-reports/` by default. It accepts at most two concurrent analyses and never uses a user-supplied filename as a report path.

### Optional React demonstration client

The separate `bio_sentinel_ui/` app presents the M1–M5 decision path and consumes the FastAPI service at `http://127.0.0.1:8000` by default. Start the API first, then:

```bash
cd bio_sentinel_ui
bun install
VITE_API_BASE_URL=http://127.0.0.1:8000 bun run dev
```

This client is a visualisation layer; it does not replace the API's controls or make the models more clinically validated.

## Commands and repository layout

| Location | Contents |
| --- | --- |
| `src/genome_firewall/` | Python package: CLI, data acquisition/QC, annotation, splitting, models, targets, decision logic, inference, and API |
| `configs/` | _S. aureus_ QC/splitting/model parameters and drug-target registry |
| `data/` | Supplied BV-BRC export, initial manifests, frozen v1 data tables, and compressed demo FASTAs |
| `artifacts/` | Versioned metadata, metric reports, and frozen example analyses; generated model binaries are ignored |
| `tests/` | Contract tests for FASTA QC, target detection, model/schema validation, calibration thresholds, decisions, and external-report scoring |
| `docs/` | Phenotype policy, implementation rationale, model/API design, and presentation material |
| `bio_sentinel_ui/` | Optional Bun/React/TanStack interface branded “GyraseX” |

Useful CLI groups:

```text
genome-firewall data summarize | audit | select | prepare-training-v1 | download | cluster-split | split-support
genome-firewall amrfinder doctor | run | batch
genome-firewall model train | lineage | evaluate-reports
genome-firewall predict
genome-firewall api
```

Run the automated checks after a code or configuration change:

```bash
uv run pytest
uv run ruff check .
```

## Limitations and next validation steps

- The committed v1 evaluation is provisional grouped-development evidence from 495 QC-passing genomes, not external validation or clinical validation.
- The project has an external-report scorer (`model evaluate-reports`), but a deduplicated external AST cohort has not yet been acquired and evaluated in this repository.
- Related genomes are grouped before splitting, but the current group allocation is deterministic by cluster size/hash; the roadmap identifies direct class-support constraints in grouped assignment as future work.
- The phenotype policy preserves submitter-provided categories and records mixed/missing testing standards; it does not harmonise raw MIC values or breakpoints.
- Only the three v1-bundle antibiotics are supported by the default service. The five-drug development metadata must not be mistaken for a deployable five-drug bundle.
- M2 target presence is a gate for `likely_to_work`, not evidence that an antibiotic will work; it is not included as a feature in the learned resistance probability.
- Tool/database versions and target references matter. Rebuild or compare bundles only with their recorded schema, AMRFinderPlus/database, and split provenance.

For the full technical rationale and acceptance policy, see [`docs/model-api-integration-plan.md`](docs/model-api-integration-plan.md), [`docs/implementation-roadmap.md`](docs/implementation-roadmap.md), and [`docs/phenotype-policy.md`](docs/phenotype-policy.md).
