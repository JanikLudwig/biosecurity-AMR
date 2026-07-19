# M1 Integration Guide — AMRFinderPlus to the Genome Firewall predictor

This guide defines how a teammate's **M1 Genome Reader** connects to Genome
Firewall. M1 is the only source of features for the trained resistance
predictor. M2 target detection is a separate branch and is **not** part of M1
training or the logistic-regression feature matrix.

## Boundary between modules

```text
Genome ID + assembly FASTA
        |
        +--> M1: AMRFinderPlus --> binary feature vector --> M3 logistic model
        |                                                  --> P(resistant)
        |
        +--> M2: target detector --> target evidence

M4 decision = P(resistant) + target evidence
```

M1/M3 decides the statistical resistance probability. M2 does not alter this
probability. M4 may convert an otherwise `likely to work` M3 result to
`no-call` when M2 cannot prove that the molecular target is present.

## 1. Input: one AMRFinderPlus TSV per genome

Place each completed AMRFinderPlus result here:

```text
data/amrfinder/<genome_id>.tsv
```

The filename stem must exactly match the genome ID used by:

- the FASTA filename in `genomes/<genome_id>.fna`;
- the laboratory-label tables;
- the training split; and
- later prediction requests.

For example, the AMRFinder output for `genomes/1280.10000.fna` must be named:

```text
data/amrfinder/1280.10000.tsv
```

The adapter accepts standard AMRFinderPlus files with either `Gene symbol` or
`Element symbol` columns. It retains AMR / AMR-susceptible element rows and
ignores non-AMR calls such as virulence entries.

## 2. Convert TSVs into the M1 feature matrix

Run this from the repository root after all desired TSVs are available:

```bash
.venv/bin/python -c '
from gfw.config import FEATURES_PARQUET
from gfw.m1_adapter import fold_amrfinder_dir, save_features

features = fold_amrfinder_dir()
if features.empty:
    raise SystemExit("No AMRFinderPlus TSVs found in data/amrfinder/")
save_features(features, FEATURES_PARQUET, synthetic=False)
print(f"Wrote {features.shape[0]} genomes x {features.shape[1]} features")
'
```

This writes:

```text
data/artifacts/features.parquet
```

The matrix contract is:

| Part | Required form |
|---|---|
| Rows | unique `genome_id` values |
| Columns | AMRFinder gene or mutation symbols |
| Values | binary presence/absence: `1` or `0` |
| Metadata | `__synthetic__ = 0` for real M1 output |

Do not add M2 target calls, phenotypic labels, split membership, probabilities,
or report fields to this matrix. It contains AMRFinder-derived M1 features only.

## 3. Validate M1 before training

```bash
.venv/bin/python -c '
from gfw.m1_adapter import load_features
features, synthetic = load_features()
assert not synthetic, "Expected real M1 features, not the placeholder"
assert features.index.is_unique, "Genome IDs must be unique"
assert set(features.stack().unique()) <= {0, 1}, "Features must be binary"
print(features.shape)
'
```

Also check that the feature-matrix IDs overlap the labelled genome IDs. A
genome without a feature row can still be reported, but M3 will see an
all-zero vector and its result should be treated as incomplete.

## 4. Train M3 on M1 only

```bash
.venv/bin/python scripts/train.py
```

`train.py` performs the following:

1. loads the real M1 matrix from `features.parquet`;
2. loads laboratory AMR labels and the leakage-safe split;
3. fits one calibrated logistic-regression model per modelable antibiotic;
4. saves each model under `models/`.

Each model records the exact M1 feature-column order it expects. At inference,
an incoming genome's AMRFinder feature row is aligned to that order; missing
features are filled with zero. M2 must not be fitted, calibrated, or appended
to this feature vector.

## 5. Inference and decision joining

For a genome with both a FASTA and an M1 feature row:

1. M1/M3 maps its AMRFinder features into each model and returns calibrated
   `P(resistant)` per antibiotic.
2. Independently, M2 reads the FASTA and returns target-presence evidence.
3. M4 joins the two outputs:
   - high `P(resistant)` → `likely to fail`;
   - low `P(resistant)` + M2 target present → `likely to work`;
   - low `P(resistant)` + M2 target absent/not proven → `no-call`;
   - uncertain probabilities → `no-call`.

M2 is therefore a downstream safety gate and evidence source; it never changes
the trained M3 probability.

## Handoff checklist

- [ ] AMRFinder TSV filename stems exactly match genome IDs.
- [ ] TSVs have been folded into `data/artifacts/features.parquet`.
- [ ] The matrix is marked `__synthetic__ = 0`.
- [ ] Rows are unique and feature values are binary.
- [ ] `scripts/train.py` was rerun after replacing placeholder features.
- [ ] Evaluation is rerun before presenting performance metrics.
- [ ] M2 target results were not added to the M1 feature matrix.

