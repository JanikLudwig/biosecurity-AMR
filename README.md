# 🧬 Genome Firewall — *S. aureus* prototype

**An AI defense system against superbugs.** Turn a reconstructed *Staphylococcus
aureus* genome into an *earlier* antibiotic-response prediction — **likely to
fail / likely to work / no-call** per drug — with a calibrated confidence, an
evidence category, the genes behind the call, and a mandatory
"confirm with standard laboratory testing" notice.

> ⚠️ **Research prototype — decision support only.** Every result must be confirmed
> with standard laboratory testing before any treatment decision. The system is
> **strictly defensive**: it only predicts and explains resistance that already
> exists. It never designs, modifies, or optimizes an organism.

This build follows the **Gemini design pipeline** (BV-BRC → AMRFinderPlus features
→ per-antibiotic predictor → confidence → UI) **plus the piece that was missing:
a Drug-Target Detector (M2)** that proves the antibiotic's molecular target is
physically present in the genome — so *likely to work* is **never** asserted from
the mere absence of resistance markers.

> The `genome_firewall/` package is the earlier zero-shot v0 and is **superseded**;
> the current system lives in `gfw/`.

---

## Architecture — deterministic modules, not LLM agents

A biosecurity decision tool must be auditable, so **no autonomous LLM sits in the
decision path**. The system is a pipeline of deterministic modules:

```
Branch A — predictor                         Branch B — target evidence

genome ID + FNA                              genome ID + FNA
   │                                             │
   ▼                                             ▼
[M1 Genome Reader] teammate-owned             [M2 Target Detector]
AMRFinderPlus TSVs → binary AMR matrix        pyrodigal ORFs → pyhmmer vs curated
   │                                             targets → per-drug evidence
   ▼                                             │
[M3 Predictor] per-drug logistic regression     │
+ calibration → P(resistant)                    │
   │                                             │
   └──────────────────────┬──────────────────────┘
                          ▼
                 [M4 Decision — only join]
                 M3 P(resistant) + M2 evidence
                 → likely to fail / likely to work / no-call
                          ▼
                 [M5 Report] JSON / Markdown → dashboard
```

M3 is trained, calibrated, and scored **only** from M1 AMRFinder-derived
feature matrices and laboratory labels. M2 is an independent inference branch:
it is never a feature column or input to training, calibration, or the
logistic-regression probability. M4 uses M2 only to withhold a candidate
*likely to work* call when the target is absent or cannot be proven present;
M3's *likely to fail* call remains M3-driven and shows M2 as context.

| Module | File | Owner | Status |
|---|---|---|---|
| M1 Genome Reader (AMRFinderPlus → features) | `gfw/m1_adapter.py` (contract only) | **teammates** | consumed via contract |
| **M2 Drug-Target Detector** | `gfw/targets/` | **this work** | ✅ real |
| M3 Predictor (LR + calibration) | `gfw/predict.py` | this work | ✅ real |
| M4 Decision layer (3-way call) | `gfw/decide.py` | this work | ✅ real |
| M5 Report | `gfw/report.py`, `gfw/engine.py` | this work | ✅ real |
| Web dashboard | `web/` | this work | ✅ real |

### What's real vs. placeholder
Everything is real **except the M1 feature matrix**, which teammates own. Until
their AMRFinderPlus output lands, a **clearly-labelled synthetic** matrix
(`scripts/make_placeholder_features.py`, tagged `__synthetic__=1`) stands in so
M3–M5 and the UI can be exercised. Its metrics are **illustrative only** and every
report/plot says so.

### M1 integration contract

Teammates produce one standard AMRFinderPlus TSV per genome at
`data/amrfinder/<genome_id>.tsv`. `gfw.m1_adapter.fold_amrfinder_dir()` folds
those TSVs into `data/artifacts/features.parquet`: `genome_id` rows, AMR gene or
mutation-symbol columns, and binary presence/absence values. `scripts/train.py`
fits M3 **only** on that matrix plus lab labels. At inference, the genome's M1 row
is mapped into each saved model's fixed feature-column order; missing symbols are
zero. Neither this integration nor this repository executes AMRFinderPlus.

To replace the placeholder: have teammates supply the TSVs, run the fold step,
and re-run `scripts/train.py`. **No other code changes** — M2, split, decision
logic, and UI stay identical.

---

## The M2 Drug-Target Detector (the new piece)

For each antibiotic, `gfw/targets/specs.py` records its molecular target as
detectable **protein genes** (e.g. fluoroquinolones→`gyrA`+`grlA`,
co-trimoxazole→`folA`+`folP`, β-lactams→`pbpB`+`pbpA`, macrolides→`rplD`+`rplV`),
or as a **membrane / cell-wall** target with no single ORF (daptomycin,
vancomycin → gate *not applicable*).

`gfw/targets/detector.py` then:
1. calls ORFs from the assembly with **pyrodigal**,
2. searches the curated *S. aureus* target references (`data/references/targets/`)
   against that proteome with **pyhmmer**, and
3. returns per drug `present / absent / not_applicable` with the matched ORF,
   its **contig**, **percent identity**, and **coverage** — fully auditable.

A **likely to work** call requires `target_status == present` and cites the
detected protein; otherwise the system returns **no-call** rather than a false
"works". `~1.5 s/genome` on CPU.

---

## Honest generalization (leakage-safe split)

`gfw/split.py`: collapse near-identical genomes by **cgMLST hc10**, then hold out
whole **MLST sequence-type lineages** for calibration and test (hc10 is too fine —
nearly every genome is its own cluster). Reported metrics are on clonal groups
unseen in training. Split is deterministic (seeded).

## Data-driven panel

`gfw/panel.py` includes **every antibiotic in the lab TSV**, tiered by balanced
lab evidence: **Tier A** (train+calibrate), **Tier B** (low-power, biased to
no-call), **Tier C** (structural no-call — "insufficient lab evidence"). Uses only
the organizer-pinned **Laboratory-Method** results, never computational phenotypes.

---

## Quickstart

```bash
python -m venv --system-site-packages .venv
.venv/bin/pip install -r requirements.txt

# One-time: curated target references (needs network to UniProt)
.venv/bin/python scripts/fetch_references.py

# M1 placeholder until teammates' AMRFinderPlus output is available
.venv/bin/python scripts/make_placeholder_features.py

# Train (saves panel + split + one calibrated model per modelable drug)
.venv/bin/python scripts/train.py

# Score on the grouped hidden test → reports/metrics.{json,csv} + plots
.venv/bin/python scripts/evaluate.py

# Predict one genome end-to-end (Markdown or --json)
.venv/bin/python -m gfw.cli predict --genome 1280.10000

# Precompute reports for the dashboard, then serve it
.venv/bin/python scripts/precompute_reports.py --n 150 --partition test
.venv/bin/python -m uvicorn web.api:app --host 0.0.0.0 --port 8000
#   → open http://localhost:8000  (Predict · Performance · How it works)
```

## Success criteria (per drug, grouped hidden test)

`scripts/evaluate.py` reports balanced accuracy, **recall for resistant and
susceptible separately**, F1, AUROC, **PR-AUC**, **Brier score + reliability
curve**, and the **no-call rate** with accuracy on the called subset.

## Responsibility

Defensive by construction · calibrated confidence with an explicit no-call ·
evidence categories (known determinant (i) / statistical (ii) / none (iii)) kept
distinct · human oversight required (lab-confirmation banner on every report).

## Layout

```
gfw/            io/ (fasta+QC, labels)  m1_adapter.py  targets/ (M2)
                predict.py (M3)  decide.py (M4)  report.py engine.py (M5)  cli.py
                panel.py  split.py  evaluate.py  config.py
scripts/        fetch_references  make_placeholder_features  train  evaluate  precompute_reports
web/            api.py (FastAPI)  static/ (SPA dashboard)
data/           references/targets/  amrfinder/ (teammates drop TSVs)  artifacts/
models/  reports/
```
