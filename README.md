# 🧬 Genome Firewall — v0 (zero-shot)

**An AI defense system against superbugs.** Turn a reconstructed bacterial
genome into an *earlier* antibiotic-response prediction — *likely to fail /
likely to work / no-call* — with a calibrated-style confidence, an evidence
category, and a mandatory "confirm with standard lab testing" message.

This is the **v0 MVP**: a transparent, **zero-shot rule engine** built on curated
AMR knowledge. **No machine-learning training. No deep learning.** It runs on a
CPU in milliseconds and every call is fully auditable to the gene that caused it.

> ⚠️ **Research prototype — decision support only.** Every result must be
> confirmed with standard laboratory testing before any treatment decision. The
> system is *strictly defensive*: it only predicts and explains resistance that
> already exists. It never designs, modifies, or optimizes an organism.

---

## Why "zero-shot"?

The recommended baseline in the brief is one logistic-regression model per
antibiotic. That already needs a labeled training set, a grouped train/test
split, and calibration. **v0 skips training entirely**: it maps the AMR genes and
mutations that an annotation tool detects directly to antibiotics using a curated
knowledge base (gene/mutation → drug). Consequences:

* **No train/test leakage — by construction.** The brief's main "weak submission"
  failure mode (near-identical genomes in both train and test) *cannot happen*
  because there is no training set.
* **Fully explainable.** Every "likely to fail" points at the exact determinant
  (`blaCTX-M-15`, `gyrA_S83L`, …). We only ever emit evidence category *(i) known
  determinant* or *(iii) no known signal* — never *(ii) statistical association*,
  because a rules engine has no statistical model to hide behind.
* **A real, defensible MVP** you can extend into the logistic-regression baseline
  later (the feature layer and evaluation dedup are already here).

---

## The three modules (per the challenge brief)

| # | Module | Where |
|---|--------|-------|
| 01 | **Genome Reader** — FASTA → AMR features (pluggable annotator) | [`annotate.py`](genome_firewall/annotate.py), [`fasta.py`](genome_firewall/fasta.py) |
| 02 | **Predictor** — features → per-drug call + target gate + no-call | [`predict.py`](genome_firewall/predict.py), [`knowledge.py`](genome_firewall/knowledge.py) |
| 03 | **Decision Report** — calibrated-style confidence, evidence, safety banner | [`report.py`](genome_firewall/report.py), [`app.py`](genome_firewall/app.py) |

### Module 01 — Genome Reader (pluggable annotator)

A thin annotator with three interchangeable backends, all emitting one normalized
schema (`AmrHit`):

* **`amrfinderplus`** *(default)* — runs the NCBI `amrfinder` CLI (public-domain).
* **`camrah`** — runs **cAMRah**, the curated six-tool consensus workflow
  (AMRFinderPlus + ResFinder + RGI/CARD + Abricate/NCBI + Abricate/ARG-ANNOT +
  BV-BRC). Richest signal; used when installed.
* **`tsv`** — ingests a *precomputed* AMRFinderPlus / cAMRah / Abricate table, so
  the whole pipeline runs **today with zero bioinformatics install** (the brief
  notes organizers may ship precomputed AMRFinderPlus results).

**On cAMRah:** we wire it in because its broader consensus improves recall of
resistance markers *and* justifies a higher confidence when *no* marker is found
(see `screening_completeness`). But it needs all six tools + databases installed,
so v0 **defaults to AMRFinderPlus** — the single public-domain tool the brief
recommends — and treats cAMRah as an opt-in upgrade (`--annotator camrah`).

### Module 02 — Predictor (zero-shot rules + target gate)

For each drug in the panel:

1. **Sample gates** → out-of-scope species or failed assembly QC ⇒ `no-call`.
2. **Determinant match** → any known resistance gene/mutation for the drug?
   Weights combine by **noisy-OR**, scaled by alignment identity/coverage.
   * combined *p(fail)* ≥ 0.60 ⇒ **likely to fail** (evidence *(i)*)
   * 0.40–0.60, or a lone low-level marker ⇒ **no-call** (weak evidence)
3. **Target-presence gate** (deterministic) → we never say *likely to work* from
   the *absence* of markers alone; the drug's molecular target must be present.
4. **No marker + target present** ⇒ **likely to work** (evidence *(iii)*), with
   confidence **bounded by screening breadth** — absence of evidence is weaker
   than evidence of absence, and a single-tool screen is bounded lower than a
   cAMRah consensus.

### Module 03 — Decision Report

A Streamlit app (and a CLI) that shows, per drug: the call, a confidence score +
band, the evidence category, the supporting genes, the target-gate status, and an
always-on **lab-confirmation** banner.

---

## Quickstart

```bash
# 1) Core engine has NO dependencies — try it immediately on a precomputed table:
python -m genome_firewall.cli predict \
    --tsv genome_firewall/examples/ecoli_resistant_amrfinder.tsv \
    --tsv-source amrfinderplus --species "Escherichia coli"

# JSON or Markdown output:
python -m genome_firewall.cli predict --tsv <table.tsv> --format json
python -m genome_firewall.cli predict --tsv <table.tsv> --format md

# 2) From an assembly, auto-selecting AMRFinderPlus/cAMRah if installed:
python -m genome_firewall.cli predict --fasta assembly.fasta --annotator auto \
    --species "Escherichia coli" --organism Escherichia

# 3) The demo app (needs: pip install -r requirements.txt):
streamlit run genome_firewall/app.py

# 4) De-duplicate a genome collection for leakage-free EVALUATION:
python -m genome_firewall.cli dedup a.fasta b.fasta c.fasta --threshold 0.9
```

Run the tests (stdlib only, no pytest needed):

```bash
python -m unittest discover -s tests -v
```

---

## Supported scope (stated honestly)

* **Species:** *Escherichia coli* only. Any other species ⇒ `no-call`.
* **Antibiotics (5):** Ampicillin, Ceftriaxone, Ciprofloxacin, Gentamicin,
  Trimethoprim-sulfamethoxazole — spanning distinct mechanisms.
* **Out of scope:** sample collection, basecalling, species ID, genome assembly,
  or de-mixing multiple organisms. Input starts at a quality-checked assembly.

The panel and gene→drug rules live in editable data files
([`data/antibiotics.json`](genome_firewall/data/antibiotics.json),
[`data/gene_drug_rules.json`](genome_firewall/data/gene_drug_rules.json)) — extend
coverage without touching engine code.

---

## How v0 addresses the Responsibility Requirement

* **Defensive by construction** — read-only prediction of existing resistance;
  no organism design/modification anywhere in the codebase.
* **Honest generalization** — narrow, explicit scope; out-of-scope ⇒ no-call. No
  training means no over-fitting to near-identical genomes.
* **Calibrated confidence + no-call** — explicit no-call for weak/conflicting
  evidence, shallow screening, absent/unknown target, or out-of-scope input.
  Confidence is clearly labeled **rule-based and uncalibrated** in v0 (see
  Limitations) rather than dressed up as a trained probability.
* **Honest explanations** — only categories *(i)* known determinant and *(iii)*
  no known signal; never presents a statistical score as biological cause.
* **Human oversight** — a mandatory lab-confirmation banner on every report; the
  tool is framed as decision support, never a treatment decision.

## Limitations (v0)

* Confidence values are **literature-informed heuristics, not empirically
  calibrated** — v0 does not fit to labels, so there is no Brier score /
  reliability curve yet. The `PredictConfig` thresholds and rule weights are the
  knobs. Calibration against a labeled BV-BRC split is the first v1 step.
* The target-presence gate assumes the panel's **essential** targets (gyrase,
  PBPs, ribosome, DHFR/DHPS) are present in any viable in-scope genome; the
  interface accepts real target-detection evidence when a future backend supplies
  it.
* Co-trimoxazole is modeled as a documented simplification (trimethoprim `dfr` as
  the primary driver; `sul` as supporting evidence).

## Roadmap to v1

1. Fit one regularized logistic regression per antibiotic on the AMR feature
   matrix (the brief's recommended baseline) using the grouped BV-BRC split.
2. Empirically **calibrate** confidence (isotonic/Platt) and report Brier score +
   reliability plots and per-group generalization.
3. Use `genome_firewall.dedup` (or Mash/XTree) to enforce a leakage-free split.

## Project layout

```
genome_firewall/
  fasta.py        # FASTA parsing + assembly QC
  annotate.py     # Module 01: pluggable annotator (amrfinderplus|camrah|tsv)
  knowledge.py    # loads the knowledge base; matches determinants -> drugs
  predict.py      # Module 02: zero-shot rule engine + target gate + no-call
  report.py       # Module 03: text/markdown report rendering
  app.py          # Module 03: Streamlit demo UI
  dedup.py        # MinHash k-mer dedup for leakage-free evaluation
  cli.py          # command-line entry point
  data/*.json     # editable drug panel + gene->drug rules
  examples/*.tsv  # precomputed AMRFinderPlus tables (resistant/susceptible/weak)
tests/            # stdlib unittest suite
```
