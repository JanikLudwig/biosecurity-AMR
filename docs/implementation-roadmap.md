# Genome Firewall implementation roadmap

This roadmap turns "likely to work" into a conservative prediction of **likely
susceptible in laboratory testing**. It does not claim patient-level efficacy and every
result requires confirmation with standard antimicrobial-susceptibility testing (AST).

## Decision contract

A genome may be called likely susceptible only when all of these conditions hold:

1. It is a quality-passing, reconstructed *Staphylococcus aureus* genome.
2. The antibiotic's expected molecular target is present and evaluable.
3. No strong known resistance gene or mutation conflicts with the call.
4. A model trained on laboratory AST predicts a sufficiently low probability of resistance.
5. That probability is calibrated on a genetically group-disjoint calibration set.
6. The genome is sufficiently similar to the validated population and evidence is not conflicting.

Any failed or uncertain gate produces `no_call`. A known high-confidence resistance mechanism
may support `likely_to_fail`; absence of such a mechanism never establishes susceptibility alone.

## Delivery plan

- [x] Filter explicit laboratory evidence for one species and five antibiotics.
- [x] Download assemblies reproducibly and apply assembly QC.
- [x] Run AMRFinderPlus and build interpretable binary features.
- [x] Cluster genomes by estimated ANI and prevent homology leakage across splits.
- [x] Train per-antibiotic logistic baselines with calibration fail-closed behavior.
- [x] Audit raw AST measurements, standards, years, and laboratory methods before label collapse.
- [ ] Define a pinned phenotype policy for EUCAST/CLSI, intermediate results, duplicate tests,
      conflicting standards, and censored MIC values such as `<=0.5`.
- [x] Acquire a 500-genome development cohort; 495 assemblies passed QC.
- [x] Tune the estimated-ANI threshold to 99.92% before model fitting. The former 99% threshold
      produced only nine single-linkage groups and no gentamicin-resistant calibration examples;
      99.92% produces 103 groups and retains both classes in held-out partitions.
- [ ] Expand toward 1,000 QC-passing genomes if held-out confidence intervals remain too wide.
- [ ] Add class-support constraints directly to grouped split assignment instead of relying only
      on a post-split support check.
- [x] Emit a per-drug/per-split phenotype-support report after every clustering run.
- [x] Learn per-drug call thresholds from calibration outcomes and return all `no_call` when the
      configured wrong-class call limits cannot be met.
- [x] Evaluate susceptible and resistant boundaries independently; a missing boundary disables only
      that call direction instead of weakening the error constraint or disabling the safe direction.
- [x] Emit plot-ready held-out reliability data and accuracy/error metrics specifically among
      non-no-call results.
- [ ] Add a versioned, cited drug-target registry and deterministic target-presence gate.
- [ ] Add an out-of-distribution/feature-novelty gate and explicit evidence-conflict handling.
- [ ] Validate once on a deduplicated external NCBI Pathogen Detection AST cohort.
- [ ] Compare the AMRFinder baseline with statistical genomic features such as sparse k-mers;
      describe these only as associations, not biological causes.

## Data-source roles

- **BV-BRC:** primary laboratory phenotype and genome source.
- **EUCAST/CLSI:** versioned breakpoint interpretation; never silently mix standards.
- **NCBI Pathogen Detection AST:** deduplicated external validation, not a lookup oracle.
- **AMRFinderPlus and optionally ResFinder:** known resistance evidence.
- **ChEMBL plus reviewed primary references:** drug-to-molecular-target registry.
- **EUCAST MIC distributions/ECOFFs:** wild-type context, not a substitute for clinical breakpoints.

## Frontier-model boundary

A frontier model may retrieve and summarize cited guidelines, mechanisms, publications, and label
provenance. It may help generate a human-readable explanation or flag inconsistencies for review.
It must not convert internet search results into a patient-level efficacy verdict, override measured
AST, invent missing evidence, or design/modify an organism. Retrieved evidence should be stored with
source, date, applicability, and exact claim, then pass deterministic validation before display.

## Immediate commands

```bash
uv run genome-firewall data audit
uv run genome-firewall data split-support \
  --splits data/processed/splits-500/genome-splits.csv
uv run genome-firewall model train \
  --splits data/processed/splits-500/genome-splits.csv
uv run streamlit run sandbox/app.py
```
