# Phenotype label policy

Version: `2026-07-19.1`

This prototype predicts laboratory antimicrobial susceptibility for *Staphylococcus aureus*.
It does not predict patient-level treatment success.

## Accepted training labels

- Evidence must equal `Laboratory Method` in the pinned BV-BRC export.
- Only explicit `Susceptible` and `Resistant` labels enter model fitting.
- `Intermediate`, missing, non-susceptible, and other labels are excluded, not coerced.
- Duplicate rows with the same genome, antibiotic, and label collapse to one observation.
- A genome-antibiotic pair containing both susceptible and resistant labels is excluded as a conflict.
- Computational phenotype predictions never enter the ground truth.

## Testing standards and measurements

The current baseline preserves submitter-provided categorical laboratory labels. It records EUCAST,
CLSI, BSAC, method, year, MIC/disk-diffusion measurement, source, and publication provenance in the
phenotype audit, but does not reinterpret raw measurements or silently convert between standards.

Consequently, mixed or missing standards remain a documented label-quality limitation. A future
MIC-derived label release must pin a breakpoint table version and keep that label set separate from
this baseline.

## Evaluation boundary

- Genetically related groups may not cross train, calibration, and test partitions.
- Calibration learns probability mappings and no-call boundaries only.
- The test partition is not used to adjust features, thresholds, target rules, or novelty floors.
- An external cohort must be deduplicated against development genomes before final validation.
