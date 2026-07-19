# Three-Slide Technical Storyboard and Claude Design Prompt

## Slide content

### Slide 1 — From genome to auditable features

**Headline:** Two independent evidence branches

**Primary visual:** One assembled FASTA file enters two parallel horizontal branches that merge at
the decision gate.

**M1 — Resistance feature extraction**

`FASTA → AMRFinderPlus → 158 binary gene/mutation features → per-drug resistance model`

- Detects catalogued resistance genes and point mutations.
- Produces an auditable presence/absence feature vector.
- Separates known biological evidence from statistical association.

**M2 — Molecular-target verification**

`FASTA → target search → target present / not verified`

- Pyrodigal + PyHMMER for protein targets.
- BLASTN for ribosomal RNA targets.
- Required before a **likely-to-work** result; target presence does not prove susceptibility.

**Merge rule:** Model probability + known-marker check + target gate + QC/generalization gates.

---

### Slide 2 — One policy, three honest outcomes

**Headline:** The evidence must clear calibrated boundaries

Use three equal cards connected to a central horizontal probability scale. Use green, amber, and
red only for the final status accents.

**Likely to work — Ciprofloxacin**

- Genome: `1280.51926`
- `P(resistant) = 9.1%`
- No relevant known resistance marker.
- GyrA, GyrB, GrlA, and GrlB targets verified.
- Below work boundary: **likely to work**.

**No-call — Erythromycin**

- Genome: `1280.9342`
- `P(resistant) = 34.1%`
- No relevant known resistance marker.
- 23S rRNA target verified.
- Between the 25.0% work and 89.4% fail boundaries: **no-call**.

**Likely to fail — Ciprofloxacin**

- Genome: `1280.51872`
- `P(resistant) = 96.1%`
- Known mutations: `GyrA S84L` and `ParC S80F`.
- Above the 56.3% fail boundary: **likely to fail**.

**Footer:** No marker is not proof of susceptibility. Every result requires laboratory
confirmation.

---

### Slide 3 — Training for generalization and calibrated uncertainty

**Headline:** Related genomes never cross the split

**Primary visual:** A funnel or left-to-right sequence:

`495 quality-checked genomes → 103 homology groups → train 347 | calibrate 74 | test 74`

Then show:

`158 features → one regularized logistic model per antibiotic → sigmoid probability calibration → error-constrained work/fail thresholds`

**Threshold policy**

- Learned only from the calibration partition.
- At least 5 calibration calls at a candidate boundary.
- At most 10% wrong-class calls under the prototype policy.
- Everything between valid boundaries becomes no-call.

**Evaluation panel**

- Balanced accuracy.
- Resistant recall and susceptible recall separately.
- PR-AUC for class imbalance.
- Brier score for probability quality.
- No-call rate and accuracy among called cases.

Add one small, clearly labelled snapshot:

`Ciprofloxacin · provisional grouped test: balanced accuracy 92.3% · resistant recall 90.9% · susceptible recall 93.8% · Brier 0.047`

**Footer:** Provisional grouped-development evaluation; not clinical validation.

## Copy-paste prompt for Claude Design

```text
Create exactly three polished 16:9 presentation slides for a 60-second technical hackathon video.
The product is called GyraseX. The audience is technically literate judges, not bioinformatics
specialists.

Visual direction:
- Clean scientific editorial style on a warm white background.
- Dark navy/slate typography with restrained green, amber, and red status accents.
- Use simple vector pipeline diagrams, small evidence chips, and strong information hierarchy.
- Keep all text large enough for a video; prefer diagrams over paragraphs.
- Use a subtle DNA motif only as decoration. Do not show organism engineering, mutation design,
  laboratory synthesis, pills, doctors, or photorealistic bacteria.
- Do not use the names Genome Firewall or Bio Sentinel.
- Use the exact title casing "GyraseX".
- Show "likely to work", "no-call", and "likely to fail"; do not shorten these to clinical claims.

SLIDE 1 — TWO INDEPENDENT EVIDENCE BRANCHES
Title: "From genome to auditable features"
Show one assembled FASTA entering two parallel branches that merge at a final decision gate.
Branch M1: "AMRFinderPlus" → "known AMR genes + mutations" → "158 binary features" →
"per-drug resistance probability".
Branch M2: split into "Pyrodigal + PyHMMER: protein targets" and "BLASTN: ribosomal RNA targets",
then output "target present / not verified".
At the merge, show four compact gates: calibrated probability, known-marker conflict, target
presence, and QC/generalization.
Include this concise note: "Target presence is required for likely to work; it does not prove
susceptibility."

SLIDE 2 — THREE HONEST OUTCOMES
Title: "The evidence must clear calibrated boundaries"
Use a central probability scale plus three equal case cards:
1. Green — "Likely to work". Ciprofloxacin, genome 1280.51926, P(resistant) 9.1%, zero relevant
   known markers, GyrA/GyrB/GrlA/GrlB verified.
2. Amber — "No-call". Erythromycin, genome 1280.9342, P(resistant) 34.1%, zero relevant known
   markers, 23S rRNA verified, between work ≤25.0% and fail ≥89.4%.
3. Red — "Likely to fail". Ciprofloxacin, genome 1280.51872, P(resistant) 96.1%, known GyrA S84L
   and ParC S80F mutations, above fail ≥56.3%.
Footer: "No marker is not proof of susceptibility · Confirm every result with laboratory testing."

SLIDE 3 — LEAKAGE-AWARE TRAINING
Title: "Related genomes never cross the split"
Show 495 quality-checked genomes clustered into 103 homology groups, then allocated as intact groups
to train 347, calibration 74, and untouched test 74.
Below it show: "158 features" → "regularized logistic model per drug" → "sigmoid calibration" →
"error-constrained thresholds".
Add a compact threshold callout: "Calibration only · minimum 5 calls · ≤10% wrong class · otherwise
no-call".
Add five metric chips: balanced accuracy, resistant/susceptible recall, PR-AUC, Brier score, and
called coverage.
Include one small result panel labelled "Provisional grouped test": "Ciprofloxacin — balanced
accuracy 92.3%, resistant recall 90.9%, susceptible recall 93.8%, Brier 0.047".
Footer: "Grouped-development evidence, not clinical validation."

Return the three slides as a visually consistent sequence. Do not add claims, metrics, antibiotics,
or biological mechanisms beyond the supplied content.
```
