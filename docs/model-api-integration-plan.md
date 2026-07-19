# Model and API Integration Plan

## Goal

Turn the existing *Staphylococcus aureus* research pipeline into a versioned service that:

1. accepts one assembled, quality-checkable FASTA file;
2. runs AMRFinderPlus and deterministic molecular-target checks;
3. builds the exact feature vector expected by the deployed models;
4. returns calibrated `likely_to_fail`, `likely_to_work`, or `no_call` results;
5. exposes the result through FastAPI and a small demonstration UI; and
6. documents the technical and evaluation choices in a judge-facing video transcript.

The service remains defensive decision support. It must not identify species, process raw
clinical samples, recommend organism modifications, or replace laboratory susceptibility
testing.

## Branch review and integration decision

Do not merge either branch wholesale. Both branches contain useful work, but they implement
different package layouts and incompatible model/feature contracts.

| Source | Reuse | Do not reuse as authoritative |
|---|---|---|
| Current `models` branch | FASTA/assembly QC, native AMRFinderPlus installation and execution, homology/lineage novelty checks, decision report, explicit no-call logic, current Streamlit sandbox | Current 104-feature local model artifacts as the final 3k models |
| `origin/feat/model-integration` | Canonical 158-feature parser, `gene::<symbol>` and `mutation::<gene>::<change>` naming, calibrated model artifact loader, training metadata, native/Docker backend abstraction, inference/reporting patterns | Its forced 0.5 R/S output, random split as a final evaluation, hard-coded three-model assumptions, per-request model loading |
| `origin/version1` | Pyrodigal + PyHMMER target detection, drug-target specification structure, strict separation between statistical prediction and deterministic target gate, useful API/UI ideas | The committed models and headline metrics: their manifest states `synthetic_features: true`; old generated datasets and reports; the branch as a full package replacement |

The deployed application should retain the current `genome_firewall` package and port only the
selected components into it.

## Verified 3k data contract

`data/training-data-v1` currently contains a consistent training dataset:

- 3,356 unique genomes and 158 binary AMRFinderPlus features;
- 58 gene features and 100 mutation features;
- 7,500 unique genome-antibiotic laboratory labels;
- cefoxitin: 1,029 resistant / 660 susceptible;
- ciprofloxacin: 1,105 resistant / 1,768 susceptible;
- erythromycin: 1,270 resistant / 1,668 susceptible;
- no conflicting genome-antibiotic label pairs;
- all listed checksums and file sizes pass; and
- AMRFinderPlus 4.2.7 with database 2026-05-15.1 is recorded.

The directory does **not** currently contain serialized trained models, training reports, a fixed
train/calibration/test assignment, or homology-group IDs. Those artifacts are required before the
models can be represented as a grouped held-out evaluation.

## Target architecture

```text
FASTA upload
  -> safe upload handling + assembly QC + supported-species scope check
  -> AMRFinderPlus
  -> canonical 158-feature adapter
  -> calibrated per-antibiotic logistic models -> P(resistant)
  -> feature/lineage novelty checks
  -> deterministic target detector
       - Pyrodigal + PyHMMER for protein targets
       - nucleotide search for RNA targets such as erythromycin's 23S rRNA target
  -> decision policy + explicit no-call
  -> versioned JSON report
  -> FastAPI response and small browser UI
```

Model probability and target detection remain separate. Target evidence may gate a candidate
`likely_to_work` result, but it must never be inserted into the trained probability after the fact.
Known AMR determinants and purely statistical associations must also remain distinct evidence
categories.

## Implementation phases

### Phase 1 — Freeze the inference and artifact contracts

- Make the 158-column feature dictionary the versioned canonical schema for the three new
  models.
- Replace the current mutation key format (`mutation::<symbol>`) at this boundary with
  `mutation::<gene>::<change>`, matching `training-data-v1`.
- Introduce an artifact manifest containing schema version, antibiotic, exact ordered feature
  list and hash, estimator/calibrator type, thresholds, training-data hash, AMRFinderPlus/database
  versions, scikit-learn/joblib versions, split method, evaluation status, and metrics path.
- Validate every model and manifest during application startup. Fail startup on a schema or
  compatibility mismatch instead of silently producing a misaligned prediction.
- Load models once during the FastAPI lifespan, not once per request.

### Phase 2 — Establish honest training and evaluation artifacts

- Obtain the precomputed homology groups/fixed split for the 3,356 genomes, or rebuild them from
  the matching FASTA assemblies using the existing sourmash-based grouping workflow. The existing
  BV-BRC `genome_id` values are sufficient to download those assemblies. Preserve them as strings
  (for example, `1280.10000`) rather than allowing CSV readers to coerce them to floating-point
  numbers.
- Rebuilding homology groups does not require rerunning AMRFinderPlus or recreating the 158-feature
  matrix. Downloaded FASTAs are used only to compute fast whole-genome sketches and can then remain
  in ignored local storage.
- Freeze genome-level train, calibration, and held-out test partitions before model selection.
- Train one regularized logistic-regression model for each supported antibiotic. Fit the base
  model only on training groups and the probability calibrator only on calibration groups.
- Select `likely_to_fail` and `likely_to_work` thresholds on calibration data. Preserve a no-call
  interval instead of deriving calls from a fixed 0.5 boundary.
- Package the three model files, one shared schema, one model manifest, split provenance, and
  evaluation reports as a single versioned model bundle.
- If the existing 3k models used a random row split, retain them only as development artifacts and
  retrain before presenting held-out performance as evidence of generalization.

### Phase 3 — Integrate the molecular-target gate

- Add `pyrodigal` and `pyhmmer` to the `uv` project and port the focused target detector from
  `version1` into `src/genome_firewall/targets/`.
- Import only small, curated target reference sequences with provenance, checksums, license/source,
  and documented identity/coverage/E-value thresholds.
- Initially support targets needed by the three deployed drugs:
  - cefoxitin: native penicillin-binding proteins;
  - ciprofloxacin: DNA gyrase and topoisomerase IV; and
  - erythromycin: 23S rRNA, using a nucleotide detector rather than presenting nearby ribosomal
    proteins as the literal drug target.
- Treat detector failure, missing references, fragmented target evidence, and unsupported target
  types as `unknown`, causing a no-call for a candidate `likely_to_work` result.
- Keep detected target proteins/RNA out of the statistical model feature vector.

### Phase 4 — Build the FastAPI service

- Add an application factory and lifespan-managed service container.
- Provide these initial endpoints:
  - `GET /api/v1/health`: process, AMRFinderPlus database, target-reference, and model readiness;
  - `GET /api/v1/models`: supported species/drugs, bundle version, feature schema, calibration,
    thresholds, evaluation status, and limitations;
  - `POST /api/v1/analyses`: multipart FASTA upload and analysis;
  - `GET /api/v1/analyses/{analysis_id}`: retrieve a completed or in-progress analysis if the
    endpoint uses background jobs.
- Use generated Pydantic request/response contracts, bounded upload size, safe generated filenames,
  per-analysis temporary directories, SHA-256 analysis IDs, cleanup, subprocess timeouts, and a
  concurrency semaphore around AMRFinderPlus/target analysis.
- Return structured failures without exposing host paths or raw subprocess commands.
- Include input QC, annotation provenance, recognized/unknown features, P(resistant), final call,
  calibrated confidence, no-call reason, evidence category, supporting determinants, target
  evidence, warnings, timings, and the mandatory laboratory-confirmation notice.
- Preserve JSON reports on disk for the demo, while keeping storage behind an interface that can
  later move to object storage or a database.

### Phase 5 — Add a deliberately small frontend

- Serve a dependency-light HTML/CSS/JavaScript page from FastAPI initially; avoid introducing a
  separate frontend build unless the product direction requires it.
- Include a FASTA upload, progress/error state, model/provenance summary, and one row per antibiotic.
- Make call, confidence, evidence category, target status, and no-call reason readable without
  expanding raw AMRFinderPlus rows.
- Keep the standard-lab confirmation banner permanently visible.
- Retain the Streamlit app as a data/model sandbox rather than making it the API client contract.

### Phase 6 — Evaluate and freeze demo evidence

For each drug, report on the untouched grouped test set:

- balanced accuracy;
- resistant and susceptible recall separately;
- F1, AUROC, and PR-AUC;
- Brier score and reliability-bin data;
- no-call rate and accuracy/balanced accuracy on called cases; and
- results by held-out genetic group where sample support is sufficient.

Also report coverage and failure behavior for assembly QC, feature novelty, lineage novelty, and the
target gate. Any development-only or random-split result must be labeled as such in the API, UI,
README, and transcript.

### Phase 7 — Produce the video transcript

Create `docs/video-transcript.md` after the final held-out reports are frozen so the narration uses
real results rather than placeholders. The transcript must fit 60 seconds (roughly 125–140 spoken
words). Recommended judge-facing structure:

1. **Problem and scope (0:00–0:08):** defensive decision support from an assembled
   *S. aureus* genome, never a replacement for laboratory testing.
2. **Pipeline (0:08–0:25):** AMRFinderPlus produces known gene/mutation features; three calibrated
   logistic models estimate resistance; target detection prevents a works call based only on the
   absence of resistance markers.
3. **Trust and evaluation (0:25–0:43):** homology-grouped train/calibration/test partitions,
   class-sensitive and calibration metrics, and an explicit no-call for uncertain or novel inputs.
4. **Demo and limitation (0:43–1:00):** upload one FASTA, show three evidence-backed results, name
   the supported scope, and end with mandatory laboratory confirmation.

## Lightweight verification

Keep tests focused on failure-prone contracts rather than broad unit coverage:

- one feature-parser fixture proving exact parity with the 158-feature dictionary;
- one startup test rejecting a model/schema hash mismatch;
- one API smoke test using a tiny fixture and mocked external executables;
- one target-gate test per target type; and
- one decision-policy test covering fail, work, and each important no-call path.

In addition, run one real local end-to-end FASTA analysis before the demo and verify its JSON report
against the CLI output.

## Completion criteria

- A fresh `uv sync` installs all Python dependencies.
- FastAPI loads a versioned, validated three-model bundle at startup.
- A valid assembled *S. aureus* FASTA can be submitted and produces the documented JSON contract.
- Feature generation exactly matches the training schema.
- A `likely_to_work` call requires calibrated low resistance probability and positive target
  evidence; weak/conflicting/out-of-distribution evidence returns no-call.
- Grouped held-out metrics, reliability, and no-call coverage are generated reproducibly.
- The minimal frontend demonstrates the complete workflow.
- The final video transcript contains the frozen metrics and makes the clinical and scope limits
  explicit.

## Implementation resolution

The missing model and grouping artifacts are reproducibly generated by the implemented pipeline:

- `data prepare-training-v1` normalizes the 3,356-genome feature and label contract;
- `data download` retrieves only missing assemblies and records QC;
- `data cluster-split` creates the homology groups and fixed partitions;
- `model train` writes the three validated joblib artifacts and bundle manifest; and
- `model lineage` packages the train-genome novelty reference.

The FastAPI-served frontend is intentionally dependency-light. A provisional bundle built from
the already available 500-genome grouped cohort enables immediate API testing while the resumable
full-cohort sequence acquisition completes. Its manifest explicitly labels it as provisional.
