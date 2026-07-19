# GyraseX integration

## Scope

This frontend is a report viewer for registered GyraseX genomes. It does not
upload or analyse arbitrary sequence files yet. That boundary is intentional: the
current API only accepts a genome ID, and M1 requires a paired AMRFinder feature
row before M3 can produce a valid resistance probability.

The UI accepts a selected FASTA filename only to derive a registered genome ID
(for example, `1280.10033.fna` becomes `1280.10033`). It then calls the
existing report endpoint. It never sends file content to a different pipeline.

## Running locally

In one terminal, from the repository root:

```bash
.venv/bin/python -m uvicorn web.api:app --host 127.0.0.1 --port 8000
```

In another terminal:

```bash
cd "Bio Sentinel"
# install dependencies with the package manager used by your frontend branch
bun install
bun run dev
```

The Vite dev proxy forwards `/api/*` to `http://127.0.0.1:8000` by default. Override that target for a different local backend with `GENOME_FIREWALL_API_URL=http://host:port bun run dev`. This keeps the browser same-origin and avoids a backend CORS change. In production, serve the
frontend behind a reverse proxy that routes `/api/*` to the FastAPI service, or
provide an equivalent same-origin proxy.

## Contract used by the frontend

| Frontend need | GyraseX API source | Notes |
| --- | --- | --- |
| Registered samples | `GET /api/genomes` | Used for autocomplete only. |
| Sample/species and calls | `GET /api/report/{genome_id}` | The dashboard is driven entirely by this report. |
| Antibiotic directory | `report.decisions[]` | Status is `call`; score is `confidence * 100`. |
| AMR evidence | `supporting_determinants`, `evidence_category` | M5 currently serializes symbols, not coordinates. |
| Target mapping | `target_evidence[]` | Shows detected gene, ORF/contig, coverage, and identity. |
| Calibration metrics | `GET /api/metrics` | Matched to a decision by canonical `drug`. |
| Reliability plot | `metric.reliability` | Uses real held-out probability/frequency points. |
| Phage fallback | `summary["likely to work"] === 0` | Generic suggestion only; there is no phage matcher API. |

## Output gaps surfaced honestly

- Raw FASTA/FASTQ upload is not available from `web/api.py`. FASTQ is not
  accepted by the current core parser.
- An unknown genome ID must not be treated as a valid prediction input: M1 looks
  up a persisted feature row by ID. Without it, M3 sees an all-zero vector.
- Species in current reports is within the supported S. aureus scope; it is not
  a general sequence-species classifier.
- The M1 adapter intentionally discards AMR hit coordinates and mutation
  annotations when it creates the binary feature matrix. The UI labels those as
  not emitted rather than inventing values.
- M2 emits detected target proteins, not a guarantee that every target in a
  multi-target drug specification is present.
- M3/M5 does not expose per-sample SHAP values or feature attributions. The
  frontend displays an explicit unavailable state rather than mock bars.
- Tier-C/no-call drugs may have no held-out metric row. The UI displays
  "Not available" for those values and renders the decision rationale.

## Future upload path

A true upload flow should be an additive backend adapter, separate from this
frontend branch:

1. persist an uploaded FASTA under a controlled server-side identifier;
2. run the team-owned M1/AMRFinder handoff to produce its feature row;
3. invoke the existing `Engine.predict_genome(genome_id, fasta_path=...)`;
4. return the unchanged M5 report schema.

That preserves the M1 -> M3 -> M4/M2 -> M5 boundaries and avoids changing the
decision or model code. It should not be added until the M1 execution and
storage contract is agreed by the pipeline owners.
