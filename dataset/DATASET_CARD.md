# Genome Firewall — S. aureus AMR training dataset (v0 build)

Built from BV-BRC (bv-brc.org), following the "For organizers" checklist in the
challenge brief (Appendix): fixed species/antibiotic panel, laboratory-only
labels, one final label per genome-antibiotic pair, genome quality flags,
source/accession records, genetic-group-aware fixed splits, checksums.

## 1. Species & scope

- **Taxon:** *Staphylococcus aureus*, NCBI taxon **1280** (note: the originating
  BV-BRC URL used taxon **1279**, the *genus* Staphylococcus — narrowed to the
  single dominant species per the brief's "ONE bacterial species" requirement).
- **Antibiotic panel (5 drugs, 5 distinct resistance mechanisms):**

| Antibiotic | Mechanism class | Genomes labeled | R | S | Intermediate | Conflicting | R:S balance |
|---|---|---|---|---|---|---|---|
| Cefoxitin | Beta-lactam (mecA/mecC surrogate, MRSA) | 2,096 | 1,308 | 770 | 0 | 18 | 0.59 |
| Ciprofloxacin | Fluoroquinolone (gyrA/grlA) | 3,320 | 1,381 | 1,864 | 73 | 2 | 0.74 |
| Erythromycin | Macrolide (ermA/ermC/msrA) | 3,361 | 1,452 | 1,899 | 8 | 2 | 0.77 |
| Gentamicin | Aminoglycoside (aac(6')-aph(2'')) | 3,385 | 375 | 2,934 | 76 | 0 | 0.13 |
| Trimethoprim-sulfamethoxazole | Folate-pathway inhibitor (dfrA/dfrG, sul) | 1,536 | 290 | 1,242 | 3 | 1 | 0.23 |

Panel selection: from 14 antibiotic candidates with ≥1,000 labeled genomes,
picked 5 that (a) span mechanistically distinct drug classes — so the model
can't cheat by learning one mega-feature — and (b) have workable class
balance. Vancomycin/linezolid/rifampin were **excluded**: clinically vital but
R:S balance ≤0.02 (resistance is rare-to-near-absent in this cohort), which
would make "no-call" the only honest output for most of the confidence range —
a reasonable v1 stretch goal, not a v0 training target. Full comparison table
reproducible via `dataset/build_labels.py`.

## 2. Labels — provenance and the resolution rule

- **Source:** BV-BRC `genome_amr` REST API resource, `taxon_id=1280`,
  `evidence == "Laboratory Method"` **only** — per the brief: *"Use the
  organizer-pinned, laboratory-measured test results - NOT general phenotype
  fields, which may contain model-generated predictions."* All
  `Computational Method` / `Computational Prediction` rows are excluded.
- **Raw pull:** `bvbrc_data/staph_aureus_amr_lab_method.tsv` — 45,876 rows,
  4,859 unique genomes, downloaded 2026-07-18. See
  `bvbrc_data/README.md` for the exact reproducible `curl` query and license
  note.
- **One final label per (genome, antibiotic) pair** (`dataset/build_labels.py`):
  1. Collapse duplicate lab records for the same pair.
  2. All non-blank phenotypes agree (Resistant/Susceptible) → that label.
  3. Both Resistant *and* Susceptible occur for the same pair → `conflicting`
     (excluded from training, kept in the table for transparency — rare:
     114/39,147 R/S-bearing pairs in the raw data, ≤0.9% within any one
     panel drug — see table above).
  4. Only Intermediate/Nonsusceptible present → `intermediate` (excluded from
     the binary train set; retained as a category for future no-call/
     confidence-calibration studies, per the brief's emphasis on honest
     uncertainty).
  5. No usable phenotype at all → dropped.
- **Output:** `dataset/labels_long.csv` — 13,698 genome-antibiotic rows,
  3,893 unique genomes, columns include `pubmed_ids`,
  `laboratory_typing_method`, `testing_standard` for full traceability back to
  the source lab study.

## 3. Genome quality flags & genetic grouping — **no FASTA download needed**

BV-BRC's `genome` API resource carries precomputed per-genome QC and
**core-genome MLST hierarchical clustering** (`cgmlst_hc0/2/5/10/20/50/100`),
retrieved for all 3,893 genomes with metadata-only API calls (`dataset/
genome_metadata.csv`, `dataset/build_splits.py`) — no sequence data needed for
the split.

- **QC rule** (`qc_pass`, mirrors `genome_firewall.fasta.compute_qc`
  thresholds): `genome_quality == "Good"` AND `checkm_completeness >= 90` AND
  `checkm_contamination <= 5` AND `contigs <= 500`. **3,866/3,893 genomes pass
  (99.3%)**.
- **Genetic grouping:** `cgmlst_hc10` — genomes within 10 core-genome allele
  differences are one group. This is the conventional resolution for calling
  outbreak/transmission clusters in genomic epidemiology (tighter thresholds
  like hc0/hc2/hc5 are near-identical-only; looser ones like hc50/hc100 start
  merging unrelated lineages — at hc100 the single largest cluster already
  holds 633/3,893 genomes, which would swallow the "hidden groups" property).
  At hc10: **2,685 groups**, 2,500 of them singletons. 683 genomes (17.5%)
  lack a computed cgMLST call — each gets its own singleton pseudo-group
  (`nogroup:<genome_id>`), a conservative fallback that can only over-split,
  never leak.
- **Split:** groups shuffled (seed `20260718`) and assigned whole to one of
  `train` (70%) / `calibration` (15%) / `hidden_test` (15%) by genome count —
  **no genetic group ever appears in more than one split** (asserted in
  `build_splits.py`). Result: 2,725 / 584 / 584 genomes; 2,306 / 524 / 538
  distinct genetic groups respectively — so `hidden_test` is guaranteed to
  contain groups never seen in `train`, as the brief requires.
- **Output:** `dataset/genome_splits.csv`.

## 4. AMR gene features — also pulled via API, no FASTA/local tools needed

BV-BRC's `sp_gene` resource (`property == "Antibiotic Resistance"`) stores
precomputed AMR gene-presence calls **already run by BV-BRC using multiple
tools** — evidence values include `AMRFinderPlus: BLASTP/EXACTP/POINTP/...`,
`RGI: protein homolog model` (CARD), and PATRIC's own k-mer search — i.e. a
multi-tool consensus in the same spirit as cAMRah, retrieved with zero local
installs (`dataset/fetch_amr_genes.py`).

- Allele-numbered gene variants (RGI/CARD naming, e.g. `tetM_1`, `norB_3`) are
  collapsed to one gene-family feature (`tetM`, `norB`) — documented
  simplification; a v1 could keep per-allele resolution.
- Gene families seen in fewer than 5 genomes are dropped (documented
  threshold, `MIN_SUPPORT` in `dataset/build_features.py`) — negligible
  statistical signal at n≈3,900, mainly noise/overfit risk.
- **Output:** `dataset/features_gene_presence.csv` (one row per genome, one
  binary column per gene family) and `dataset/model_table.csv` (labels +
  QC + split + gene features joined — ready for the brief's recommended
  baseline: one regularized logistic regression per antibiotic).

- **Raw pull:** 345,413 gene-hit rows across all 3,893 genomes (every genome
  has ≥1 AMR-relevant hit). 197,980 rows resolve to a usable gene symbol
  (the rest are CARD/RGI regulatory or homolog rows with no clean gene/product
  match, dropped).
- **Naming-convention fix:** AMRFinderPlus, RGI/CARD and PATRIC's k-mer method
  format the *same* gene differently — e.g. raw data contained `ermA`,
  `Erm(A)`, `ERM(A)` as three separate strings. `_canonicalize()` in
  `build_features.py` lowercases, strips `" family"`, and removes
  parentheses/apostrophes before collapsing — this took the erm-family alone
  from 6 fragmented columns down to a clean set (`erm`, `erm33`, `erma`,
  `ermb`, `ermc`, `ermt`).
- **Result:** 167 distinct gene families → **133 kept** (support ≥5 genomes,
  34 rare ones dropped) → `dataset/features_gene_presence.csv`,
  **3,891 genomes × 133 binary gene-presence features**. Only 6/13,698
  genome-antibiotic label rows have zero gene features at all.
- **Biology sanity check** (resistant vs. susceptible mean gene-presence
  rate, spot-checked against known S. aureus resistance mechanisms):

  | Drug | Gene | Resistant | Susceptible |
  |---|---|---|---|
  | Cefoxitin | `mecA` (PBP2a, MRSA) | **95.3%** | 13.4% |
  | Erythromycin | `ermC` (23S methylase) | **47.5%** | 2.0% |
  | Erythromycin | `ermA` (23S methylase) | **26.9%** | 0.6% |
  | Gentamicin | `aac6-ie-aph2-ia` (bifunctional AME) | **45.3%** | 2.7% |
  | Co-trimoxazole | `dfrG` (DHFR) | **42.4%** | 5.0% |

  Every marker shows a large, mechanism-consistent gap — the labels and the
  gene-presence features agree with textbook *S. aureus* resistance biology,
  which is the main sanity check available before actual model fitting.

## 5. Known limitations

- Labels come from heterogeneous studies (different labs, testing standards,
  years) pooled by BV-BRC — `laboratory_typing_method` / `testing_standard`
  columns in `labels_long.csv` let you filter to a single standard (e.g. only
  CLSI) if you want a more homogeneous training set.
- `sp_gene` AMR calls are BV-BRC's own precomputed run, not one we executed
  ourselves — for full auditability/reproducibility, running AMRFinderPlus
  directly on the genome assemblies (`genome_firewall.annotate`) remains the
  gold-standard path; this dataset trades that off for zero-install speed.
- Vancomycin/linezolid/rifampin resistance is rare in this cohort — excluded
  from v0's training panel but real, high-stakes cases; worth a dedicated
  no-call-heavy treatment in v1.
- No genome sequences (FASTA) are included in this dataset — only phenotype
  labels, QC metadata, and gene-presence features, all sourced via the BV-BRC
  API. Assemblies (if needed for a from-scratch AMRFinderPlus run, or for a
  deep-learning stretch) are retrievable per genome_id from the `genome`
  resource / BV-BRC FTP.

## 6. License & attribution

BV-BRC data are NIAID-funded public data, freely available for reuse; cite
BV-BRC (bv-brc.org) as the source. See `bvbrc_data/README.md`. This dataset
build is a research artifact — not validated for clinical or regulatory use.

## 7. Reproducing this dataset

```bash
# raw label pull — see bvbrc_data/README.md for the exact curl commands
python dataset/build_labels.py      # -> labels_long.csv, panel_summary.csv
# genome QC + cgMLST metadata pull is inline in this card's step 3
python dataset/build_splits.py      # -> genome_splits.csv
python dataset/fetch_amr_genes.py   # -> bvbrc_data/staph_aureus_sp_gene_amr.csv
python dataset/build_features.py    # -> features_gene_presence.csv, model_table.csv
```

Checksums for the raw/derived label files: `dataset/CHECKSUMS.sha256`.
