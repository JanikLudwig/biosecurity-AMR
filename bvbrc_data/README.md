# BV-BRC raw data pull — Staphylococcus aureus AMR (laboratory-measured)

Raw reference data downloaded from the BV-BRC (Bacterial and Viral Bioinformatics
Resource Center, bv-brc.org) public REST API. This mirrors the exact filter used
in the BV-BRC web UI at:
`https://www.bv-brc.org/view/Taxonomy/1279#view_tab=amr&filter=eq(evidence,"Laboratory Method")`

**Note:** taxon **1279 is the genus *Staphylococcus*** (multi-species). The
challenge brief scopes v0 to **ONE bacterial species**, so this pull was narrowed
to taxon **1280 = *Staphylococcus aureus*** (the dominant species in that genus
and the one with by far the richest lab-confirmed AMR record set).

## File

`staph_aureus_amr_lab_method.tsv` — 45,876 genome-antibiotic AMR records across
**4,859 unique genomes**, evidence == `Laboratory Method` only (per the brief:
*"Use the organizer-pinned, laboratory-measured test results - NOT general
phenotype fields, which may contain model-generated predictions"* — this
excludes all `Computational Method` / `Computational Prediction` rows).

Phenotype breakdown: Susceptible 30,005 · Resistant 11,453 · Intermediate 456 ·
Nonsusceptible 18 · blank/unspecified 3,944.

Top antibiotics by record count: ciprofloxacin, gentamicin, erythromycin,
tetracycline, vancomycin, fusidic acid, rifampin, penicillin, cefoxitin,
oxacillin, clindamycin, linezolid, daptomycin, mupirocin,
trimethoprim/sulfamethoxazole, chloramphenicol, methicillin, trimethoprim,
teicoplanin, fosfomycin, levofloxacin, moxifloxacin, tigecycline.

## How it was retrieved (reproducible)

```bash
BASE='https://www.bv-brc.org/api/genome_amr/?and(eq(taxon_id,1280),eq(evidence,%22Laboratory+Method%22))'
curl -s "${BASE}&limit(25000,0)"     -H "accept: text/tsv" >  staph_aureus_amr_lab_method.tsv
curl -s "${BASE}&limit(25000,25000)" -H "accept: text/tsv" | tail -n +2 >> staph_aureus_amr_lab_method.tsv
```
(BV-BRC's Solr-backed API caps a single response at 25,000 rows; paginated with
`limit(rows,offset)`.)

- **Downloaded:** 2026-07-18
- **Source:** BV-BRC `genome_amr` API resource, https://www.bv-brc.org/api/
- **License/attribution:** BV-BRC data are NIAID-funded public data, freely
  available for reuse; cite BV-BRC (bv-brc.org) as the source. Verify current
  terms at bv-brc.org before redistribution beyond this research prototype.
- This TSV contains **AMR phenotype labels only** — no genome sequences. Genome
  FASTA/assembly files must be pulled separately (BV-BRC FTP or `genome` API
  resource, keyed by `Genome ID`) for the Module 01 annotator to run on.

## Status

This is **raw, unfiltered reference data** — not yet wired into
`genome_firewall`'s knowledge base, which currently targets *E. coli* (5
antibiotics). See project README for the species/antibiotic-panel decision
needed before this data can drive predictions.
