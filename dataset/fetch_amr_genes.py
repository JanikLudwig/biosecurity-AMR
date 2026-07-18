"""Fetch precomputed AMR gene-presence calls from BV-BRC's sp_gene resource.

For every genome in dataset/labels_long.csv, pulls all rows where
property == "Antibiotic Resistance". BV-BRC's sp_gene table already merges
hits from AMRFinderPlus, RGI/CARD, and PATRIC's own k-mer method (mirrors the
spirit of cAMRah's multi-tool consensus) — so this gives model-ready AMR gene
features with no FASTA download and no local bioinformatics install.

Output: bvbrc_data/staph_aureus_sp_gene_amr.csv (raw, one row per gene hit;
many-to-one with genome_id since multiple tools/genes hit per genome).
"""

from __future__ import annotations

import json
import time
import urllib.request

import pandas as pd

FIELDS = ["genome_id", "gene", "product", "antibiotics", "antibiotics_class",
          "evidence", "source", "classification"]
BATCH = 50


def fetch_batch(ids: list[str]) -> list[dict]:
    id_list = ",".join(ids)
    q = (f"and(in(genome_id,({id_list})),eq(property,%22Antibiotic%20Resistance%22))"
         f"&select({','.join(FIELDS)})&limit(10000)")
    url = "https://www.bv-brc.org/api/sp_gene/?" + q
    req = urllib.request.Request(url, headers={"accept": "application/json"})
    with urllib.request.urlopen(req, timeout=90) as resp:
        return json.loads(resp.read().decode())


def main() -> None:
    labels = pd.read_csv("dataset/labels_long.csv", dtype=str)
    genome_ids = sorted(labels["genome_id"].unique())
    print(f"genomes: {len(genome_ids)}", flush=True)

    all_rows: list[dict] = []
    t0 = time.time()
    for i in range(0, len(genome_ids), BATCH):
        chunk = genome_ids[i:i + BATCH]
        for attempt in range(3):
            try:
                recs = fetch_batch(chunk)
                break
            except Exception as exc:  # noqa: BLE001
                print(f"  retry batch {i} ({attempt+1}/3): {exc}", flush=True)
                time.sleep(2)
        else:
            recs = []
            print(f"  GAVE UP on batch {i}", flush=True)
        all_rows.extend(recs)
        print(f"  ...{i+len(chunk)}/{len(genome_ids)} genomes, "
              f"{len(all_rows)} rows so far, {time.time()-t0:.0f}s", flush=True)

    print(f"DONE: {len(all_rows)} rows, {time.time()-t0:.0f}s total", flush=True)
    df = pd.DataFrame(all_rows)
    out = "bvbrc_data/staph_aureus_sp_gene_amr.csv"
    df.to_csv(out, index=False)
    print(f"wrote {out}  shape={df.shape}", flush=True)
    print(f"genomes with >=1 AMR gene hit: {df['genome_id'].nunique()}/{len(genome_ids)}", flush=True)


if __name__ == "__main__":
    main()
