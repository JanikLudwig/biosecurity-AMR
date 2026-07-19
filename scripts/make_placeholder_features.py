#!/usr/bin/env python3
"""Write a **clearly-labelled synthetic** AMR feature matrix (M1 placeholder).

Teammates own the real AMRFinderPlus feature extractor. Until it lands, this
generates a feature matrix with the *same schema* so M3 (predictor), M4
(decision) and evaluation can be built and demonstrated end-to-end.

The synthesis is biologically structured — a determinant gene is present mostly
in genomes actually resistant to a drug it explains, and mostly absent otherwise,
plus noise accessory genes — so a per-drug logistic regression recovers real-ish
signal. It is **not** real annotation: the parquet is tagged ``__synthetic__``
and every downstream report/metric surfaces that. Swap in real features by
dropping AMRFinderPlus TSVs into ``data/amrfinder/`` and running
``fold_amrfinder_dir`` (see ``gfw.m1_adapter``).
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gfw.config import FEATURES_PARQUET, ensure_dirs
from gfw.io.labels import load_lab_labels
from gfw.m1_adapter import save_features

# AMRFinderPlus-style determinant symbols per drug (primary first).
DETERMINANTS = {
    "methicillin": ["mecA"], "oxacillin": ["mecA"], "cefoxitin": ["mecA"],
    "ceftaroline": ["mecA", "pbp4_promoter"],
    "penicillin": ["blaZ"],
    "erythromycin": ["erm(C)", "erm(A)", "msr(A)"],
    "clindamycin": ["erm(A)", "erm(C)", "lnu(A)"],
    "ciprofloxacin": ["gyrA_S84L", "grlA_S80F"],
    "moxifloxacin": ["gyrA_S84L", "grlA_S80F"], "levofloxacin": ["gyrA_S84L", "grlA_S80F"],
    "gentamicin": ["aac(6')-aph(2'')", "aph(3')-III"],
    "tobramycin": ["aac(6')-aph(2'')"], "kanamycin": ["aph(3')-III"],
    "streptomycin": ["ant(6)-Ia", "str"],
    "tetracycline": ["tet(K)", "tet(M)"], "doxycycline": ["tet(M)"],
    "minocycline": ["tet(M)"],
    "trimethoprim_sulfamethoxazole": ["dfrG", "dfrA", "sul1"],
    "trimethoprim": ["dfrG", "dfrA"],
    "chloramphenicol": ["catA", "fexA"],
    "fusidic_acid": ["fusB", "fusC"], "mupirocin": ["mupA"],
    "rifampin": ["rpoB_H481Y"], "linezolid": ["cfr", "optrA"],
    "daptomycin": ["mprF_S295L"], "vancomycin": ["walK_A468T"],
    "tigecycline": ["tet(M)_mut"],
}

P_POS = 0.82   # determinant present given resistant to a drug it explains
P_NEG = 0.04   # determinant present given susceptible (background carriage)
N_BACKGROUND = 30
SEED = 1280


def main() -> int:
    ensure_dirs()
    rng = np.random.default_rng(SEED)
    lab = load_lab_labels()
    genomes = sorted(lab["genome_id"].unique())

    # gene -> set of drugs it explains
    gene_drugs: dict = {}
    for drug, genes in DETERMINANTS.items():
        for g in genes:
            gene_drugs.setdefault(g, set()).add(drug)
    det_genes = sorted(gene_drugs)

    # Per (genome, drug) resistance lookup.
    res = {(r.genome_id, r.drug): r.resistant for r in lab.itertuples()}

    all_genes = det_genes + [f"acc_{i:02d}" for i in range(N_BACKGROUND)]
    bg_prev = rng.uniform(0.08, 0.45, size=N_BACKGROUND)  # per-accessory prevalence

    data = np.zeros((len(genomes), len(all_genes)), dtype=np.int8)
    for gi, gid in enumerate(genomes):
        for ci, gene in enumerate(det_genes):
            drugs = gene_drugs[gene]
            labels = [res[(gid, d)] for d in drugs if (gid, d) in res]
            if not labels:
                p = P_NEG
            elif any(l == 1 for l in labels):
                p = P_POS
            else:
                p = P_NEG
            data[gi, ci] = int(rng.random() < p)
        for j in range(N_BACKGROUND):
            data[gi, len(det_genes) + j] = int(rng.random() < bg_prev[j])

    mat = pd.DataFrame(data, index=pd.Index(genomes, name="genome_id"), columns=all_genes)
    save_features(mat, FEATURES_PARQUET, synthetic=True)
    print(f"Wrote SYNTHETIC feature matrix: {mat.shape[0]} genomes x {mat.shape[1]} genes")
    print(f"  determinant features: {len(det_genes)} | background noise: {N_BACKGROUND}")
    print(f"  -> {FEATURES_PARQUET}  (tagged __synthetic__=1)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
