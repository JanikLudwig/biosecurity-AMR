"""Consume **M1** (the teammate-owned AMRFinderPlus Genome Reader).

I do not run AMRFinderPlus. I define and consume the contract:

* **Input A** — one AMRFinderPlus TSV per genome at
  ``data/amrfinder/<genome_id>.tsv`` (standard AMRFinderPlus columns).
* **Input B** — a combined **feature matrix** parquet at
  ``data/artifacts/features.parquet``: rows indexed by ``genome_id``, columns =
  AMR gene / mutation symbols, values ∈ {0,1} presence/absence. A boolean
  ``__synthetic__`` attribute in the parquet metadata marks placeholder data.

:func:`fold_amrfinder_dir` turns a directory of AMRFinderPlus TSVs (Input A) into
the feature matrix (Input B), so teammates only have to drop files. Until their
real output lands, ``scripts/make_placeholder_features.py`` writes a clearly
labelled synthetic matrix with the same schema so M3–M5 can be exercised.
"""

from __future__ import annotations

import csv
import glob
import os
from typing import Dict, List, Optional, Tuple

import pandas as pd

from .config import AMRFINDER_DIR, FEATURES_PARQUET

# AMRFinderPlus column names vary slightly by version; resolve by trying aliases.
_GENE_COLS = ["Gene symbol", "Element symbol", "GENE", "gene"]
_TYPE_COLS = ["Element type", "type"]


def _resolve(headers: List[str], candidates: List[str]) -> Optional[str]:
    lut = {h.strip().lower(): h for h in headers}
    for c in candidates:
        if c.strip().lower() in lut:
            return lut[c.strip().lower()]
    return None


def parse_amrfinder_tsv(path: str) -> List[str]:
    """Return the AMR gene/mutation symbols in one AMRFinderPlus TSV."""
    genes: List[str] = []
    with open(path, "rt") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        delim = "\t" if "\t" in sample else ","
        reader = csv.DictReader(fh, delimiter=delim)
        headers = reader.fieldnames or []
        gcol = _resolve(headers, _GENE_COLS)
        tcol = _resolve(headers, _TYPE_COLS)
        if not gcol:
            return genes
        for row in reader:
            if tcol:
                et = (row.get(tcol) or "").strip().upper()
                if et and et not in ("AMR", "AMR-SUSCEPTIBLE"):
                    continue
            g = (row.get(gcol) or "").strip()
            if g:
                genes.append(g)
    return genes


def fold_amrfinder_dir(amr_dir: str = AMRFINDER_DIR) -> pd.DataFrame:
    """Fold ``<genome_id>.tsv`` files into a binary feature matrix (Input A→B)."""
    rows: Dict[str, Dict[str, int]] = {}
    for path in sorted(glob.glob(os.path.join(amr_dir, "*.tsv"))):
        gid = os.path.splitext(os.path.basename(path))[0]
        rows[gid] = {g: 1 for g in parse_amrfinder_tsv(path)}
    if not rows:
        return pd.DataFrame()
    mat = pd.DataFrame.from_dict(rows, orient="index").fillna(0).astype("int8")
    mat.index.name = "genome_id"
    return mat.sort_index()


def save_features(mat: pd.DataFrame, path: str = FEATURES_PARQUET,
                  synthetic: bool = False) -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq
    os.makedirs(os.path.dirname(path), exist_ok=True)
    table = pa.Table.from_pandas(mat.reset_index())
    md = dict(table.schema.metadata or {})
    md[b"__synthetic__"] = b"1" if synthetic else b"0"
    table = table.replace_schema_metadata(md)
    pq.write_table(table, path)


def load_features(path: str = FEATURES_PARQUET) -> Tuple[pd.DataFrame, bool]:
    """Return (feature_matrix indexed by genome_id, is_synthetic)."""
    import pyarrow.parquet as pq
    table = pq.read_table(path)
    synthetic = (table.schema.metadata or {}).get(b"__synthetic__", b"0") == b"1"
    df = table.to_pandas()
    if "genome_id" in df.columns:
        df["genome_id"] = df["genome_id"].astype(str)
        df = df.set_index("genome_id")
    return df, bool(synthetic)
