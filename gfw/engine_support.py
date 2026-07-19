"""Scope + QC gate helpers (kept separate so :mod:`gfw.engine` stays readable)."""

from __future__ import annotations

from typing import Dict, Optional, Tuple

import pandas as pd

from .io.fasta import assembly_stats, qc_assembly
from .io.labels import load_metadata


def _metadata_row(genome_id: str) -> Optional[pd.Series]:
    meta = load_metadata()
    hit = meta[meta["genome_id"] == str(genome_id)]
    return hit.iloc[0] if len(hit) else None


def scope_and_qc(genome_id: str, contigs) -> Tuple[bool, Dict[str, object], str]:
    """Return (species_in_scope, qc_dict, species_name).

    Species scope uses the dataset metadata when available; a genome whose
    metadata is missing is treated as in-scope but flagged via QC.
    """
    stats = assembly_stats(contigs)
    row = _metadata_row(genome_id)
    species = str(row["species"]) if row is not None and "species" in row else ""
    comp = float(row["checkm_completeness"]) if row is not None and pd.notna(
        row.get("checkm_completeness")) else None
    cont = float(row["checkm_contamination"]) if row is not None and pd.notna(
        row.get("checkm_contamination")) else None

    scope_ok = ("aureus" in species.lower()) if species else True
    qc = qc_assembly(stats, checkm_completeness=comp, checkm_contamination=cont)
    return scope_ok, qc.as_dict(), (species or "Staphylococcus aureus")
