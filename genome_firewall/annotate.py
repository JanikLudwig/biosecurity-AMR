"""Module 01 — The Genome Reader: FASTA -> normalized AMR features.

The brief's default annotation tool is **AMRFinderPlus** (NCBI, public-domain).
It also lists **cAMRah**, a curated workflow that runs six AMR tools/databases
(AMRFinderPlus, ResFinder, RGI/CARD, Abricate/NCBI, Abricate/ARG-ANNOT, BV-BRC)
and returns a consensus. Both are heavy to install, so this module is a thin,
**pluggable annotator** with three backends that all emit the same normalized
schema (:class:`AmrHit`):

    * ``amrfinderplus`` — run the ``amrfinder`` CLI on the FASTA (default when installed)
    * ``camrah``        — run cAMRah and read its consensus AMR table (richest signal)
    * ``tsv``           — ingest a *precomputed* AMRFinderPlus/cAMRah/Abricate TSV

The ``tsv`` backend is what makes the whole prototype runnable today with zero
bioinformatics install — the brief itself notes organizers may ship precomputed
AMRFinderPlus results. ``auto`` picks the best available backend.

Design choice — *why default to AMRFinderPlus, not cAMRah:* cAMRah gives a
broader, consensus call set (better recall of resistance markers, and a higher
justified confidence for "no marker found"), but requires all six tools + their
databases. For a v0 MVP we default to the single public-domain tool the brief
recommends and treat cAMRah as an opt-in upgrade (``--annotator camrah``).
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class AmrHit:
    """One AMR gene or resistance-associated mutation, normalized across tools."""

    gene: str                       # gene symbol, e.g. blaCTX-M-15 or gyrA_S83L
    element_type: str = "AMR"       # AMR / STRESS / VIRULENCE
    element_subtype: str = "AMR"    # AMR (acquired gene) / POINT (mutation)
    drug_class: str = ""            # broad class, e.g. BETA-LACTAM
    subclass: str = ""              # specific subclass, e.g. CEPHALOSPORIN
    method: str = ""                # detection method, e.g. EXACTX, POINTX, BLASTX
    identity: Optional[float] = None
    coverage: Optional[float] = None
    contig: Optional[str] = None
    source: str = "unknown"         # which tool/db produced this hit
    raw: Dict[str, str] = field(default_factory=dict)

    @property
    def is_point_mutation(self) -> bool:
        return (self.element_subtype or "").upper() == "POINT" or "_" in self.gene

    def as_dict(self) -> Dict[str, object]:
        return {
            "gene": self.gene,
            "element_type": self.element_type,
            "element_subtype": self.element_subtype,
            "drug_class": self.drug_class,
            "subclass": self.subclass,
            "method": self.method,
            "identity": self.identity,
            "coverage": self.coverage,
            "contig": self.contig,
            "source": self.source,
        }


@dataclass
class AnnotationResult:
    """Everything Module 02 needs from the annotation step."""

    hits: List[AmrHit]
    backend: str                    # backend actually used
    # 0..1 estimate of how completely resistance was screened. Higher = a
    # "no marker found" result is more trustworthy (feeds "likely to work"
    # confidence in the predictor). cAMRah screens with 6 tools -> highest.
    screening_completeness: float = 0.6
    warnings: List[str] = field(default_factory=list)

    def as_dict(self) -> Dict[str, object]:
        return {
            "backend": self.backend,
            "screening_completeness": self.screening_completeness,
            "n_hits": len(self.hits),
            "warnings": list(self.warnings),
            "hits": [h.as_dict() for h in self.hits],
        }


# Relative trust that "we found no marker" truly means susceptible, per backend.
_COMPLETENESS = {
    "camrah": 0.85,          # six-tool consensus screens most known determinants
    "amrfinderplus": 0.65,   # single curated tool, the recommended default
    "tsv": 0.6,              # precomputed table from an unspecified single tool
    "tsv:camrah": 0.85,
    "tsv:amrfinderplus": 0.65,
}


# --------------------------------------------------------------------------- #
# Header handling — AMRFinderPlus / Abricate columns vary across versions, so we
# resolve each field we need by trying a list of accepted header names.
# --------------------------------------------------------------------------- #
_COLUMN_ALIASES = {
    "gene": ["Gene symbol", "Element symbol", "GENE", "gene", "resistance_gene", "best_hit"],
    "element_type": ["Element type", "type"],
    "element_subtype": ["Element subtype", "subtype"],
    "drug_class": ["Class", "class", "resistance", "amr_class"],
    "subclass": ["Subclass", "subclass"],
    "method": ["Method", "method"],
    "identity": ["% Identity to reference sequence", "% Identity to reference",
                 "%Identity", "%IDENTITY", "pident", "identity"],
    "coverage": ["% Coverage of reference sequence", "% Coverage of reference",
                 "%Coverage", "%COVERAGE", "coverage"],
    "contig": ["Contig id", "Contig", "SEQUENCE", "contig", "sequence"],
}


def _resolve(headers: List[str], field: str) -> Optional[str]:
    lut = {h.strip().lower(): h for h in headers}
    for alias in _COLUMN_ALIASES[field]:
        if alias.strip().lower() in lut:
            return lut[alias.strip().lower()]
    return None


def _to_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    value = value.strip()
    if value in ("", "NA", "N/A", "-", "."):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_amr_table(path: str, source: str = "tsv") -> List[AmrHit]:
    """Parse an AMRFinderPlus / cAMRah / Abricate style TSV into normalized hits.

    Only ``element_type == AMR`` rows are kept (STRESS/VIRULENCE are ignored for
    antibiotic prediction). Tab-delimited is assumed; comma is a fallback.
    """
    hits: List[AmrHit] = []
    with open(path, "rt") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        delimiter = "\t" if "\t" in sample else ("," if "," in sample else "\t")
        reader = csv.DictReader(fh, delimiter=delimiter)
        headers = reader.fieldnames or []
        col = {f: _resolve(headers, f) for f in _COLUMN_ALIASES}
        for row in reader:
            gene = (row.get(col["gene"]) or "").strip() if col["gene"] else ""
            if not gene:
                continue
            etype = (row.get(col["element_type"]) or "AMR").strip() if col["element_type"] else "AMR"
            # If the table distinguishes element types, keep only AMR rows.
            if etype and etype.upper() not in ("AMR", "AMR-SUSCEPTIBLE", ""):
                continue
            hits.append(
                AmrHit(
                    gene=gene,
                    element_type=etype or "AMR",
                    element_subtype=((row.get(col["element_subtype"]) or "AMR").strip()
                                     if col["element_subtype"] else "AMR"),
                    drug_class=((row.get(col["drug_class"]) or "").strip().upper()
                                if col["drug_class"] else ""),
                    subclass=((row.get(col["subclass"]) or "").strip().upper()
                              if col["subclass"] else ""),
                    method=((row.get(col["method"]) or "").strip()
                            if col["method"] else ""),
                    identity=_to_float(row.get(col["identity"]) if col["identity"] else None),
                    coverage=_to_float(row.get(col["coverage"]) if col["coverage"] else None),
                    contig=((row.get(col["contig"]) or "").strip()
                            if col["contig"] else None),
                    source=source,
                    raw={k: v for k, v in row.items() if k},
                )
            )
    return hits


# --------------------------------------------------------------------------- #
# Live backends (used only when the CLI is installed).
# --------------------------------------------------------------------------- #
def _have(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def run_amrfinderplus(fasta_path: str, organism: Optional[str] = None,
                      extra_args: Optional[List[str]] = None) -> List[AmrHit]:
    """Run the ``amrfinder`` CLI on an assembled FASTA and parse its TSV output."""
    if not _have("amrfinder"):
        raise RuntimeError("amrfinder CLI not found on PATH")
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as tmp:
        out_path = tmp.name
    cmd = ["amrfinder", "-n", fasta_path, "-o", out_path]
    if organism:
        cmd += ["--organism", organism]
    if extra_args:
        cmd += extra_args
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    try:
        return parse_amr_table(out_path, source="amrfinderplus")
    finally:
        try:
            os.unlink(out_path)
        except OSError:
            pass


def run_camrah(fasta_path: str, workdir: Optional[str] = None) -> List[AmrHit]:
    """Run cAMRah and read its consensus AMR table.

    cAMRah's exact CLI/paths depend on the local install; we look for its
    consensus table and parse it with :func:`parse_amr_table`. If cAMRah is not
    installed this raises, and the caller falls back to another backend.
    """
    if not _have("camrah"):
        raise RuntimeError("camrah CLI not found on PATH")
    workdir = workdir or tempfile.mkdtemp(prefix="camrah_")
    subprocess.run(["camrah", "-i", fasta_path, "-o", workdir],
                   check=True, capture_output=True, text=True)
    # cAMRah writes a merged/consensus table; try common names.
    for name in ("consensus.tsv", "camrah_consensus.tsv", "summary.tsv", "results.tsv"):
        candidate = os.path.join(workdir, name)
        if os.path.exists(candidate):
            return parse_amr_table(candidate, source="camrah")
    raise RuntimeError(f"cAMRah ran but no consensus table found in {workdir}")


def annotate(fasta_path: Optional[str] = None,
             backend: str = "auto",
             tsv_path: Optional[str] = None,
             organism: Optional[str] = None,
             tsv_source: str = "tsv") -> AnnotationResult:
    """Produce normalized AMR features from a genome, via the chosen backend.

    Parameters
    ----------
    fasta_path : path to a reconstructed assembly (for live backends).
    backend    : ``auto`` | ``amrfinderplus`` | ``camrah`` | ``tsv``.
    tsv_path   : precomputed AMR table (required for ``tsv``; optional otherwise).
    organism   : e.g. ``Escherichia`` — improves AMRFinderPlus point-mutation calls.
    tsv_source : label for the tool that produced ``tsv_path`` (``amrfinderplus``/``camrah``).
    """
    warnings: List[str] = []

    if backend == "tsv" or (backend == "auto" and tsv_path and not _have("amrfinder") and not _have("camrah")):
        if not tsv_path:
            raise ValueError("backend 'tsv' requires tsv_path")
        hits = parse_amr_table(tsv_path, source=tsv_source)
        key = f"tsv:{tsv_source}" if f"tsv:{tsv_source}" in _COMPLETENESS else "tsv"
        return AnnotationResult(hits, backend=f"tsv:{tsv_source}",
                                screening_completeness=_COMPLETENESS[key], warnings=warnings)

    if backend in ("auto", "camrah") and _have("camrah") and fasta_path:
        try:
            hits = run_camrah(fasta_path)
            return AnnotationResult(hits, backend="camrah",
                                    screening_completeness=_COMPLETENESS["camrah"])
        except Exception as exc:  # noqa: BLE001 — fall through to next backend
            warnings.append(f"cAMRah backend failed: {exc}")
            if backend == "camrah":
                raise

    if backend in ("auto", "amrfinderplus") and _have("amrfinder") and fasta_path:
        try:
            hits = run_amrfinderplus(fasta_path, organism=organism)
            return AnnotationResult(hits, backend="amrfinderplus",
                                    screening_completeness=_COMPLETENESS["amrfinderplus"],
                                    warnings=warnings)
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"AMRFinderPlus backend failed: {exc}")
            if backend == "amrfinderplus":
                raise

    # Last resort: a precomputed table if we were given one.
    if tsv_path:
        hits = parse_amr_table(tsv_path, source=tsv_source)
        key = f"tsv:{tsv_source}" if f"tsv:{tsv_source}" in _COMPLETENESS else "tsv"
        warnings.append("no annotation CLI available; used precomputed TSV")
        return AnnotationResult(hits, backend=f"tsv:{tsv_source}",
                                screening_completeness=_COMPLETENESS[key], warnings=warnings)

    raise RuntimeError(
        "No usable annotation backend. Install AMRFinderPlus or cAMRah, or pass a "
        "precomputed AMR table via tsv_path (backend='tsv')."
    )
