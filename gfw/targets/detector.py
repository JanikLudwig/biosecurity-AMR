"""M2 — Drug-Target Detector (deterministic, auditable).

Given a genome assembly, prove whether each antibiotic's molecular target is
physically present by:

  1. calling protein-coding genes from the FASTA with **pyrodigal** (ORF finder),
  2. searching the curated *S. aureus* target proteins against that predicted
     proteome with **pyhmmer** (phmmer), and
  3. mapping the per-gene protein hits onto each drug's target spec.

The output is the deterministic gate the decision layer needs: a drug may only be
called *likely to work* when its target is provably ``present`` — never from the
mere absence of resistance markers. Every call is auditable down to the reference
gene, the matched ORF, its contig, percent identity and coverage.
"""

from __future__ import annotations

import glob
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import pyrodigal
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence
from pyhmmer.plan7 import Background, Pipeline

from ..config import REFERENCES_DIR
from .specs import DRUG_TARGETS, TargetSpec, spec_for

# Presence thresholds. Essential orthologs are highly conserved within a species
# (~>90% identity), so these comfortably separate a real target from noise while
# flagging fragmented / out-of-scope genomes where the ORF is missing or partial.
MIN_IDENTITY = 0.80
MIN_REF_COVERAGE = 0.60
MAX_EVALUE = 1e-10

_ALPHABET = Alphabet.amino()


# --------------------------------------------------------------------------- #
# Reference loading.
# --------------------------------------------------------------------------- #
def _read_fasta(path: str) -> Tuple[str, str]:
    header, seq = "", []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].strip()
            else:
                seq.append(line.strip())
    return header, "".join(seq)


@dataclass
class TargetReference:
    gene: str
    header: str
    seq: str

    @property
    def length(self) -> int:
        return len(self.seq)


def load_references(refs_dir: str = REFERENCES_DIR) -> Dict[str, TargetReference]:
    """Load one reference protein per target gene from ``refs_dir``."""
    refs: Dict[str, TargetReference] = {}
    for path in sorted(glob.glob(os.path.join(refs_dir, "*.fasta"))):
        gene = os.path.splitext(os.path.basename(path))[0]
        if gene == "targets":  # skip the combined file
            continue
        header, seq = _read_fasta(path)
        if seq:
            refs[gene] = TargetReference(gene=gene, header=header, seq=seq)
    return refs


# --------------------------------------------------------------------------- #
# Gene calling + homology search.
# --------------------------------------------------------------------------- #
@dataclass
class ProteinHit:
    gene: str
    present: bool
    identity: float
    ref_coverage: float
    evalue: float
    bitscore: float
    orf_id: Optional[str] = None
    contig: Optional[str] = None

    def as_dict(self) -> Dict[str, object]:
        return {"gene": self.gene, "present": self.present,
                "identity": round(self.identity, 4),
                "ref_coverage": round(self.ref_coverage, 4),
                "evalue": self.evalue, "bitscore": round(self.bitscore, 1),
                "orf_id": self.orf_id, "contig": self.contig}


@dataclass
class TargetDetection:
    """All per-gene target-protein evidence for one genome."""

    genome_id: str
    n_proteins: int
    hits: Dict[str, ProteinHit] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)

    def present_genes(self) -> List[str]:
        return [g for g, h in self.hits.items() if h.present]

    def as_dict(self) -> Dict[str, object]:
        return {"genome_id": self.genome_id, "n_proteins": self.n_proteins,
                "present_genes": self.present_genes(),
                "hits": {g: h.as_dict() for g, h in self.hits.items()},
                "warnings": list(self.warnings)}


def call_proteins(contigs) -> List[Tuple[str, str]]:
    """Predict ORFs and return (orf_id, protein_seq). ``contigs`` = list of
    objects with ``.name`` and ``.seq`` (see :mod:`gfw.io.fasta`)."""
    seqs = [c.seq for c in contigs if c.seq]
    total = sum(len(s) for s in seqs)
    finder = pyrodigal.GeneFinder(meta=False)
    if total >= 20_000:
        # Single-genome mode: train on the whole assembly for best accuracy.
        finder.train(*[s.encode() if isinstance(s, str) else s for s in seqs])
    else:
        finder = pyrodigal.GeneFinder(meta=True)
    proteins: List[Tuple[str, str]] = []
    for c in contigs:
        if not c.seq:
            continue
        genes = finder.find_genes(c.seq.encode() if isinstance(c.seq, str) else c.seq)
        for i, gene in enumerate(genes, start=1):
            proteins.append((f"{c.name}_{i}", gene.translate().rstrip("*")))
    return proteins


def _identity_and_coverage(hit, ref_len: int) -> Tuple[float, float]:
    """Percent identity over the aligned region and coverage of the reference."""
    dom = hit.best_domain
    aln = dom.alignment
    idseq = aln.identity_sequence or ""
    aligned = len(idseq)
    identical = sum(1 for ch in idseq if ch.isalpha())
    identity = identical / aligned if aligned else 0.0
    coverage = (aln.hmm_to - aln.hmm_from + 1) / ref_len if ref_len else 0.0
    return identity, coverage


def search_targets(proteins: List[Tuple[str, str]],
                   references: Dict[str, TargetReference]) -> Dict[str, ProteinHit]:
    """phmmer each reference target against the predicted proteome; keep the best
    hit per gene and decide presence from identity / coverage / E-value."""
    # Build a digital database of the predicted proteome once.
    db_seqs = []
    for orf_id, aa in proteins:
        if not aa:
            continue
        ts = TextSequence(name=orf_id.encode(), sequence=aa)
        db_seqs.append(ts.digitize(_ALPHABET))
    database = DigitalSequenceBlock(_ALPHABET, db_seqs)

    background = Background(_ALPHABET)
    hits_by_gene: Dict[str, ProteinHit] = {}
    for gene, ref in references.items():
        query = TextSequence(name=gene.encode(), sequence=ref.seq).digitize(_ALPHABET)
        pipeline = Pipeline(_ALPHABET, background=background)
        top = pipeline.search_seq(query, database)
        if len(top) == 0:
            hits_by_gene[gene] = ProteinHit(gene, False, 0.0, 0.0, float("inf"), 0.0)
            continue
        best = top[0]
        identity, coverage = _identity_and_coverage(best, ref.length)
        name = best.name
        orf_id = name.decode() if isinstance(name, (bytes, bytearray)) else str(name)
        contig = orf_id.rsplit("_", 1)[0]
        present = (best.evalue <= MAX_EVALUE and identity >= MIN_IDENTITY
                   and coverage >= MIN_REF_COVERAGE)
        hits_by_gene[gene] = ProteinHit(
            gene=gene, present=present, identity=identity, ref_coverage=coverage,
            evalue=float(best.evalue), bitscore=float(best.score),
            orf_id=orf_id, contig=contig)
    return hits_by_gene


def detect(genome_id: str, contigs,
           references: Optional[Dict[str, TargetReference]] = None) -> TargetDetection:
    """Full M2 pass for one genome: ORFs -> target search -> per-gene evidence."""
    references = references if references is not None else load_references()
    proteins = call_proteins(contigs)
    warnings: List[str] = []
    if not proteins:
        warnings.append("no ORFs predicted; genome may be empty or malformed")
    hits = search_targets(proteins, references) if proteins else {}
    return TargetDetection(genome_id=genome_id, n_proteins=len(proteins),
                           hits=hits, warnings=warnings)


# --------------------------------------------------------------------------- #
# Per-drug target gate.
# --------------------------------------------------------------------------- #
@dataclass
class DrugTarget:
    drug: str
    target_status: str                 # present | absent | not_applicable
    target_kind: str
    detected: List[Dict[str, object]] = field(default_factory=list)
    missing: List[str] = field(default_factory=list)
    rationale: str = ""

    def as_dict(self) -> Dict[str, object]:
        return {"drug": self.drug, "target_status": self.target_status,
                "target_kind": self.target_kind, "detected": self.detected,
                "missing": self.missing, "rationale": self.rationale}


def drug_target_status(drug: str, detection: TargetDetection,
                       spec: Optional[TargetSpec] = None) -> DrugTarget:
    """Map a genome's protein hits onto a drug's molecular target.

    ``present`` iff at least one of the drug's target genes is detected; the
    detected proteins (with identity / contig) are the positive evidence a
    *likely to work* call must cite.
    """
    spec = spec or spec_for(drug)
    if spec is None or spec.target_kind not in ("protein",):
        kind = spec.target_kind if spec else "unknown"
        return DrugTarget(
            drug=drug, target_status="not_applicable", target_kind=kind,
            rationale=(f"{kind} target has no single detectable ORF; the "
                       "presence gate does not apply to this drug."))
    detected, missing = [], []
    for gene in spec.target_genes:
        h = detection.hits.get(gene)
        if h and h.present:
            detected.append(h.as_dict())
        else:
            missing.append(gene)
    if detected:
        names = ", ".join(d["gene"] for d in detected)
        return DrugTarget(
            drug=drug, target_status="present", target_kind="protein",
            detected=detected, missing=missing,
            rationale=(f"Molecular target detected in genome ({names}); the drug "
                       "has something to act on."))
    return DrugTarget(
        drug=drug, target_status="absent", target_kind="protein",
        detected=detected, missing=missing,
        rationale=(f"Target gene(s) {spec.target_genes} not detected — assembly "
                   "may be incomplete or out-of-scope; withhold a 'works' call."))
