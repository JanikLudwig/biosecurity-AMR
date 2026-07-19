"""Curated *Staphylococcus aureus* drug → molecular-target knowledge.

Two knowledge maps, kept as data so a domain expert can extend them:

* ``TARGET_PROTEINS`` — the reference target proteins we detect in a genome
  (the fetch script pulls these from UniProt into ``data/references/targets/``).
* ``DRUG_TARGETS`` — for each antibiotic: its class, its molecular target
  *kind* (a detectable **protein**, or a **membrane / cell-wall** target that no
  single ORF represents), the specific **target genes** M2 should look for, and
  the **known-determinant patterns** M4 uses only to *label evidence category (i)*
  and explain a call (never to make the prediction — that is the trained model).

Target-gene choices reflect what is (a) essential and (b) representable as a
protein ORF, so pyrodigal+pyhmmer can prove its presence:
  * fluoroquinolones → DNA gyrase (gyrA) + topoisomerase IV (grlA)
  * β-lactams        → penicillin-binding proteins (PBP2=pbpB, PBP1=pbpA)
  * folate pathway   → DHFR (folA) + DHPS (folP)   [cleanest two-protein gate]
  * 23S-ribosome     → ribosomal proteins L4 (rplD) / L22 (rplV) / L3 (rplC)
  * 30S-ribosome     → ribosomal proteins S12 (rpsL) / S10 (rpsJ)
  * EF-G / IleRS / RNAP / MurA → fusA / ileS / rpoB / murA
  * membrane / cell wall (daptomycin, glycopeptides) → no single-protein gate.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple


# --------------------------------------------------------------------------- #
# Reference target proteins to detect (gene symbol -> UniProt fetch spec).
# organism_id 1280 = Staphylococcus aureus.
# --------------------------------------------------------------------------- #
TARGET_PROTEINS: Dict[str, Dict[str, str]] = {
    "gyrA": {"query": "gene:gyrA AND organism_id:1280", "desc": "DNA gyrase subunit A"},
    "gyrB": {"query": "gene:gyrB AND organism_id:1280", "desc": "DNA gyrase subunit B"},
    "grlA": {"query": "gene:grlA AND organism_id:1280", "desc": "Topoisomerase IV subunit A (ParC)"},
    "grlB": {"query": "gene:grlB AND organism_id:1280", "desc": "Topoisomerase IV subunit B (ParE)"},
    "pbpB": {"query": "gene:pbpB AND organism_id:1280", "desc": "Penicillin-binding protein 2"},
    "pbpA": {"query": "gene:pbpA AND organism_id:1280", "desc": "Penicillin-binding protein 1"},
    "pbpC": {"query": "gene:pbpC AND organism_id:1280", "desc": "Penicillin-binding protein 3"},
    "pbpD": {"query": "gene:pbpD AND organism_id:1280", "desc": "Penicillin-binding protein 4"},
    "folA": {"query": "gene:folA AND organism_id:1280", "desc": "Dihydrofolate reductase (DHFR)"},
    "folP": {"query": "gene:folP AND organism_id:1280", "desc": "Dihydropteroate synthase (DHPS)"},
    "fusA": {"query": "gene:fusA AND organism_id:1280", "desc": "Elongation factor G (EF-G)"},
    "ileS": {"query": "gene:ileS AND organism_id:1280", "desc": "Isoleucine--tRNA ligase"},
    "rpoB": {"query": "gene:rpoB AND organism_id:1280", "desc": "DNA-directed RNA polymerase subunit beta"},
    "murA": {"query": "gene:murA AND organism_id:1280", "desc": "UDP-GlcNAc 1-carboxyvinyltransferase (MurA)"},
    "rplD": {"query": "gene:rplD AND organism_id:1280", "desc": "50S ribosomal protein L4"},
    "rplC": {"query": "gene:rplC AND organism_id:1280", "desc": "50S ribosomal protein L3"},
    "rplV": {"query": "gene:rplV AND organism_id:1280", "desc": "50S ribosomal protein L22"},
    "rpsL": {"query": "gene:rpsL AND organism_id:1280", "desc": "30S ribosomal protein S12"},
    "rpsJ": {"query": "gene:rpsJ AND organism_id:1280", "desc": "30S ribosomal protein S10"},
}


@dataclass(frozen=True)
class TargetSpec:
    drug_class: str
    mechanism: str
    target_kind: str                     # "protein" | "membrane" | "cell_wall" | "unknown"
    target_genes: Tuple[str, ...]        # protein symbols in TARGET_PROTEINS to detect
    determinant_patterns: Tuple[str, ...]  # regex vs AMR gene symbols -> evidence (i)
    note: str = ""


# Regex fragments reused across drugs (case-insensitive, matched on gene symbol).
_MEC = r"^mec[ABC]"
_ERM = r"^erm"
_TET = r"^tet"
_DFR = r"^dfr"
_SUL = r"^sul[0-9]"
_CAT = r"^cat|^cml"
_AAC = r"aac\(6.\)|aph\(2|ant\(|aph\(3|aph\(9|arm[aA]|rmt|spc"


# --------------------------------------------------------------------------- #
# Drug → target spec. Keys are canonical drug ids (see io.labels.normalize_drug).
# --------------------------------------------------------------------------- #
DRUG_TARGETS: Dict[str, TargetSpec] = {
    # β-lactams -----------------------------------------------------------
    "penicillin": TargetSpec(
        "beta-lactam (penicillin)", "Inhibits PBP transpeptidation (cell wall).",
        "protein", ("pbpB", "pbpA"),
        (r"^blaZ", _MEC),
        "Penicillinase (blaZ) hydrolyses penicillin; mecA/mecC (PBP2a) bypasses the target."),
    "oxacillin": TargetSpec(
        "beta-lactam (anti-staph penicillin)", "Inhibits PBP transpeptidation.",
        "protein", ("pbpB", "pbpA"),
        (_MEC,),
        "Oxacillin/methicillin resistance is driven by mecA/mecC (acquired PBP2a)."),
    "methicillin": TargetSpec(
        "beta-lactam (anti-staph penicillin)", "Inhibits PBP transpeptidation.",
        "protein", ("pbpB", "pbpA"),
        (_MEC,),
        "Defines MRSA: mecA/mecC encode PBP2a, a low-affinity target bypass."),
    "cefoxitin": TargetSpec(
        "beta-lactam (cephamycin)", "Inhibits PBP transpeptidation; MRSA surrogate.",
        "protein", ("pbpB", "pbpA"),
        (_MEC,),
        "Standard phenotypic surrogate for mecA-mediated methicillin resistance."),
    # Fluoroquinolones ----------------------------------------------------
    "ciprofloxacin": TargetSpec(
        "fluoroquinolone", "Inhibits DNA gyrase and topoisomerase IV.",
        "protein", ("gyrA", "grlA"),
        (r"^gyrA", r"^grlA|^parC", r"^gyrB", r"^parE", r"^qnr", r"norA"),
        "High-level resistance from QRDR point mutations in gyrA + grlA(parC)."),
    "moxifloxacin": TargetSpec(
        "fluoroquinolone", "Inhibits DNA gyrase and topoisomerase IV.",
        "protein", ("gyrA", "grlA"),
        (r"^gyrA", r"^grlA|^parC", r"^gyrB", r"^parE"), ""),
    "levofloxacin": TargetSpec(
        "fluoroquinolone", "Inhibits DNA gyrase and topoisomerase IV.",
        "protein", ("gyrA", "grlA"),
        (r"^gyrA", r"^grlA|^parC", r"^gyrB", r"^parE"), ""),
    # Macrolide / lincosamide / phenicol (ribosome, 50S/23S) --------------
    "erythromycin": TargetSpec(
        "macrolide", "Binds 23S rRNA in the 50S subunit; blocks elongation.",
        "protein", ("rplD", "rplV"),
        (_ERM, r"^msr", r"^mph", r"^ere", r"^lsa"),
        "erm* rRNA methylases (MLS_B) and msr(A) efflux are the main determinants."),
    "clindamycin": TargetSpec(
        "lincosamide", "Binds 23S rRNA peptidyl-transferase centre (50S).",
        "protein", ("rplD", "rplC"),
        (_ERM, r"^lnu", r"^vga", r"^lsa", r"^sal", r"^cfr"),
        "erm* confers MLS_B (inducible/constitutive); lnu/vga/lsa also implicated."),
    "chloramphenicol": TargetSpec(
        "phenicol", "Binds 23S rRNA peptidyl-transferase centre (50S).",
        "protein", ("rplD", "rplV"),
        (_CAT, r"^cfr", r"^fex", r"^optrA"),
        "Chloramphenicol acetyltransferases (cat) inactivate the drug."),
    "linezolid": TargetSpec(
        "oxazolidinone", "Binds 23S rRNA; blocks initiation complex (50S).",
        "protein", ("rplC", "rplD"),
        (r"^cfr", r"^optrA", r"^poxtA", r"^rrl|23S"),
        "Rare in S. aureus; cfr/optrA/poxtA or 23S G2576T mutations."),
    # Aminoglycoside / tetracycline / glycylcycline (ribosome, 30S) -------
    "gentamicin": TargetSpec(
        "aminoglycoside", "Binds 16S rRNA in the 30S subunit; mistranslation.",
        "protein", ("rpsL", "rpsJ"),
        (_AAC,),
        "aac(6')-aph(2'') bifunctional enzyme is the dominant staph determinant."),
    "tetracycline": TargetSpec(
        "tetracycline", "Binds 30S subunit (16S rRNA); blocks tRNA binding.",
        "protein", ("rpsL", "rpsJ"),
        (_TET,),
        "tet(K)/tet(L) efflux and tet(M)/tet(O) ribosomal protection."),
    "tigecycline": TargetSpec(
        "glycylcycline", "Binds 30S subunit; higher affinity than tetracyclines.",
        "protein", ("rpsL", "rpsJ"),
        (r"^tet\(M\)", r"^mepA", r"^rpsJ"), ""),
    # Folate pathway ------------------------------------------------------
    "trimethoprim_sulfamethoxazole": TargetSpec(
        "folate-pathway inhibitor (combination)",
        "Trimethoprim inhibits DHFR (folA); sulfamethoxazole inhibits DHPS (folP).",
        "protein", ("folA", "folP"),
        (_DFR, _SUL, r"folA|F98Y"),
        "Two-target combination; dfr (DHFR bypass) drives clinical failure, sul supports."),
    "trimethoprim": TargetSpec(
        "folate-pathway inhibitor", "Inhibits dihydrofolate reductase (folA).",
        "protein", ("folA",),
        (_DFR, r"folA|F98Y"), ""),
    # Single-protein targets ---------------------------------------------
    "fusidic_acid": TargetSpec(
        "fusidane", "Inhibits elongation factor G (fusA); blocks translocation.",
        "protein", ("fusA",),
        (r"^fus[BCD]", r"^fusA"),
        "fusB/fusC protection proteins or fusA point mutations."),
    "mupirocin": TargetSpec(
        "pseudomonic acid", "Inhibits isoleucyl-tRNA synthetase (ileS).",
        "protein", ("ileS",),
        (r"^mup[AB]", r"^ileS"),
        "High-level resistance from mupA (ileS2); low-level from ileS mutations."),
    "rifampin": TargetSpec(
        "rifamycin", "Inhibits RNA polymerase β-subunit (rpoB).",
        "protein", ("rpoB",),
        (r"^rpoB|^arr", ),
        "Resistance from rpoB RRDR point mutations."),
    "fosfomycin": TargetSpec(
        "fosfomycin", "Inhibits MurA (UDP-GlcNAc enolpyruvyl transferase).",
        "protein", ("murA",),
        (r"^fos", r"^murA", r"^glpT|^uhpT"), ""),
    "phosphomycin": TargetSpec(  # alternate spelling seen in the TSV
        "fosfomycin", "Inhibits MurA (UDP-GlcNAc enolpyruvyl transferase).",
        "protein", ("murA",),
        (r"^fos", r"^murA"), ""),
    "ceftaroline": TargetSpec(
        "beta-lactam (anti-MRSA cephalosporin)",
        "Binds PBP2a as well as native PBPs; active against MRSA.",
        "protein", ("pbpB", "pbpA"),
        (r"^mecA", r"^pbp4|pbpA_"),
        "5th-gen cephalosporin; resistance from PBP2a/PBP4 target mutations."),
    "ampicillin": TargetSpec(
        "beta-lactam (aminopenicillin)", "Inhibits PBP transpeptidation.",
        "protein", ("pbpB", "pbpA"), (r"^blaZ", _MEC), ""),
    "amoxicillin_clavulanic_acid": TargetSpec(
        "beta-lactam + β-lactamase inhibitor", "Inhibits PBP transpeptidation.",
        "protein", ("pbpB", "pbpA"), (_MEC,),
        "Clavulanate restores activity vs blaZ; mecA still confers resistance."),
    "doxycycline": TargetSpec(
        "tetracycline", "Binds 30S subunit (16S rRNA); blocks tRNA binding.",
        "protein", ("rpsL", "rpsJ"), (_TET,), ""),
    "minocycline": TargetSpec(
        "tetracycline", "Binds 30S subunit (16S rRNA).",
        "protein", ("rpsL", "rpsJ"), (_TET,), ""),
    "tobramycin": TargetSpec(
        "aminoglycoside", "Binds 16S rRNA in the 30S subunit.",
        "protein", ("rpsL", "rpsJ"), (_AAC,), ""),
    "kanamycin": TargetSpec(
        "aminoglycoside", "Binds 16S rRNA in the 30S subunit.",
        "protein", ("rpsL", "rpsJ"), (_AAC,), ""),
    "streptomycin": TargetSpec(
        "aminoglycoside", "Binds 16S rRNA / S12 (rpsL) in the 30S subunit.",
        "protein", ("rpsL", "rpsJ"), (_AAC, r"^str[AB]", r"^ant\(6"), ""),
    "neomycin": TargetSpec(
        "aminoglycoside", "Binds 16S rRNA in the 30S subunit.",
        "protein", ("rpsL", "rpsJ"), (_AAC,), ""),
    # Membrane / cell-wall targets — no single-protein gate ---------------
    "daptomycin": TargetSpec(
        "lipopeptide", "Disrupts the cytoplasmic membrane (Ca2+-dependent).",
        "membrane", (),
        (r"^mprF", r"^cls", r"^pgsA", r"^walK", r"^rpoB"),
        "Target is the membrane, not one ORF; presence gate is not applicable."),
    "vancomycin": TargetSpec(
        "glycopeptide", "Binds D-Ala-D-Ala of lipid II (cell-wall precursor).",
        "cell_wall", (),
        (r"^van[ABCDEGH]", r"^walK", r"^graSR"),
        "VISA is a wall-thickening phenotype (no single gene); vanA/vanB rare in staph."),
    "teicoplanin": TargetSpec(
        "glycopeptide", "Binds D-Ala-D-Ala of lipid II (cell-wall precursor).",
        "cell_wall", (),
        (r"^van[ABCDEGH]", r"^walK"), ""),
}


def spec_for(drug: str) -> Optional[TargetSpec]:
    return DRUG_TARGETS.get(drug)
