"""Genome Firewall (v0, zero-shot) — a defensive AMR decision-support prototype.

Turns a reconstructed bacterial genome (FASTA) into a per-antibiotic prediction
(*likely to fail* / *likely to work* / *no-call*) using a transparent, rule-based
("zero-shot") knowledge base — no model training, no deep learning.

Pipeline (see the challenge brief modules):
    Module 01  genome_firewall.annotate  — FASTA -> AMR gene/mutation features
    Module 02  genome_firewall.predict    — features -> per-drug prediction + gate
    Module 03  genome_firewall.report     — prediction -> auditable decision report

This is a research prototype. Every result must be confirmed with standard
laboratory testing. The system is strictly defensive: it only predicts and
explains resistance that already exists. It never designs or modifies organisms.
"""

__version__ = "0.1.0"

SAFETY_NOTICE = (
    "Research prototype — decision support only. Every antibiotic-response result "
    "MUST be confirmed with standard laboratory testing before any treatment "
    "decision. This tool predicts and explains existing resistance to support "
    "clinicians and public-health tracking; it never designs or modifies organisms."
)
