"""Module 03 — The Decision Report: render a SamplePrediction for humans.

Every result is shown as: drug -> call, calibrated-style confidence, evidence
category, and the supporting genes/mutations. A mandatory "confirm with standard
lab testing" banner is always present. This module renders to plain text and
Markdown; the Streamlit app (``genome_firewall/app.py``) reuses the same objects.
"""

from __future__ import annotations

from typing import List

from .predict import (EV_KNOWN, EV_NONE, EV_STATISTICAL, NO_CALL, RESISTANT,
                      SUSCEPTIBLE, DrugPrediction, SamplePrediction)

_EVIDENCE_LABEL = {
    EV_KNOWN: "Known resistance gene/mutation detected",
    EV_STATISTICAL: "Statistical association only",
    EV_NONE: "No known resistance signal found",
}

_CALL_ICON = {RESISTANT: "✗", SUSCEPTIBLE: "✓", NO_CALL: "?"}


def confidence_band(conf: float) -> str:
    if conf >= 0.8:
        return "high"
    if conf >= 0.6:
        return "moderate"
    if conf >= 0.4:
        return "low"
    return "very low"


def render_text(sample: SamplePrediction, width: int = 78) -> str:
    """Plain-text report suitable for a terminal / log."""
    line = "=" * width
    out: List[str] = []
    out.append(line)
    out.append("GENOME FIREWALL — Antibiotic-Response Report (v0, zero-shot)".center(width))
    out.append(line)
    out.append(f"Species          : {sample.species}"
               f"{'' if sample.species_supported else '  [OUT OF SUPPORTED SCOPE]'}")
    out.append(f"Annotation       : {sample.annotation_backend} "
               f"(screening completeness {sample.screening_completeness:.2f})")
    if sample.qc:
        out.append(f"Assembly QC      : {'PASS' if sample.qc.get('passed') else 'FAIL'} "
                   f"({sample.qc.get('n_contigs')} contigs, "
                   f"{sample.qc.get('total_length')} bp)")
    out.append(f"No-call rate     : {sample.no_call_rate:.0%} "
               f"({sum(1 for p in sample.predictions if p.call == NO_CALL)}"
               f"/{len(sample.predictions)} drugs)")
    out.append("")

    for p in sample.predictions:
        out.append(_render_drug_text(p, width))
        out.append("-" * width)

    if sample.warnings:
        out.append("WARNINGS:")
        for w in sample.warnings:
            out.append(f"  ! {w}")
        out.append("")

    out.append("!! " + sample.safety_notice)
    out.append(line)
    return "\n".join(out)


def _render_drug_text(p: DrugPrediction, width: int) -> str:
    icon = _CALL_ICON.get(p.call, "?")
    head = f"[{icon}] {p.drug_name} ({p.drug_class})  ->  {p.call.upper()}"
    lines = [head]
    lines.append(f"      confidence : {p.confidence:.2f} ({confidence_band(p.confidence)}, "
                 f"rule-based/uncalibrated)")
    lines.append(f"      evidence   : {_EVIDENCE_LABEL.get(p.evidence_category, p.evidence_category)}")
    lines.append(f"      target     : {p.target_status}")
    if p.supporting_markers:
        genes = ", ".join(sorted({m['gene'] for m in p.supporting_markers}))
        lines.append(f"      markers    : {genes}")
    if p.no_call_reason:
        lines.append(f"      no-call because: {p.no_call_reason}")
    lines.append(f"      note       : {p.rationale}")
    return "\n".join(lines)


def render_markdown(sample: SamplePrediction) -> str:
    """Markdown report (used by CLI --format md and shareable outputs)."""
    md: List[str] = []
    md.append("# Genome Firewall — Antibiotic-Response Report")
    md.append("*v0 zero-shot rule-based prototype — research use only.*\n")
    scope = "" if sample.species_supported else " **⚠️ outside supported scope**"
    md.append(f"- **Species:** {sample.species}{scope}")
    md.append(f"- **Annotation backend:** {sample.annotation_backend} "
              f"(screening completeness {sample.screening_completeness:.2f})")
    if sample.qc:
        md.append(f"- **Assembly QC:** {'PASS ✅' if sample.qc.get('passed') else 'FAIL ❌'}")
    md.append(f"- **No-call rate:** {sample.no_call_rate:.0%}\n")

    md.append("| Antibiotic | Prediction | Confidence | Evidence | Markers |")
    md.append("|---|---|---|---|---|")
    for p in sample.predictions:
        genes = ", ".join(sorted({m['gene'] for m in p.supporting_markers})) or "—"
        md.append(f"| {p.drug_name} | **{p.call}** | {p.confidence:.2f} "
                  f"({confidence_band(p.confidence)}) | "
                  f"{_EVIDENCE_LABEL.get(p.evidence_category, p.evidence_category)} | {genes} |")
    md.append("")

    md.append("### Rationale per drug")
    for p in sample.predictions:
        md.append(f"- **{p.drug_name} → {p.call}:** {p.rationale}")
    md.append("")

    if sample.warnings:
        md.append("### Warnings")
        for w in sample.warnings:
            md.append(f"- {w}")
        md.append("")

    md.append("> ⚠️ **" + sample.safety_notice + "**")
    return "\n".join(md)
