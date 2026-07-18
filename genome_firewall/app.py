"""Streamlit demo for Genome Firewall v0 (Module 03 — Decision Report UI).

Run with:
    streamlit run genome_firewall/app.py

Upload a reconstructed assembly (FASTA) or a precomputed AMR table (TSV), pick a
backend, and get a per-antibiotic report with calibrated-style confidence, an
evidence category, an explicit no-call, and a mandatory lab-confirmation banner.
"""

from __future__ import annotations

import os
import tempfile

import streamlit as st

from genome_firewall import SAFETY_NOTICE, __version__
from genome_firewall import knowledge as kb
from genome_firewall.annotate import annotate
from genome_firewall.fasta import compute_qc
from genome_firewall.predict import NO_CALL, RESISTANT, SUSCEPTIBLE, predict_sample
from genome_firewall.report import confidence_band

_HERE = os.path.dirname(__file__)
_EXAMPLES = {
    "E. coli — multidrug-resistant (example)":
        os.path.join(_HERE, "examples", "ecoli_resistant_amrfinder.tsv"),
    "E. coli — susceptible / no markers (example)":
        os.path.join(_HERE, "examples", "ecoli_susceptible_amrfinder.tsv"),
    "E. coli — weak/ambiguous markers (example)":
        os.path.join(_HERE, "examples", "ecoli_weak_amrfinder.tsv"),
}

_CALL_STYLE = {
    RESISTANT: ("#b3261e", "✗", "LIKELY TO FAIL"),
    SUSCEPTIBLE: ("#1e7d32", "✓", "LIKELY TO WORK"),
    NO_CALL: ("#8a6d00", "?", "NO-CALL"),
}


def _run(fasta_path, tsv_path, backend, tsv_source, species, organism):
    qc = compute_qc(fasta_path).as_dict() if fasta_path else None
    annotation = annotate(fasta_path=fasta_path, backend=backend, tsv_path=tsv_path,
                          organism=organism or None, tsv_source=tsv_source)
    return predict_sample(annotation, species=species or None, qc=qc), annotation


def main() -> None:
    st.set_page_config(page_title="Genome Firewall v0", page_icon="🧬", layout="wide")

    # ---- persistent safety banner (mandatory) ----
    st.warning("⚠️ " + SAFETY_NOTICE, icon="⚠️")

    st.title("🧬 Genome Firewall")
    st.caption(f"v{__version__} · zero-shot, rule-based AMR decision support · "
               "strictly defensive research prototype")

    with st.sidebar:
        st.header("About")
        st.markdown(
            "Turns a **reconstructed bacterial genome** into a per-antibiotic "
            "prediction — *likely to fail / likely to work / no-call* — from known "
            "resistance genes and mutations. **No machine-learning training** (v0)."
        )
        st.subheader("Supported scope")
        st.markdown(f"**Species:** {kb.supported_species()}")
        st.markdown("**Antibiotics:** " +
                    ", ".join(d["name"] for d in kb.panel_drugs()))
        st.info("Anything outside this scope is returned as **no-call** — the tool "
                "does not guess beyond what it covers.")
        st.subheader("Responsibility")
        st.markdown(
            "- Defensive by construction — predicts existing resistance only.\n"
            "- Honest evidence — separates *known determinant* from *statistical*.\n"
            "- Calibrated-style confidence + explicit **no-call**.\n"
            "- Human oversight required — confirm with the lab."
        )

        st.header("Input")
        mode = st.radio("Backend", ["tsv", "auto", "amrfinderplus", "camrah"], index=0,
                        help="'tsv' ingests a precomputed AMR table and needs no "
                             "bioinformatics install. 'auto' uses cAMRah/AMRFinderPlus "
                             "if installed.")
        example_choice = st.selectbox("Load an example (precomputed TSV)",
                                      ["— none —"] + list(_EXAMPLES.keys()))
        species = st.text_input("Species", value="Escherichia coli")
        organism = st.text_input("AMRFinderPlus --organism (optional)", value="Escherichia")
        tsv_source = st.selectbox("Precomputed table source",
                                  ["amrfinderplus", "camrah", "abricate", "resfinder"], index=0)

    up_fasta = st.file_uploader("Reconstructed assembly (FASTA / .gz)",
                                type=["fasta", "fa", "fna", "gz"])
    up_tsv = st.file_uploader("Precomputed AMR table (TSV)", type=["tsv", "txt", "csv"])

    run = st.button("Predict antibiotic response", type="primary")
    if not run:
        st.stop()

    fasta_path = tsv_path = None
    tmp_paths = []
    try:
        if up_fasta is not None:
            suffix = ".gz" if up_fasta.name.endswith(".gz") else ".fasta"
            fasta_path = _spool(up_fasta.getvalue(), suffix)
            tmp_paths.append(fasta_path)
        if up_tsv is not None:
            tsv_path = _spool(up_tsv.getvalue(), ".tsv")
            tmp_paths.append(tsv_path)
        if tsv_path is None and example_choice != "— none —":
            tsv_path = _EXAMPLES[example_choice]

        if fasta_path is None and tsv_path is None:
            st.error("Upload a FASTA and/or a TSV, or pick an example.")
            st.stop()

        with st.spinner("Annotating and predicting…"):
            sample, annotation = _run(fasta_path, tsv_path,
                                      backend=mode, tsv_source=tsv_source,
                                      species=species, organism=organism)
        _render(sample, annotation)
    finally:
        for p in tmp_paths:
            try:
                os.unlink(p)
            except OSError:
                pass


def _spool(data: bytes, suffix: str) -> str:
    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp:
        tmp.write(data)
        return tmp.name


def _render(sample, annotation) -> None:
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Species", sample.species,
              "in scope" if sample.species_supported else "OUT OF SCOPE")
    c2.metric("Backend", sample.annotation_backend)
    c3.metric("AMR determinants", len(annotation.hits))
    c4.metric("No-call rate", f"{sample.no_call_rate:.0%}")

    if sample.qc:
        with st.expander("Assembly QC"):
            st.json(sample.qc)
    if sample.warnings:
        for w in sample.warnings:
            st.warning(w)

    st.subheader("Antibiotic-response report")
    for p in sample.predictions:
        color, icon, label = _CALL_STYLE[p.call]
        with st.container(border=True):
            top = st.columns([3, 2, 2, 3])
            top[0].markdown(f"### {p.drug_name}")
            top[0].caption(p.drug_class)
            top[1].markdown(
                f"<span style='color:{color};font-weight:700;font-size:1.1rem'>"
                f"{icon} {label}</span>", unsafe_allow_html=True)
            top[2].metric("Confidence", f"{p.confidence:.2f}",
                          confidence_band(p.confidence))
            top[3].caption("Evidence"); top[3].write(_evidence_label(p.evidence_category))
            st.progress(min(1.0, max(0.0, p.confidence)))
            if p.supporting_markers:
                genes = ", ".join(sorted({m["gene"] for m in p.supporting_markers}))
                st.markdown(f"**Supporting markers:** `{genes}`")
                with st.expander("Marker details"):
                    st.dataframe(p.supporting_markers, width="stretch")
            st.caption(("🎯 target: " + p.target_status) +
                       (f"  ·  ⚠️ no-call: {p.no_call_reason}" if p.no_call_reason else ""))
            st.write(p.rationale)

    st.divider()
    st.error("⚠️ " + sample.safety_notice)


def _evidence_label(cat: str) -> str:
    return {
        "known_resistance_determinant": "🧬 Known resistance gene/mutation",
        "statistical_association": "📊 Statistical association only",
        "no_known_resistance_signal": "— No known resistance signal",
    }.get(cat, cat)


if __name__ == "__main__":
    main()
