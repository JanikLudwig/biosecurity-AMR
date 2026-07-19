from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

from genome_firewall.annotation.amrfinder import parse_output
from genome_firewall.config import DEFAULT_CONFIG, load_config
from genome_firewall.data.phenotypes import load_and_clean, summarize_labels

st.set_page_config(page_title="Genome Firewall Sandbox", layout="wide")


def render_table(frame: pd.DataFrame, *, max_rows: int = 250, height: int = 360) -> None:
    """Render without Streamlit's pandas-to-Arrow conversion, which can crash on reruns."""
    visible_rows = frame.head(max_rows)
    table_html = visible_rows.to_html(index=False, escape=True, border=0, classes="gf-table")
    st.markdown(
        f'<div class="gf-table-wrap" style="max-height:{height}px">{table_html}</div>',
        unsafe_allow_html=True,
    )
    if len(frame) > max_rows:
        st.caption(f"Showing the first {max_rows:,} of {len(frame):,} rows.")


st.markdown(
    """
    <style>
    .gf-table-wrap { overflow: auto; border: 1px solid rgba(128,128,128,.25); border-radius: .4rem; }
    .gf-table { width: 100%; border-collapse: collapse; font-size: .85rem; }
    .gf-table th { position: sticky; top: 0; background: var(--secondary-background-color); z-index: 1; }
    .gf-table th, .gf-table td { padding: .35rem .55rem; border-bottom: 1px solid rgba(128,128,128,.18); text-align: left; white-space: nowrap; }
    </style>
    """,
    unsafe_allow_html=True,
)
st.title("Genome Firewall data sandbox")
st.warning(
    "Research prototype only. Predictions must be confirmed by standard laboratory testing."
)

config = load_config(DEFAULT_CONFIG)
dataset = config["dataset"]


@st.cache_data(show_spinner="Loading laboratory phenotype data...")
def load_sandbox_data(
    source: str,
    source_mtime_ns: int,
    species: str,
    taxon_id: int,
    evidence: str,
    antibiotics: tuple[str, ...],
):
    """Read the source once; widget reruns reuse this cleaned dataset."""
    del source_mtime_ns  # cache key invalidates automatically when the source changes
    return load_and_clean(
        Path(source),
        species=species,
        taxon_id=taxon_id,
        evidence=evidence,
        antibiotics=antibiotics,
    )


source_path = Path(dataset["source_csv"])
result = load_sandbox_data(
    str(source_path),
    source_path.stat().st_mtime_ns,
    dataset["species"],
    dataset["taxon_id"],
    dataset["evidence"],
    tuple(dataset["antibiotics"]),
)
summary = summarize_labels(result.labels)

left, middle, right = st.columns(3)
left.metric("Usable genomes", f"{result.labels['genome_id'].nunique():,}")
middle.metric("Binary AST labels", f"{len(result.labels):,}")
right.metric("Configured drugs", len(dataset["antibiotics"]))

st.subheader("Laboratory label coverage")
long_counts = summary.melt(
    id_vars="antibiotic",
    value_vars=["Resistant", "Susceptible"],
    var_name="label",
    value_name="genomes",
)
figure = px.bar(
    long_counts,
    x="antibiotic",
    y="genomes",
    color="label",
    barmode="group",
    color_discrete_map={"Resistant": "#d1495b", "Susceptible": "#2a9d8f"},
)
st.plotly_chart(figure, width="stretch")
render_table(summary, height=260)

st.subheader("Cleaned observations")
drug = st.selectbox("Antibiotic", dataset["antibiotics"])
label = st.multiselect("Label", ["Resistant", "Susceptible"], default=["Resistant", "Susceptible"])
visible = result.labels.loc[
    result.labels["antibiotic"].eq(drug) & result.labels["label"].isin(label)
]
render_table(visible)

with st.expander("Exclusions and conflicts"):
    st.json(result.excluded_counts)
    render_table(result.conflicts, height=240)

amr_outputs = sorted(Path("data/interim/amrfinder").glob("*.tsv"))
st.subheader("AMRFinderPlus output browser")
if not amr_outputs:
    st.info("No AMRFinderPlus TSV files exist yet. This panel will populate after FASTA annotation.")
else:
    selected = st.selectbox("AMRFinder result", amr_outputs, format_func=lambda path: path.name)
    raw_amr = pd.read_csv(selected, sep="\t", dtype=object, keep_default_na=False)
    st.caption("Raw AMRFinderPlus TSV")
    render_table(raw_amr, max_rows=100)
    st.caption("Normalized evidence used by Genome Firewall")
    normalized = parse_output(selected, genome_id=selected.stem)
    render_table(normalized, max_rows=100)

st.divider()
st.header("Development pipeline outputs")

features_path = Path("data/processed/features/amr-features.csv")
splits_path = Path("data/processed/splits-500/genome-splits.csv")
if not splits_path.is_file():
    splits_path = Path("data/processed/splits/genome-splits.csv")
if features_path.is_file():
    feature_matrix = pd.read_csv(features_path, dtype=object, keep_default_na=False)
    feature_columns = [column for column in feature_matrix if column != "genome_id"]
    prevalence = (
        feature_matrix[feature_columns]
        .sum()
        .sort_values(ascending=False)
        .head(20)
        .rename_axis("feature")
        .reset_index(name="genomes")
    )
    st.subheader("Most common AMRFinder features in the development cohort")
    st.plotly_chart(
        px.bar(prevalence, x="genomes", y="feature", orientation="h").update_layout(
            yaxis={"categoryorder": "total ascending"}
        ),
        width="stretch",
    )
else:
    st.info("Run the AMRFinder batch command to create a cohort feature matrix.")

if splits_path.is_file():
    splits = pd.read_csv(splits_path, dtype=object, keep_default_na=False)
    split_summary = (
        splits.groupby("split")
        .agg(genomes=("genome_id", "size"), clusters=("cluster_id", "nunique"))
        .reset_index()
    )
    st.subheader("Homology-aware dataset split")
    split_left, split_right = st.columns([1, 2])
    with split_left:
        render_table(split_summary, height=220)
    cluster_summary = (
        splits.groupby(["split", "cluster_id"])["genome_id"]
        .size()
        .rename("genomes")
        .reset_index()
    )
    split_right.plotly_chart(
        px.bar(cluster_summary, x="cluster_id", y="genomes", color="split"),
        width="stretch",
    )
    st.caption(
        "Each genetic-similarity cluster belongs to exactly one partition; no cluster crosses "
        "training, calibration, and test."
    )

model_summary_path = Path("artifacts/models/model-summary.csv")
predictions_path = Path("artifacts/models/test-predictions.csv")
reliability_path = Path("artifacts/models/test-reliability.csv")
if model_summary_path.is_file():
    model_summary = pd.read_csv(model_summary_path, dtype=object, keep_default_na=False)
    st.subheader("Baseline model readiness")
    readiness_columns = [
        "antibiotic",
        "status",
        "calibration_status",
        "no_call_rate",
        "class_counts.train.susceptible",
        "class_counts.train.resistant",
        "class_counts.calibration.susceptible",
        "class_counts.calibration.resistant",
        "class_counts.test.susceptible",
        "class_counts.test.resistant",
        "called_test_metrics.coverage",
        "called_test_metrics.accuracy",
        "called_test_metrics.resistant_fraction_in_susceptible_calls",
        "called_test_metrics.susceptible_fraction_in_resistant_calls",
    ]
    readiness_columns = [column for column in readiness_columns if column in model_summary]
    render_table(model_summary[readiness_columns], height=300)
    unavailable = ~model_summary["calibration_status"].str.startswith("sigmoid")
    if unavailable.any():
        st.warning(
            "At least one development model could not establish a safe calibration boundary; "
            "its decisions remain no-call."
        )
    partial = model_summary["calibration_status"].eq("sigmoid_partial_thresholds")
    if partial.any():
        drugs = ", ".join(model_summary.loc[partial, "antibiotic"])
        st.info(
            f"{drugs}: only one call direction met the configured calibration-error limit; "
            "the other direction remains no-call."
        )

if predictions_path.is_file():
    predictions = pd.read_csv(predictions_path, dtype=object, keep_default_na=False)
    predictions["probability_resistant"] = pd.to_numeric(
        predictions["probability_resistant"], errors="coerce"
    )
    model_drug = st.selectbox(
        "Model diagnostics antibiotic",
        dataset["antibiotics"],
        key="model_diagnostics_antibiotic",
    )
    selected_predictions = predictions.loc[predictions["antibiotic"].eq(model_drug)]
    st.plotly_chart(
        px.histogram(
            selected_predictions,
            x="probability_resistant",
            color="label",
            nbins=10,
            range_x=[0, 1],
        ),
        width="stretch",
    )
    if reliability_path.is_file():
        reliability = pd.read_csv(reliability_path)
        selected_reliability = reliability.loc[reliability["antibiotic"].eq(model_drug)]
        reliability_figure = px.scatter(
            selected_reliability,
            x="mean_probability_resistant",
            y="observed_resistant_fraction",
            size="samples",
            range_x=[0, 1],
            range_y=[0, 1],
            title="Held-out reliability",
        )
        reliability_figure.add_shape(
            type="line", x0=0, y0=0, x1=1, y1=1, line={"dash": "dash"}
        )
        st.plotly_chart(reliability_figure, width="stretch")
    render_table(selected_predictions)
