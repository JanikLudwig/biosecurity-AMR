"""Central paths, constants, and tunable thresholds for the pipeline.

Everything that a reviewer might want to re-tune (decision thresholds, the
minimum lab evidence needed to model a drug, split proportions) lives here so
the engine code stays declarative.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Tuple

# --------------------------------------------------------------------------- #
# Filesystem layout (all relative to the repository root).
# --------------------------------------------------------------------------- #
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

GENOMES_DIR = os.path.join(ROOT, "genomes")
BVBRC_DIR = os.path.join(ROOT, "bvbrc_data")

DATA_DIR = os.path.join(ROOT, "data")
REFERENCES_DIR = os.path.join(DATA_DIR, "references", "targets")
AMRFINDER_DIR = os.path.join(DATA_DIR, "amrfinder")     # teammates drop TSVs here
ARTIFACTS_DIR = os.path.join(DATA_DIR, "artifacts")     # feature matrix, target calls, split

MODELS_DIR = os.path.join(ROOT, "models")
REPORTS_DIR = os.path.join(ROOT, "reports")

# Pre-defined inputs shipped in the repo.
LABELS_TSV = os.path.join(BVBRC_DIR, "aureus_all_amr_labels.tsv")
MANIFEST_CSV = os.path.join(BVBRC_DIR, "genome_manifest.csv")
METADATA_CSV = os.path.join(BVBRC_DIR, "aureus_genome_metadata.csv")

# Derived artifacts (produced by scripts/).
FEATURES_PARQUET = os.path.join(ARTIFACTS_DIR, "features.parquet")
TARGET_CALLS_PARQUET = os.path.join(ARTIFACTS_DIR, "target_calls.parquet")
SPLIT_CSV = os.path.join(ARTIFACTS_DIR, "split.csv")
PANEL_JSON = os.path.join(ARTIFACTS_DIR, "panel.json")


def ensure_dirs() -> None:
    for d in (DATA_DIR, REFERENCES_DIR, AMRFINDER_DIR, ARTIFACTS_DIR,
              MODELS_DIR, REPORTS_DIR):
        os.makedirs(d, exist_ok=True)


# --------------------------------------------------------------------------- #
# Panel construction (data-driven from the lab TSV).
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class PanelConfig:
    """How many laboratory-measured labels a drug needs to be modelled."""

    tier_a_min_per_class: int = 100   # both R and S ≥ this  -> train + calibrate
    tier_b_min_per_class: int = 20    # 20..100 minority     -> low-power, flagged
    # below tier_b_min_per_class in the minority class      -> Tier C structural no-call
    min_groups: int = 8               # need ≥ this many MLST groups to split honestly


# --------------------------------------------------------------------------- #
# Leakage-safe split.
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class SplitConfig:
    # Proportions of MLST groups (weighted by genome count) for each partition.
    train: float = 0.70
    calibration: float = 0.15
    test: float = 0.15
    seed: int = 1280
    dedup_by: str = "hc10"            # collapse near-identical genomes first
    group_by: str = "mlst"            # whole MLST lineages held out at test


# --------------------------------------------------------------------------- #
# Assembly QC gate (uses the CheckM columns already in the metadata).
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class QCConfig:
    min_completeness: float = 90.0
    max_contamination: float = 5.0
    max_contigs: int = 300
    min_length: int = 2_400_000
    max_length: int = 3_300_000


# --------------------------------------------------------------------------- #
# Decision thresholds (M4). Confidence is the calibrated probability of the
# called class; these bands convert it into likely-to-fail/work/no-call.
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class DecisionConfig:
    p_fail_hi: float = 0.65           # p(R) ≥ hi           -> likely to fail
    p_work_lo: float = 0.35           # p(R) ≤ lo (+target) -> likely to work
    # The band (lo, hi) is the honest no-call zone.
    tier_b_nocall_widen: float = 0.10  # push Tier-B drugs further toward no-call
    # A "likely to work" call is only allowed when the target is provably present.
    require_target_for_work: bool = True


PANEL = PanelConfig()
SPLIT = SplitConfig()
QC = QCConfig()
DECISION = DecisionConfig()
