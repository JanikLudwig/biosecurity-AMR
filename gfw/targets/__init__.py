"""M2 — the Drug-Target Detector.

Proves, deterministically, that an antibiotic's molecular target is present in a
genome, so the decision layer may only say *likely to work* when the drug
actually has something to act on — never from the mere absence of resistance
markers (`# design notes.md`).
"""

from .specs import DRUG_TARGETS, TARGET_PROTEINS, TargetSpec, spec_for  # noqa: F401
