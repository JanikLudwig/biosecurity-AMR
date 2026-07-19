"""Genome Firewall — an AI defense system against superbugs (S. aureus prototype).

A deterministic, auditable pipeline that turns a reconstructed *Staphylococcus
aureus* genome into an *earlier* antibiotic-response prediction — one of
**likely to fail / likely to work / no-call** per drug — with a calibrated
confidence, an evidence category, and a mandatory lab-confirmation notice.

Module map (see ``plans`` / README):

* **M1** Genome Reader — *teammate-owned* AMRFinderPlus feature extractor;
  consumed here via :mod:`gfw.m1_adapter`.
* **M2** Drug-Target Detector — :mod:`gfw.targets`; proves the drug's molecular
  target is present in the genome so *likely to work* is never asserted from the
  mere absence of resistance markers.
* **M3** Predictor — :mod:`gfw.predict`; one calibrated logistic regression per
  antibiotic on the AMR feature matrix.
* **M4** Decision layer — :mod:`gfw.decide`; fuses the model probability, the
  target gate, and the determinant evidence into the three-way call.
* **M5** Report — :mod:`gfw.report`; JSON / Markdown decision report.

This tool is **strictly defensive**: it only predicts and explains resistance
that already exists. It never designs, modifies, or optimizes an organism.
"""

from __future__ import annotations

__version__ = "1.0.0-saureus"

SAFETY_NOTICE = (
    "Research prototype — decision support only. Every antibiotic-response "
    "result MUST be confirmed with standard laboratory testing before any "
    "treatment decision. This system predicts and explains resistance that "
    "already exists; it never designs, modifies, or optimizes an organism."
)

SPECIES = "Staphylococcus aureus"
SPECIES_TAXON = 1280
