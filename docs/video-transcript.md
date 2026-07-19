# 60-Second Technical Video Transcript

The narration is 139 words. Read at a measured technical-demo pace while showing the three slides
specified in [`three-slide-design-prompt.md`](three-slide-design-prompt.md).

## Slide 1 — Feature extraction (0:00–0:18)

> GyraseX reads a *Staphylococcus aureus* genome. M1 turns known resistance genes and
> mutations into 158 auditable binary features. M2 independently verifies molecular targets using
> PyHMMER for proteins and BLASTN for ribosomal RNA. This prevents marker absence alone from
> triggering likely to work.

## Slide 2 — Three possible decisions (0:18–0:40)

> For genome 1280.51926, ciprofloxacin has no marker, verified targets, and nine-point-one
> percent resistance: likely to work. For 1280.9342, erythromycin has no marker and a
> verified target, but thirty-four-point-one percent falls inside the calibrated uncertainty band:
> no-call. In 1280.51872, GyrA S84L and ParC S80F support ninety-six-point-one percent: likely to
> fail.

## Slide 3 — Training without genome leakage (0:40–1:00)

> We train one regularized logistic model per drug on laboratory phenotypes. Homology-grouped train,
> calibration, and test splits reduce sequence leakage. Sigmoid calibration and error-constrained
> thresholds define the three decision regions. We measure balanced accuracy, class-specific recall,
> PR-AUC, Brier score, and no-call coverage. Every prediction still requires laboratory
> confirmation.

## Presenter notes

- Say **likely to work** and **likely to fail**, not simply susceptible and resistant.
- M1 marker absence is supporting evidence, never proof of susceptibility.
- M2 proves target presence; it is not itself a susceptibility prediction.
- The displayed evaluation is provisional grouped-development evidence, not clinical validation.
