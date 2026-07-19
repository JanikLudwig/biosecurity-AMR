# 60-Second Video Transcript

## Spoken narration

Antibiotic resistance testing can take days. Genome Firewall is a defensive research prototype
that analyzes an assembled *Staphylococcus aureus* genome earlier.

First, M1 runs AMRFinderPlus and converts known resistance genes and mutations into 158 auditable
features. Calibrated models estimate resistance for cefoxitin, ciprofloxacin, and erythromycin.

But missing resistance markers do not prove a drug will work. M2 calls the genome's proteins and
uses PyHMMER to verify protein targets, while a nucleotide search verifies erythromycin's 23S RNA
target.

We group closely related genomes before splitting them into training, calibration, and untouched
test sets, reducing sequence leakage. We evaluate both resistant and susceptible recall,
class-imbalance metrics, probability calibration, and performance after no-calls.

The report separates known biological evidence from statistical association and returns no-call
for uncertain, novel, conflicting, or target-unverified cases. Every result must still be confirmed
by standard laboratory testing.

## Suggested visuals

- **0:00–0:08:** FASTA upload and the defensive-use label.
- **0:08–0:25:** M1 AMRFinder genes/mutations flowing into three probability cards.
- **0:25–0:38:** M2 PyHMMER protein hits and the 23S RNA target gate.
- **0:38–0:50:** grouped train/calibration/test split and compact reliability/no-call metrics.
- **0:50–1:00:** final decision table, evidence categories, and laboratory-confirmation warning.
