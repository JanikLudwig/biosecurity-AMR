# Data Directory Structure

- `raw/`: Raw input data like original genome FASTA files. Do not modify these.
- `interim/`: Intermediate results, e.g. AMRFinderPlus TSV outputs.
- `processed/`: Final processed datasets like the binary feature matrices (`features.csv.gz`).

**Important Note**:
A missing AMR marker in these datasets does **NOT** prove susceptibility. The Genome Firewall is a research prototype, not a clinical decision support system.

**Git**:
- `raw/`, `interim/`, and `processed/` are ignored by git (see `.gitignore`). Do not commit large genome files or TSV collections.
- `tests/fixtures/` may contain small, synthetic data exclusively for unit testing.
