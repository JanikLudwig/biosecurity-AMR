#!/bin/sh
# ============================================================================
# amrfinder_merge.sh — run AFTER the array job completes.
#   1. Completeness check: every genome in the list produced a .tsv?
#   2. Fold per-genome TSVs -> real feature matrix (gfw.m1_adapter, synthetic=False).
#   3. Emit a combined TSV and append completion stats to the provenance record.
#   Usage:  sh scripts/amrfinder_merge.sh
# ============================================================================
set -eu
cd "$(dirname "$0")/.."
. scripts/amrfinder_config.sh
. "$CONDA_SH"; conda activate "$CONDA_ENV"

expected=$(wc -l < "$LIST")
produced=$(find "$OUTDIR" -maxdepth 1 -name '*.tsv' ! -name 'RUN_PROVENANCE*' | wc -l)
echo "completeness: $produced / $expected genomes have output"

# List any genomes still missing an output file.
MISSING=$LOGDIR/missing_genomes.txt
: > "$MISSING"
while IFS= read -r fna; do
    name=$(basename "$fna" .fna)
    [ -s "$OUTDIR/$name.tsv" ] || echo "$name" >> "$MISSING"
done < "$LIST"
nmiss=$(wc -l < "$MISSING")
if [ "$nmiss" -gt 0 ]; then
    echo "WARNING: $nmiss genomes missing — see $MISSING"
    echo "  Re-submit to fill gaps (resume skips finished): bsub < scripts/amrfinder_array.sh"
else
    echo "all genomes accounted for."
    rm -f "$MISSING"
fi

# Fold into the real feature matrix (replaces the __synthetic__ placeholder parquet).
python - <<'PY'
from gfw.m1_adapter import fold_amrfinder_dir, save_features, FEATURES_PARQUET
mat = fold_amrfinder_dir()
if mat.empty:
    raise SystemExit("no AMRFinder TSVs parsed — nothing to fold")
save_features(mat, FEATURES_PARQUET, synthetic=False)
print(f"real feature matrix: {mat.shape[0]} genomes x {mat.shape[1]} gene/mutation symbols")
print(f"  -> {FEATURES_PARQUET}  (__synthetic__=0)")
PY

# Combined TSV for human inspection (header from first file + all data rows w/ a genome_id col).
COMBINED=$OUTDIR/../amrfinder_all.tsv
first=$(find "$OUTDIR" -maxdepth 1 -name '*.tsv' ! -name 'RUN_PROVENANCE*' | sort | head -1)
if [ -n "${first:-}" ]; then
    { printf 'genome_id\t'; head -1 "$first"; } > "$COMBINED"
    for f in "$OUTDIR"/*.tsv; do
        case "$f" in *RUN_PROVENANCE*) continue;; esac
        gid=$(basename "$f" .tsv)
        tail -n +2 "$f" | sed "s/^/$gid\t/"
    done >> "$COMBINED"
    echo "combined table -> $COMBINED"
fi

{ echo "----------------------------------------------------------"
  echo "COMPLETED_utc : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "produced      : $produced / $expected genomes"
  echo "missing       : $nmiss"
} >> "$PROV"
echo "provenance updated: $PROV"
