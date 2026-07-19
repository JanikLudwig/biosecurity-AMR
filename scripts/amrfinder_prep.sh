#!/bin/sh
# ============================================================================
# amrfinder_prep.sh — run ONCE after installing AMRFinderPlus, before submitting.
# Builds the genome list, then freezes the four reproducibility facts
# (genomes / software version / database version / flags) into RUN_PROVENANCE.txt.
# Prints the exact block to paste to teammates.
#   Usage:  sh scripts/amrfinder_prep.sh
# ============================================================================
set -eu
cd "$(dirname "$0")/.."
. scripts/amrfinder_config.sh
. "$CONDA_SH"; conda activate "$CONDA_ENV"

mkdir -p "$OUTDIR" "$LOGDIR"

# 1) Freeze the input set (sorted, absolute paths) + a portable set fingerprint.
find "$GENOME_DIR" -name '*.fna' | sort > "$LIST"
N=$(wc -l < "$LIST")
SETHASH=$(sed 's#.*/##; s#\.fna$##' "$LIST" | sort | sha256sum | cut -d' ' -f1)

# 2) Capture software + database version (authoritative: amrfinder --version output).
VER=$(amrfinder --version 2>&1 || echo "UNKNOWN")
# DB version: read the database's own version file, with fallbacks.
DBVER=$(cat "$CONDA_PREFIX"/share/amrfinderplus/data/latest/version.txt 2>/dev/null \
        || ls "$CONDA_PREFIX"/share/amrfinderplus/data 2>/dev/null | grep -E '^[0-9]{4}-' | sort | tail -1 \
        || echo "UNKNOWN — run 'amrfinder -u' first")

# 3) Write the provenance record.
{
  echo "=========================================================="
  echo " AMRFinderPlus run — reproducibility record"
  echo "=========================================================="
  echo "generated_utc : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host          : $(hostname)"
  echo "operator      : $(whoami)"
  echo "----------------------------------------------------------"
  echo "GENOMES       : $N assemblies"
  echo "genome_dir    : $GENOME_DIR"
  echo "set_sha256    : $SETHASH   (sha256 of sorted genome_ids)"
  echo "list_file     : $LIST"
  echo "----------------------------------------------------------"
  echo "SOFTWARE      :"
  echo "$VER" | sed 's/^/    /'
  echo "DATABASE_VER  : $DBVER"
  echo "----------------------------------------------------------"
  echo "FLAGS         : $AMRFINDER_FLAGS"
  echo "----------------------------------------------------------"
  echo "conda pins    :"
  conda list -n "$CONDA_ENV" -e 2>/dev/null | grep -Ei 'amrfinder|blast|hmmer|prodigal' | sed 's/^/    /'
  echo "=========================================================="
} | tee "$PROV"

echo
echo ">>> Provenance written to: $PROV"
echo ">>> Send that block to teammates. Now submit with:"
echo "        bsub < scripts/amrfinder_array.sh"
