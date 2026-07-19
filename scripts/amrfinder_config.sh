# ============================================================================
# amrfinder_config.sh — single source of truth for the AMRFinderPlus run.
# Sourced by amrfinder_prep.sh, amrfinder_array.sh, amrfinder_merge.sh so the
# flags RECORDED in provenance are guaranteed identical to the flags RUN.
# ============================================================================

ROOT=/work3/janlud/biosecurity-AMR
GENOME_DIR=$ROOT/genomes
LIST=$ROOT/genomes.list                 # one absolute genome path per line
OUTDIR=$ROOT/data/amrfinder             # per-genome <genome_id>.tsv land here
LOGDIR=$ROOT/logs
PROV=$OUTDIR/RUN_PROVENANCE.txt         # the reproducibility record for teammates

# --- conda env holding ncbi-amrfinderplus (edit if you installed elsewhere) ---
CONDA_SH=/work3/janlud/miniforge3/etc/profile.d/conda.sh
CONDA_ENV=amrfinder

# --- REPRODUCIBILITY-CRITICAL knobs. Changing any of these = everyone re-runs ---
ORGANISM=Staphylococcus_aureus          # enables S. aureus point-mutation screen
THREADS=4                               # does NOT affect results, only speed
# Exact per-genome invocation. Both prep (for the record) and array (for the run)
# call THIS function, so the documented command == the executed command.
run_amrfinder() {   # $1=input.fna  $2=genome_id  $3=output.tsv
    amrfinder -n "$1" --name "$2" -o "$3" \
        --organism "$ORGANISM" --plus --threads "$THREADS"
}
# Human-readable flag string for the provenance file (keep in sync with above).
AMRFINDER_FLAGS="-n <genome.fna> --name <genome_id> --organism ${ORGANISM} --plus --threads ${THREADS}"
