#!/bin/sh
### AMRFinderPlus over all genomes as an LSF job array (DTU HPC).
### Run scripts/amrfinder_prep.sh first, then:   bsub < scripts/amrfinder_array.sh
### The array size in -J MUST equal NTASKS below.
#BSUB -J amrfinder[1-99]%40       # 99 tasks, at most 40 running at once (polite on shared 'hpc')
#BSUB -q hpc
#BSUB -n 4                        # cores per task (matches THREADS in config)
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2GB]"        # per-core => 8GB reserved per task
#BSUB -W 2:00                     # walltime per task (generous; ~35 min real at 46 genomes)
#BSUB -o logs/amr_%J_%I.out       # %J=jobID  %I=array index
#BSUB -e logs/amr_%J_%I.err

set -eu
cd "$LS_SUBCWD" 2>/dev/null || cd "$(dirname "$0")/.."
. scripts/amrfinder_config.sh
. "$CONDA_SH"; conda activate "$CONDA_ENV"
mkdir -p "$OUTDIR" "$LOGDIR"

NTASKS=99                          # MUST match [1-99] in -J above
FAILLOG=$LOGDIR/failures_${LSB_JOBID}_${LSB_JOBINDEX}.txt

# ---- this task's contiguous slice of the genome list ----
total=$(wc -l < "$LIST")
per=$(( (total + NTASKS - 1) / NTASKS ))          # ceil(total / NTASKS)
start=$(( (LSB_JOBINDEX - 1) * per + 1 ))
end=$(( LSB_JOBINDEX * per ))
echo "[task $LSB_JOBINDEX/$NTASKS] genomes $start..$end of $total"

# Materialize this task's slice to a file so the loop runs in THIS shell
# (a `sed | while` pipe would run the loop in a subshell and lose the counters).
SLICE=$LOGDIR/slice_${LSB_JOBID}_${LSB_JOBINDEX}.txt
sed -n "${start},${end}p" "$LIST" > "$SLICE"

done_ct=0; skip_ct=0; fail_ct=0
while IFS= read -r fna; do
    [ -z "$fna" ] && continue
    name=$(basename "$fna" .fna)
    out=$OUTDIR/$name.tsv
    [ -s "$out" ] && { skip_ct=$((skip_ct+1)); continue; }   # resume: already complete
    # Write to a temp file and rename only on success => no partial files survive a kill,
    # and one genome's failure never aborts the task (set -e is scoped out by the if).
    if run_amrfinder "$fna" "$name" "$out.tmp.$$" >>"$LOGDIR/amr_${LSB_JOBID}_${LSB_JOBINDEX}.log" 2>&1; then
        mv "$out.tmp.$$" "$out"
        done_ct=$((done_ct+1))
    else
        rm -f "$out.tmp.$$"
        echo "$name" >> "$FAILLOG"
        fail_ct=$((fail_ct+1))
    fi
done < "$SLICE"
rm -f "$SLICE"

echo "[task $LSB_JOBINDEX] done: $done_ct new, $skip_ct skipped, $fail_ct failed"
[ -s "$FAILLOG" ] && echo "[task $LSB_JOBINDEX] failures listed in $FAILLOG" || true
