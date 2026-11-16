#!/bin/bash
#SBATCH --job-name=plk-timing
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#
# Step 2: every config for both builds in ONE job on ONE node, one process at a time, so node
# to node variance (about 2x between identical runs) stays out of the comparison. All logic
# is in benchmarking/scripts/time-per-length-k.py (see its docstring). Generous wall: the job
# ends when it finishes, and each process is capped (--cap-min) inside the script.

set -euo pipefail
REPO="$HOME/projects/pepmatch/PEPMatch2.0"
: "${RUN_DIR:?submit through submit-plk-pipeline.sh, which sets RUN_DIR}"
EXTRA="${EXTRA:-}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-1},8}"

echo "=== per-length k vs one k per batch ==="
echo "job     : ${SLURM_JOB_ID:-manual}"
echo "host    : $(hostname)"
echo "run dir : $RUN_DIR"
echo "started : $(date -Is)"
source "$HOME/projects/pepmatch/env.sh"
# The environment's own pepmatch must not be what gets timed: PYTHONPATH is set per process
# by the script, and each worker checks where pepmatch came from.
unset PYTHONPATH

SCRATCH="/mnt/hpcscratch/ryoung/per-length-k-${SLURM_JOB_ID:-manual}"
mkdir -p "$SCRATCH"
cp "$REPO/benchmarking/proteomes/human.fasta" "$SCRATCH/"

python3 "$REPO/benchmarking/scripts/time-per-length-k.py" \
  --workdir "$SCRATCH" --outdir "$RUN_DIR" --proteome "$SCRATCH/human.fasta" \
  --base-site "$RUN_DIR/site-base" --head-site "$RUN_DIR/site-head" \
  --base-label "$(cat "$RUN_DIR/base-commit.txt")" --head-label "$(cat "$RUN_DIR/head-commit.txt")" \
  --threads "$THREADS" $EXTRA

echo
echo "=== paste these back to Kairo ==="
for f in run-info.txt per-length-k-summary.tsv; do
  [[ -f "$RUN_DIR/$f" ]] && { echo "----- $f"; cat "$RUN_DIR/$f"; }
done
tar -czf "$RUN_DIR.tgz" -C "$(dirname "$RUN_DIR")" --exclude='src-*' --exclude='wheel-*' --exclude='site-*' "$(basename "$RUN_DIR")"
echo
echo "full results: $RUN_DIR.tgz"
echo "finished    : $(date -Is)"
echo "scratch     : $SCRATCH  (auto-purges after 5 days)"
