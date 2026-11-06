#!/bin/bash
#SBATCH --job-name=st-timing
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#
# Step 2: every timing tier in ONE job on ONE node, sequentially, 64 threads -- the same
# rules as Figures 4/5. Identical runs vary ~2x between nodes, so splitting tiers across jobs
# would put node variance into the comparisons. All logic is in
# benchmarking/scripts/time-short-targets.py (see its docstring). Generous wall: the job
# ends when it finishes, and each config is capped (--cap-min) inside the script.

set -euo pipefail
REPO="$HOME/projects/pepmatch/PEPMatch2.0"
: "${RUN_DIR:?submit through submit-st-pipeline.sh, which sets RUN_DIR}"
EXTRA="${EXTRA:-}"

echo "=== short-target timing ==="
echo "job     : ${SLURM_JOB_ID:-manual}"
echo "host    : $(hostname)"
echo "run dir : $RUN_DIR"
echo "started : $(date -Is)"
source "$HOME/projects/pepmatch/env.sh"
python3 -c 'import polars' || { echo "FATAL: polars missing from the venv." >&2; exit 1; }

SCRATCH="/mnt/hpcscratch/ryoung/short-targets-${SLURM_JOB_ID:-manual}"
mkdir -p "$SCRATCH"
cp "$REPO/benchmarking/proteomes/human.fasta" "$SCRATCH/"

python3 "$REPO/benchmarking/scripts/time-short-targets.py" \
  --workdir "$SCRATCH" --outdir "$RUN_DIR" --proteome "$SCRATCH/human.fasta" \
  --threads "${SLURM_CPUS_PER_TASK:-1}" --expect-so-sha256 "$(cat "$RUN_DIR/engine-sha256.txt")" $EXTRA

echo
echo "=== paste these back to Kairo ==="
for f in run-info.txt short-target-summary.tsv short-target-oracle.tsv; do
  [[ -f "$RUN_DIR/$f" ]] && { echo "----- $f"; cat "$RUN_DIR/$f"; }
done
tar -czf "$RUN_DIR.tgz" -C "$(dirname "$RUN_DIR")" "$(basename "$RUN_DIR")"
echo
echo "full results: $RUN_DIR.tgz"
echo "finished    : $(date -Is)"
echo "scratch     : $SCRATCH  (auto-purges after 5 days)"
