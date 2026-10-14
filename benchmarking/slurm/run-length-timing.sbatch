#!/bin/bash
#SBATCH --job-name=length-timing
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
# 12h. PEPMatch timing (11 index builds + 2 arms x 36 sets x 2 sizes x 3 repeats x 2
# timers) projects well under 2h; the 64-query brute-force contrast is ~1h. Generous wall,
# no cost: the job ends when it finishes.
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# Step 3 of the per-query-length timing pipeline: the numbers behind Figure 5.
#
# ONE job, ONE node, everything sequential, 64 threads -- the same rules as Figure 4.
# Identical runs vary ~2x between nodes (METHODS sec 4), so splitting lengths across jobs
# or an array would put that variance into the curve and it would read as a length
# effect. Do not convert this to an array.
#
# All logic lives in scripts/time-per-length.py (see its docstring): per-length and
# fixed-k arms, engine-only and harness timers, marginal per-query cost, recall on every
# set, naive brute force on 64 queries per length.

set -euo pipefail

REPO="$HOME/projects/pepmatch/PEPMatch2.0"
OUT="$REPO/benchmarking/results/per-length"
PREFIX="${PREFIX:-synth-len}"
LENGTHS="${LENGTHS:-8-25}"
INDELS="${INDELS:-1,2}"
TIME_SIZES="${TIME_SIZES:-1000,10000}"
REPEATS="${REPEATS:-3}"
BF_N="${BF_N:-64}"
RESULTS="${RESULTS:-$OUT}"

echo "=== per-length timing ==="
echo "job        : ${SLURM_JOB_ID:-manual}"
echo "host       : $(hostname)"
echo "started    : $(date -Is)"

source "$HOME/projects/pepmatch/env.sh"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
export PEPMATCH_SEED_STRATEGY=rarest

python3 -c 'import pandas, polars' || { echo "FATAL: missing pandas/polars." >&2; exit 1; }

SCRATCH="/mnt/hpcscratch/ryoung/pepmatch-length-timing-${SLURM_JOB_ID:-manual}"
mkdir -p "$SCRATCH/proteomes" "$SCRATCH/work"
cp "$REPO/benchmarking/proteomes/human.fasta" "$SCRATCH/proteomes/"
export PEPMATCH_BENCH_PROTEOME_DIR="$SCRATCH/proteomes"
cd "$SCRATCH/work"

python3 "$REPO/benchmarking/scripts/time-per-length.py" \
  --prefix "$PREFIX" --indels "$INDELS" --lengths "$LENGTHS" --sizes "$TIME_SIZES" \
  --repeats "$REPEATS" --bf-n "$BF_N" --threads "$THREADS" \
  --workdir "$SCRATCH/work" --outdir "$RESULTS" \
  --expect-so-sha256 "$(cat "$OUT/engine-sha256.txt")"

echo
echo "=== paste these back to Kairo ==="
for f in run-info.txt index-build.tsv per-length-summary.tsv per-length-bf.tsv; do
  [[ -f "$RESULTS/$f" ]] && { echo "----- $f"; cat "$RESULTS/$f"; }
done
echo
echo "finished   : $(date -Is)"
echo "scratch    : $SCRATCH  (auto-purges after 5 days)"
