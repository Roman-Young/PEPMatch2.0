#!/bin/bash
#SBATCH --job-name=synth-scaling
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# THE SCALING FIGURE: PEPMatch 2.0 vs Brute Force across 100 -> 1,000,000 queries.
#
#   sbatch benchmarking/slurm/run-synth-scaling.sbatch
#
# Two deliberate constraints, both of which change the numbers if broken:
#
# 1. PEPMatch + Brute Force ONLY (--methods). The aligners belong to the separate
#    fixed-size table (run-synth-aligners.sbatch, 5,000 queries). Running BLAST over a
#    1M-query set would cost hundreds of hours and terabytes of RSS to produce a column
#    the figure does not use.
#
# 2. ONE job, ONE node, sizes run SEQUENTIALLY. A scaling curve only means something if
#    every point was measured on the same hardware: the max_target_seqs sweep showed the
#    SAME tool and settings varying ~2x between nodes (Brute Force COSMIC 310s vs 633s).
#    Splitting the sizes across a SLURM array would scatter that variance through the
#    curve and read as algorithmic behaviour. Do not convert this to an array job.
#
# Results land in results/synth-scaling/ so they never mix with the aligner table.

set -euo pipefail

REPO="$HOME/projects/pepmatch/PEPMatch2.0"
SCRATCH="/mnt/hpcscratch/ryoung/pepmatch-synth-scaling-${SLURM_JOB_ID:-manual}"
RESULTS="$REPO/benchmarking/results"
SIZES="${SIZES:-100 1k 10k 100k 1m}"
INDELS="${INDELS:-1}"

echo "=== synthetic ${INDELS}-indel SCALING sweep (PEPMatch vs Brute Force) ==="
echo "job        : ${SLURM_JOB_ID:-manual}"
echo "host       : $(hostname)"
echo "started    : $(date -Is)"

source "$HOME/projects/pepmatch/env.sh"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
echo "threads    : $THREADS"
echo "sizes      : $SIZES"

# No aligner binaries are needed here -- deliberately not checked, so a missing blastp
# cannot block the figure.
python3 -c 'import pandas, pyarrow, polars' || {
  echo "FATAL: missing python deps (pandas, pyarrow, polars)." >&2; exit 1; }
python3 -c 'import pepmatch; from pepmatch import _rs; print("  ENGINE OK")' || {
  echo "FATAL: PEPMatch Rust engine not importable -- run 'maturin develop --release'." >&2
  exit 1; }

# --- scratch. The CWD matters: Preprocessor.preprocess() rewrites <proteome>_<k>mers.pepidx
# unconditionally (no existence check) and Matcher mmaps it. Two PEPMatch runs sharing a
# working directory can therefore read a half-written 276 MB index -- silently wrong
# numbers, not a crash. A per-job scratch CWD is what keeps that impossible.
mkdir -p "$SCRATCH/proteomes" "$SCRATCH/work" "$RESULTS/logs"
cp "$REPO/benchmarking/proteomes/human.fasta" "$SCRATCH/proteomes/"
export PEPMATCH_BENCH_PROTEOME_DIR="$SCRATCH/proteomes"
export PEPMATCH_BENCH_LOGDIR="$RESULTS/logs"
cd "$SCRATCH/work"

echo; echo "=== pre-flight: dataset integrity ==="
python3 "$REPO/benchmarking/scripts/verify-synthetic-dataset.py" \
  --prefix "synth-${INDELS}indel" --fast || {
  echo "FATAL: synthetic dataset failed verification -- refusing to publish numbers." >&2
  exit 1; }

# Confirm every size THIS job will iterate actually exists before starting. The verifier
# above only checks the sizes the manifest happens to describe, so a generation run with a
# smaller --n-max than $SIZES would sail past it and then die mid-sweep -- hours in, after
# the cheap points have already burned a node.
echo; echo "=== pre-flight: every requested size is present ==="
missing=0
for size in $SIZES; do
  q="$REPO/benchmarking/queries/synth-${INDELS}indel-${size}.fasta"
  e="$REPO/benchmarking/expected/synth-${INDELS}indel-${size}-expected.tsv"
  for path in "$q" "$e"; do
    if [[ ! -f "$path" ]]; then
      echo "  MISSING: $path" >&2
      missing=1
    fi
  done
  [[ -f "$q" && -f "$e" ]] && echo "  ok: $size"
done
if [[ "$missing" -ne 0 ]]; then
  echo "FATAL: sizes requested in SIZES='$SIZES' are not on disk. Generate them first" >&2
  echo "       (slurm/run-synth-generate.sbatch) or narrow SIZES." >&2
  exit 1
fi

for size in $SIZES; do
  dataset="synth_${INDELS}indel_${size}"
  echo
  echo "############### SCALING: $dataset  ($(date -Is)) ###############"
  python3 "$REPO/benchmarking/benchmarking.py" \
    -b "$dataset" \
    --methods PEPMatch,brute_force \
    -p "$THREADS" \
    --results-subdir synth-scaling
done

echo; echo "=== SCALING RESULTS ==="
for size in $SIZES; do
  f="$RESULTS/synth-scaling/synth_${INDELS}indel_${size}_benchmarking.tsv"
  echo; echo "--- ${size} ---"
  if [[ -f "$f" ]]; then
    column -t -s $'\t' "$f" 2>/dev/null || cat "$f"
    grep -q "FAILED\|SKIPPED" "$f" && echo "  !! FAILED/SKIPPED -- see $PEPMATCH_BENCH_LOGDIR"
  else
    echo "  MISSING -- the run did not reach this size"
  fi
done

echo
echo "finished   : $(date -Is)"
echo "Next: python3 $REPO/benchmarking/scripts/plot-indel-scaling.py"
echo "scratch    : $SCRATCH  (auto-purges after 5 days)"
