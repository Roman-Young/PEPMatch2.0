#!/bin/bash
#SBATCH --job-name=pepmatch-final
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# FINAL 1-indel benchmark: timing AND memory, all in ONE job on ONE node at 32 cores,
# with max_target_seqs=100 (the sweep winner, now the config default). Producing every
# table from a single node keeps the numbers internally consistent -- the sweep showed
# fixed-tool timing varies ~2x between nodes, so timing and memory should not be stitched
# together across jobs.
#
#   sbatch benchmarking/slurm/run-indel-final.sbatch
#
# Writes the deliverable tables straight into benchmarking/results/ (no subdir), so
# scripts/build-results-doc.py picks them up directly:
#   {cosmic,cedar}_indel_benchmarking.tsv   (timing)
#   {cosmic,cedar}_indel_memory.tsv         (memory, + _memory_detail.tsv)
# These OVERWRITE the earlier 16-core / max_target_seqs=100000 numbers, which they
# supersede.

set -euo pipefail

REPO="$HOME/projects/pepmatch/PEPMatch2.0"
SCRATCH="/mnt/hpcscratch/ryoung/pepmatch-final-${SLURM_JOB_ID:-manual}"
RESULTS="$REPO/benchmarking/results"

echo "=== PEPMatch indel benchmark -- FINAL (timing + memory, 32 cores, max_target_seqs=100) ==="
echo "job        : ${SLURM_JOB_ID:-manual}"
echo "host       : $(hostname)"
echo "started    : $(date -Is)"

source "$HOME/projects/pepmatch/env.sh"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
echo "threads    : $THREADS"

# --- fail fast if a tool or the engine is missing -----------------------------------
for tool in blastp makeblastdb diamond mmseqs python3; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "FATAL: '$tool' not on PATH after sourcing env.sh." >&2
    exit 1
  fi
done
python3 -c 'import pandas, pyarrow, polars' || {
  echo "FATAL: missing python deps (pandas, pyarrow, polars)." >&2; exit 1; }
python3 -c 'import pepmatch; from pepmatch import _rs; print("  ENGINE OK")' || {
  echo "FATAL: PEPMatch Rust engine not importable -- run 'maturin develop --release'." >&2
  exit 1; }

# --- scratch (databases + indexes off the home quota) -------------------------------
mkdir -p "$SCRATCH/proteomes" "$SCRATCH/work" "$RESULTS/logs"
cp "$REPO/benchmarking/proteomes/human.fasta" "$SCRATCH/proteomes/"
export PEPMATCH_BENCH_PROTEOME_DIR="$SCRATCH/proteomes"
export PEPMATCH_BENCH_LOGDIR="$RESULTS/logs"
cd "$SCRATCH/work"

# --- pre-flight: coordinate convention ----------------------------------------------
echo; echo "=== pre-flight: coordinate convention ==="
if ! python3 "$REPO/benchmarking/scripts/verify_coordinates.py"; then
  echo "FATAL: coordinate check failed -- recall numbers would be meaningless." >&2
  exit 1
fi

# --- TIMING (max_target_seqs=100 comes from the config default) ---------------------
for dataset in cosmic_indel cedar_indel; do
  echo
  echo "############### TIMING: $dataset  ($(date -Is)) ###############"
  python3 "$REPO/benchmarking/benchmarking.py" -b "$dataset" -p "$THREADS"
done

# --- MEMORY (separate pass; re-runs each method, so it comes after timing) -----------
for dataset in cosmic_indel cedar_indel; do
  echo
  echo "############### MEMORY: $dataset  ($(date -Is)) ###############"
  python3 "$REPO/benchmarking/benchmarking.py" -b "$dataset" --memory-only -p "$THREADS"
done

# --- summary ------------------------------------------------------------------------
echo; echo "=== FINAL RESULTS ==="
for dataset in cosmic_indel cedar_indel; do
  for kind in benchmarking memory; do
    f="$RESULTS/${dataset}_${kind}.tsv"
    echo; echo "--- ${dataset} ${kind} ---"
    if [[ -f "$f" ]]; then
      column -t -s $'\t' "$f" 2>/dev/null || cat "$f"
      grep -q "FAILED\|SKIPPED" "$f" && echo "  !! FAILED/SKIPPED present -- see $PEPMATCH_BENCH_LOGDIR"
    else
      echo "  MISSING -- the run did not reach this table"
    fi
  done
done

echo
echo "finished   : $(date -Is)"
echo "Next: python3 $REPO/benchmarking/scripts/build-results-doc.py  (regenerates the workbook)"
echo "scratch    : $SCRATCH  (auto-purges after 5 days; results are in $RESULTS)"
