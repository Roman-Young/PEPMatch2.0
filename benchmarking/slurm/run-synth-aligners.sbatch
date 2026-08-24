#!/bin/bash
#SBATCH --job-name=synth-aligners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# THE ALIGNER TABLE: every method (PEPMatch, Brute Force, BLAST short + default,
# DIAMOND, MMseqs2) on the synthetic set at ONE fixed size -- 5,000 queries -- in the
# same speed/recall/memory shape as the existing COSMIC/CEDAR tables.
#
#   sbatch benchmarking/slurm/run-synth-aligners.sbatch
#
# Why 5,000 and not the scaling ladder: the aligners are being characterised, not scaled.
# BLAST's peak RSS grows with query count (18.4 GB was measured on 3,510 queries at
# max_target_seqs=100000), so pushing it up the ladder buys nothing the figure uses.
# max_target_seqs=100 is the config default since the sweep showed recall plateaus there
# while search time drops ~46x.
#
# Datasets sharing human.fasta must run SEQUENTIALLY: BLAST/DIAMOND/MMseqs2 each name
# their database after the proteome path, so two concurrent jobs corrupt each other's DB
# and produce plausible-but-wrong recall. This job is the only one that touches aligners.
#
# Writes to results/ root (no subdir) so scripts/build-results-doc.py finds the tables.

set -euo pipefail

REPO="$HOME/projects/pepmatch/PEPMatch2.0"
SCRATCH="/mnt/hpcscratch/ryoung/pepmatch-synth-aligners-${SLURM_JOB_ID:-manual}"
RESULTS="$REPO/benchmarking/results"
SIZE="${SIZE:-5k}"
INDELS="${INDELS:-1}"
DATASET="synth_${INDELS}indel_${SIZE}"

echo "=== synthetic ${INDELS}-indel ALIGNER table ($DATASET) ==="
echo "job        : ${SLURM_JOB_ID:-manual}"
echo "host       : $(hostname)"
echo "started    : $(date -Is)"

source "$HOME/projects/pepmatch/env.sh"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
echo "threads    : $THREADS"

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

mkdir -p "$SCRATCH/proteomes" "$SCRATCH/work" "$RESULTS/logs"
cp "$REPO/benchmarking/proteomes/human.fasta" "$SCRATCH/proteomes/"
export PEPMATCH_BENCH_PROTEOME_DIR="$SCRATCH/proteomes"
export PEPMATCH_BENCH_LOGDIR="$RESULTS/logs"
cd "$SCRATCH/work"

echo; echo "=== pre-flight: coordinate convention ==="
python3 "$REPO/benchmarking/scripts/verify_coordinates.py" || {
  echo "FATAL: a tool reports 0-based coordinates -- recall would be a false 0%." >&2
  exit 1; }

echo; echo "=== pre-flight: dataset integrity ==="
python3 "$REPO/benchmarking/scripts/verify-synthetic-dataset.py" \
  --prefix "synth-${INDELS}indel" --fast || {
  echo "FATAL: synthetic dataset failed verification." >&2; exit 1; }

echo; echo "############### TIMING: $DATASET  ($(date -Is)) ###############"
python3 "$REPO/benchmarking/benchmarking.py" -b "$DATASET" -p "$THREADS"

echo; echo "############### MEMORY: $DATASET  ($(date -Is)) ###############"
python3 "$REPO/benchmarking/benchmarking.py" -b "$DATASET" --memory-only -p "$THREADS"

echo; echo "=== ALIGNER TABLE RESULTS ==="
for kind in benchmarking memory; do
  f="$RESULTS/${DATASET}_${kind}.tsv"
  echo; echo "--- ${DATASET} ${kind} ---"
  if [[ -f "$f" ]]; then
    column -t -s $'\t' "$f" 2>/dev/null || cat "$f"
    grep -q "FAILED\|SKIPPED" "$f" && echo "  !! FAILED/SKIPPED -- see $PEPMATCH_BENCH_LOGDIR"
  else
    echo "  MISSING -- the run did not reach this table"
  fi
done

echo
echo "finished   : $(date -Is)"
echo "Next: python3 $REPO/benchmarking/scripts/build-results-doc.py"
echo "scratch    : $SCRATCH  (auto-purges after 5 days)"
