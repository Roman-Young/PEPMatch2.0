#!/bin/bash
#SBATCH --job-name=pepmatch-indel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# 1-indel benchmark: PEPMatch2.0 vs Brute Force, BLAST, DIAMOND, MMseqs2.
# Produces the two TIMING tables (COSMIC + CEDAR). Memory is deliberately NOT measured
# here -- see the -m note at the bottom.
#
#   sbatch benchmarking/slurm/run-indel-benchmark.sbatch
#
# WHY THE DATASETS RUN SEQUENTIALLY, NOT AS TWO JOBS:
# cosmic_indel and cedar_indel use the SAME proteome (human.fasta), and BLAST/DIAMOND/
# MMseqs2 each name their database after the proteome path. Two concurrent jobs would
# build and delete the identical database files underneath each other -- which does not
# reliably crash, it can produce a plausible but WRONG recall number. Sequential is the
# only safe arrangement without per-job database prefixes.
#
# --mem=64G is deliberately generous: DIAMOND with `-k 0 --ultra-sensitive -e 10000`
# (report every alignment, effectively no e-value filter) exhausted an 8GB machine
# during testing. Lower this only after watching a successful run's peak usage.

set -euo pipefail

REPO="$HOME/projects/pepmatch/PEPMatch2.0"
SCRATCH="/mnt/hpcscratch/ryoung/pepmatch-indel-${SLURM_JOB_ID:-manual}"
# Generated tables land in benchmarking/results/ (see benchmarking.py:results_dir).
RESULTS="$REPO/benchmarking/results"

echo "=== PEPMatch indel benchmark ==="
echo "job        : ${SLURM_JOB_ID:-manual}"
echo "host       : $(hostname)"
echo "started    : $(date -Is)"
echo "repo       : $REPO"
echo "scratch    : $SCRATCH"

# --- environment -------------------------------------------------------------------
# A batch shell does not source the same rc files as an interactive login, so the
# module/venv/PATH setup must be explicit or every aligner silently resolves to nothing.
source "$HOME/projects/pepmatch/env.sh"

THREADS="${SLURM_CPUS_PER_TASK:-1}"
echo "threads    : $THREADS"

# --- fail fast if a tool is missing -------------------------------------------------
# Without this, a missing binary surfaces much later as a confusing pandas error.
for tool in blastp makeblastdb diamond mmseqs python3; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "FATAL: '$tool' not on PATH after sourcing env.sh -- fix before submitting." >&2
    exit 1
  fi
  echo "  found $tool -> $(command -v "$tool")"
done

python3 -c 'import pandas, pyarrow, polars' || {
  echo "FATAL: missing python deps (need pandas, pyarrow, polars): pip install pyarrow" >&2
  exit 1
}

python3 -c 'import pepmatch; from pepmatch import _rs; print("  ENGINE OK", pepmatch.__file__)' || {
  echo "FATAL: PEPMatch Rust engine not importable -- run 'maturin develop --release' first." >&2
  exit 1
}

# --- scratch layout -----------------------------------------------------------------
# Tool databases are written NEXT TO the proteome, and the MMseqs2 index alone is
# multi-GB, so the proteome is copied to scratch and every generated database follows
# it there. Only the small results tables are written back to home.
mkdir -p "$SCRATCH/proteomes" "$SCRATCH/work" "$RESULTS/logs"
cp "$REPO/benchmarking/proteomes/human.fasta" "$SCRATCH/proteomes/"

export PEPMATCH_BENCH_PROTEOME_DIR="$SCRATCH/proteomes"
export PEPMATCH_BENCH_LOGDIR="$RESULTS/logs"

# PEPMatch's index (~276MB per k) and the aligners' temp files land in the CWD.
cd "$SCRATCH/work"

# --- pre-flight: coordinate convention ----------------------------------------------
# If a tool reports 0-based positions where the expected file is 1-based, every row
# misses and that tool scores a FALSE 0%. Verified for BLAST/PEPMatch/Brute Force by
# code reading; DIAMOND and MMseqs2 pass their raw sstart/tstart through unadjusted, so
# this gate checks them against a peptide at a known position before we burn hours.
echo
echo "=== pre-flight: coordinate convention ==="
if ! python3 "$REPO/benchmarking/scripts/verify_coordinates.py"; then
  echo "FATAL: coordinate check failed -- recall numbers would be meaningless." >&2
  exit 1
fi

# --- the benchmarks -----------------------------------------------------------------
# Results are written incrementally (one row per method), so even a wall-clock kill
# leaves a partial table rather than nothing.
for dataset in cosmic_indel cedar_indel; do
  echo
  echo "############################################################"
  echo "  $dataset  ($(date -Is))"
  echo "############################################################"
  python3 "$REPO/benchmarking/benchmarking.py" -b "$dataset" -p "$THREADS"
  echo "  -> $RESULTS/${dataset}_benchmarking.tsv"
done

# --- summary ------------------------------------------------------------------------
echo
echo "=== RESULTS ==="
for dataset in cosmic_indel cedar_indel; do
  f="$RESULTS/${dataset}_benchmarking.tsv"
  echo
  echo "--- $dataset ---"
  if [[ -f "$f" ]]; then
    column -t -s $'\t' "$f" 2>/dev/null || cat "$f"
    # Surface failures loudly: a FAILED row is easy to miss in a long log.
    if grep -q "FAILED\|SKIPPED" "$f"; then
      echo "  !! this table contains FAILED/SKIPPED methods -- see $PEPMATCH_BENCH_LOGDIR"
    fi
  else
    echo "  MISSING -- the run did not reach this dataset"
  fi
done

echo
echo "finished   : $(date -Is)"
echo "scratch    : $SCRATCH  (auto-purges after 5 days; results are in $RESULTS)"
echo
echo "NOTE: memory was not measured in this run. Add -m to the benchmarking.py calls"
echo "      to collect it -- but it re-runs every method, roughly doubling wall time,"
echo "      and the current numbers are not comparable across in-process vs subprocess"
echo "      tools (the Python interpreter baseline dominates). Treat as a separate job."
