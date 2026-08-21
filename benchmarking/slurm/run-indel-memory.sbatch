#!/bin/bash
#SBATCH --job-name=pepmatch-indel-mem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# MEMORY table for the 1-indel benchmark. Run this AFTER the timing job
# (run-indel-benchmark.sbatch) -- it is a separate job on purpose.
#
#   sbatch benchmarking/slurm/run-indel-memory.sbatch
#
# WHY SEPARATE: measuring memory means running every method a second time (peak RSS can
# only be observed by actually executing the work). Folding that into the timing job
# would roughly double its wall time and put the headline deliverable at risk of a
# wall-clock kill. Timing first, memory second.
#
# WHAT IS REPORTED: peak resident memory attributable to the method, EXCLUDING the
# shared Python interpreter baseline. PEPMatch and brute force run inside this Python
# process (interpreter + pandas + polars = a fixed several-hundred-MB floor); the
# aligners exec standalone C binaries carrying no Python at all. Comparing raw peaks
# made the two in-process methods collapse onto an identical number dominated by that
# floor. Every raw component is written to the table so the subtraction can be audited.
#
# PHASES: memory is split into preprocessing and search (mirroring the timing table),
# because preprocessing is a one-time amortised cost while search is paid per query.
# Each is measured in its own process -- peak RSS is a high-water mark that cannot be
# reset within a process -- plus a third end-to-end run for the authoritative total,
# since the true peak can exceed max(phases).
#
# RUNTIME: three phases means roughly twice the work of a single pass. On COSMIC that
# is dominated by BLAST (~96 min short + ~31 min default, each searched twice), so
# budget ~5h; the 12h limit is deliberate headroom.

set -euo pipefail

REPO="$HOME/projects/pepmatch/PEPMatch2.0"
SCRATCH="/mnt/hpcscratch/ryoung/pepmatch-mem-${SLURM_JOB_ID:-manual}"
# Generated tables land in benchmarking/results/ (see benchmarking.py:results_dir).
RESULTS="$REPO/benchmarking/results"

echo "=== PEPMatch indel benchmark -- MEMORY ==="
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

python3 -c 'import pepmatch; from pepmatch import _rs' || {
  echo "FATAL: PEPMatch Rust engine not importable -- run 'maturin develop --release'." >&2
  exit 1
}

mkdir -p "$SCRATCH/proteomes" "$SCRATCH/work" "$RESULTS/logs"
cp "$REPO/benchmarking/proteomes/human.fasta" "$SCRATCH/proteomes/"

export PEPMATCH_BENCH_PROTEOME_DIR="$SCRATCH/proteomes"
export PEPMATCH_BENCH_LOGDIR="$RESULTS/logs"
cd "$SCRATCH/work"

# Sequential, same reason as the timing job: both datasets share human.fasta and the
# aligners name their databases after it, so concurrent jobs would corrupt each other.
for dataset in cosmic_indel cedar_indel; do
  echo
  echo "############################################################"
  echo "  $dataset  ($(date -Is))"
  echo "############################################################"
  python3 "$REPO/benchmarking/benchmarking.py" -b "$dataset" --memory-only -p "$THREADS"
  echo "  -> $RESULTS/${dataset}_memory.tsv"
done

echo
echo "=== MEMORY RESULTS ==="
for dataset in cosmic_indel cedar_indel; do
  f="$RESULTS/${dataset}_memory.tsv"
  echo
  echo "--- $dataset ---"
  if [[ -f "$f" ]]; then
    column -t -s $'\t' "$f" 2>/dev/null || cat "$f"
    if grep -q "FAILED\|SKIPPED" "$f"; then
      echo "  !! contains FAILED/SKIPPED methods -- see $PEPMATCH_BENCH_LOGDIR"
    fi
  else
    echo "  MISSING -- the run did not reach this dataset"
  fi
done

echo
echo "finished   : $(date -Is)"
echo "scratch    : $SCRATCH  (auto-purges after 5 days; results are in $RESULTS)"
