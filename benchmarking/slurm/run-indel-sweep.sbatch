#!/bin/bash
#SBATCH --job-name=pepmatch-sweep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# max_target_seqs SWEEP: one full benchmark (all methods, both datasets) at a single
# -max_target_seqs / --max-seqs value, so BLAST and MMseqs2 can be compared against the
# rest of the tools at that cap. Submit once per value -- they run on separate nodes
# concurrently and cannot collide, because each job gets its OWN scratch (job id) for
# databases and its OWN results/sweep-<value>/ folder for the tables.
#
#   sbatch --export=ALL,MAXSEQS=100  benchmarking/slurm/run-indel-sweep.sbatch
#   sbatch --export=ALL,MAXSEQS=1000 benchmarking/slurm/run-indel-sweep.sbatch
#
# Each writes to benchmarking/results/sweep-<MAXSEQS>/. Compare the BLAST/MMseqs2 rows
# across the sweep folders (and against the existing 100000 run) to pick the value that
# best trades recall for runtime; that value's full table is then the final timing table.
#
# Only BLAST and MMseqs2 respond to MAXSEQS; PEPMatch/DIAMOND/Brute Force are unaffected,
# but the full suite is run so every table is complete and directly comparable.

set -euo pipefail

: "${MAXSEQS:?set MAXSEQS, e.g. sbatch --export=ALL,MAXSEQS=100 ...}"

REPO="$HOME/projects/pepmatch/PEPMatch2.0"
SCRATCH="/mnt/hpcscratch/ryoung/pepmatch-sweep-${MAXSEQS}-${SLURM_JOB_ID:-manual}"
SUBDIR="sweep-${MAXSEQS}"
RESULTS="$REPO/benchmarking/results/$SUBDIR"

echo "=== PEPMatch indel sweep: max_target_seqs=$MAXSEQS ==="
echo "job        : ${SLURM_JOB_ID:-manual}"
echo "host       : $(hostname)"
echo "started    : $(date -Is)"
echo "results    : $RESULTS"

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

echo
echo "=== pre-flight: coordinate convention ==="
if ! python3 "$REPO/benchmarking/scripts/verify_coordinates.py"; then
  echo "FATAL: coordinate check failed -- recall numbers would be meaningless." >&2
  exit 1
fi

for dataset in cosmic_indel cedar_indel; do
  echo
  echo "############################################################"
  echo "  $dataset  max_target_seqs=$MAXSEQS  ($(date -Is))"
  echo "############################################################"
  python3 "$REPO/benchmarking/benchmarking.py" -b "$dataset" -p "$THREADS" \
    --max-target-seqs "$MAXSEQS" --results-subdir "$SUBDIR"
done

echo
echo "=== RESULTS (max_target_seqs=$MAXSEQS) ==="
for dataset in cosmic_indel cedar_indel; do
  f="$RESULTS/${dataset}_benchmarking.tsv"
  echo; echo "--- $dataset ---"
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
