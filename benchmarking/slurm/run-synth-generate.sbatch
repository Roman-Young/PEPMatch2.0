#!/bin/bash
#SBATCH --job-name=synth-generate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# GENERATE the synthetic 1-indel dataset (and its nested subsets) on the cluster.
#
# This is a cluster job rather than a dev-box script because ground truth is
# proteome-wide: every query is resolved against all 81,837 proteins with the committed
# oracle. At 1M queries that is ~200 core-hours -- about 6.5 h on 32 cores, but ~55 h on
# the 4-core dev box, which also has 8 GB and no swap.
#
#   sbatch benchmarking/slurm/run-synth-generate.sbatch
#   NMAX=100000 sbatch benchmarking/slurm/run-synth-generate.sbatch    # smaller run
#
# Writes into the repo (queries/ and expected/). The 100k/1m artifacts are gitignored;
# the small ones plus the manifest are meant to be committed.
#
# Cost scales with the LENGTH distribution, not just the count: a short query splits into
# short pigeonhole pieces, and a 2-mer seed hits ~70% of the proteome. Run the generator
# with --cost-estimate before changing --min-len or --indels.

set -euo pipefail

REPO="$HOME/projects/pepmatch/PEPMatch2.0"
NMAX="${NMAX:-1000000}"
INDELS="${INDELS:-1}"
SIZES="${SIZES:-100,1000,5000,10000,100000,1000000}"
AUDIT="${AUDIT:-200}"
# Length range. Unset => the generator's own default (8 for 1-indel, 12 for 2-indel).
# For a 2-indel set that matches the 1-indel length composition, pass MINLEN=8 -- but note
# 8-11mers are ~30x more expensive for 2-indel ground truth (a 2-mer seed hits ~70% of the
# proteome), so check --cost-estimate before a large 2-indel run.
MINLEN="${MINLEN:-}"
MAXLEN="${MAXLEN:-}"
LEN_ARGS=""
[[ -n "$MINLEN" ]] && LEN_ARGS="$LEN_ARGS --min-len $MINLEN"
[[ -n "$MAXLEN" ]] && LEN_ARGS="$LEN_ARGS --max-len $MAXLEN"

echo "=== synthetic ${INDELS}-indel dataset generation ==="
echo "job        : ${SLURM_JOB_ID:-manual}"
echo "host       : $(hostname)"
echo "started    : $(date -Is)"

source "$HOME/projects/pepmatch/env.sh"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
echo "threads    : $THREADS"
echo "n_max      : $NMAX"
echo "sizes      : $SIZES"

python3 -c 'import pandas' || { echo "FATAL: pandas missing." >&2; exit 1; }

# Project the cost before committing hours to it.
echo; echo "=== projected cost ==="
echo "len args   :${LEN_ARGS:- (generator default)}"
python3 "$REPO/benchmarking/scripts/build-synthetic-indel-benchmark.py" \
  --indels "$INDELS" --n-max "$NMAX" --sizes "$SIZES" $LEN_ARGS --threads "$THREADS" --cost-estimate

echo; echo "=== generating ==="
python3 "$REPO/benchmarking/scripts/build-synthetic-indel-benchmark.py" \
  --indels "$INDELS" \
  --n-max "$NMAX" \
  --sizes "$SIZES" \
  $LEN_ARGS \
  --threads "$THREADS" \
  --audit-sample "$AUDIT"

echo; echo "=== self-certifying the dataset (full gate, incl. recall) ==="
# The generation job proves its own output before any benchmark is launched: the fast
# gates (schema / coordinate literalism / nesting / 50-50 balance / manifest hashes) AND
# the real gate -- PEPMatch and Brute Force must each score EXACTLY 100.0 recall on the
# 100-query subset, scored through the harness's own recall(). A bad 1M dataset is caught
# here, in this job, instead of hours into the sweep.
#
# The recall gate builds PEPMatch's k=4 index, which Preprocessor writes to the CWD, so
# run it from a scratch work dir -- that keeps the 276 MB .pepidx out of the repo and off
# the shared home quota.
SCRATCH="/mnt/hpcscratch/ryoung/pepmatch-synth-gen-${SLURM_JOB_ID:-manual}"
mkdir -p "$SCRATCH/proteomes" "$SCRATCH/work"
cp "$REPO/benchmarking/proteomes/human.fasta" "$SCRATCH/proteomes/"
export PEPMATCH_BENCH_PROTEOME_DIR="$SCRATCH/proteomes"
(
  cd "$SCRATCH/work"
  python3 "$REPO/benchmarking/scripts/verify-synthetic-dataset.py" \
    --prefix "synth-${INDELS}indel" --run-tools-on 100 --threads "$THREADS"
) || { echo "FATAL: dataset failed self-certification -- do NOT run the benchmarks." >&2; exit 1; }

echo
echo "finished   : $(date -Is)"
echo "scratch    : $SCRATCH  (auto-purges after 5 days)"
echo "Next: sbatch benchmarking/slurm/run-synth-scaling.sbatch"
