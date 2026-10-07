#!/bin/bash
#SBATCH --job-name=figure3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
# 12h. The aligners dominate (BLAST is minutes per dataset); the memory pass roughly
# doubles the timing pass. A generous wall is free: the job ends when it finishes.
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# FIGURE 3 RERUN with the annotation speedup: PEPMatch 2.0 vs BLAST (short + default),
# DIAMOND and MMseqs2 on COSMIC, CEDAR and the 5,000-query 2-indel synthetic set.
#
#   git checkout bench/annotation-speedup
#   sbatch benchmarking/slurm/run-figure3.sbatch
#
# Same settings as the original Figure 3 runs (run-indel-final.sbatch and
# run-synth-aligners.sbatch): 32 threads, the config's max_target_seqs=100, timing pass
# then a separate memory pass. Three deliberate differences:
#   1. All three panels in ONE job on ONE node (the originals were two jobs, two nodes),
#      so every tool in the figure was measured on the same hardware.
#   2. Results go to results/figure3-<job>/, never results/ itself, which still holds the
#      original Figure 3 tables. Nothing from the earlier runs is overwritten.
#   3. Figure 3's tools only (--methods). Without the filter the harness would also run
#      the naive brute force, ~8 h of work for a column Figure 3 does not show.
#
# PROVENANCE GATE: refuses to run unless matcher.py and the installed engine are the
# exact files Phase 1 validated (job 23130188, commit a17baf1), by sha256.

set -euo pipefail

REPO="$HOME/projects/pepmatch/PEPMatch2.0"
EXPECT_MATCHER_SHA="${EXPECT_MATCHER_SHA:-6fb6bae2b1f42c3e404b9958b23282ed162cc39fd4789ae525f5b63d3ad5f63d}"
EXPECT_ENGINE_SHA="${EXPECT_ENGINE_SHA:-ee446c574d74c43e5ad560423ca114655af78a8b1bf9d8cbacd061940b74fa5b}"
METHODS="${METHODS:-PEPMatch,blast (short),blast (default),diamond,mmseqs2}"
DATASETS="${DATASETS:-cosmic_indel cedar_indel synth_2indel_5k}"
SUBDIR="figure3-${SLURM_JOB_ID:-manual}"
RESULTS="$REPO/benchmarking/results/$SUBDIR"

echo "=== Figure 3 rerun (annotation speedup) ==="
echo "job        : ${SLURM_JOB_ID:-manual}"
echo "host       : $(hostname)"
echo "started    : $(date -Is)"

cd "$REPO"
source "$HOME/projects/pepmatch/env.sh"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
echo "threads    : $THREADS"
echo "methods    : $METHODS"
echo "commit     : $(git rev-parse --short HEAD) ($(git rev-parse --abbrev-ref HEAD))"

echo; echo "=== provenance: same code Phase 1 validated ==="
SO="$(python3 -c 'from pepmatch import _rs; print(_rs.__file__)')"
PKG="$(python3 -c 'import pepmatch, os; print(os.path.dirname(pepmatch.__file__))')"
got_matcher="$(sha256sum "$PKG/matcher.py" | cut -c1-64)"
got_engine="$(sha256sum "$SO" | cut -c1-64)"
echo "matcher.py : $got_matcher"
echo "engine     : $got_engine"
if [[ "$got_matcher" != "$EXPECT_MATCHER_SHA" || "$got_engine" != "$EXPECT_ENGINE_SHA" ]]; then
  echo "FATAL: installed PEPMatch is not the build Phase 1 validated." >&2
  echo "       Check out bench/annotation-speedup and rerun the Phase 1 job (it rebuilds)," >&2
  echo "       then pass its fingerprint via EXPECT_MATCHER_SHA / EXPECT_ENGINE_SHA." >&2
  exit 1
fi
echo "  OK"

for tool in blastp makeblastdb diamond mmseqs python3; do
  command -v "$tool" >/dev/null 2>&1 || {
    echo "FATAL: '$tool' not on PATH after sourcing env.sh." >&2; exit 1; }
done
python3 -c 'import pandas, pyarrow, polars' || {
  echo "FATAL: missing python deps (pandas, pyarrow, polars)." >&2; exit 1; }

# Aligner databases are written next to the proteome: keep them on scratch, and give this
# job its own copy so no concurrent job can rebuild a database underneath it.
SCRATCH="/mnt/hpcscratch/ryoung/pepmatch-figure3-${SLURM_JOB_ID:-manual}"
mkdir -p "$SCRATCH/proteomes" "$SCRATCH/work" "$RESULTS/logs"
cp "$REPO/benchmarking/proteomes/human.fasta" "$SCRATCH/proteomes/"
export PEPMATCH_BENCH_PROTEOME_DIR="$SCRATCH/proteomes"
export PEPMATCH_BENCH_LOGDIR="$RESULTS/logs"
cd "$SCRATCH/work"

echo; echo "=== pre-flight: coordinate convention ==="
python3 "$REPO/benchmarking/scripts/verify_coordinates.py" || {
  echo "FATAL: coordinate check failed." >&2; exit 1; }

for dataset in $DATASETS; do
  echo; echo "############### TIMING: $dataset  ($(date -Is)) ###############"
  python3 "$REPO/benchmarking/benchmarking.py" -b "$dataset" -p "$THREADS" \
    --methods "$METHODS" --results-subdir "$SUBDIR"
done
for dataset in $DATASETS; do
  echo; echo "############### MEMORY: $dataset  ($(date -Is)) ###############"
  python3 "$REPO/benchmarking/benchmarking.py" -b "$dataset" -p "$THREADS" --memory-only \
    --methods "$METHODS" --results-subdir "$SUBDIR"
done

echo; echo "=== paste these back to Kairo ==="
for dataset in $DATASETS; do
  for kind in benchmarking memory; do
    f="$RESULTS/${dataset}_${kind}.tsv"
    echo "----- ${dataset}_${kind}.tsv"
    if [[ -f "$f" ]]; then cat "$f"; else echo "MISSING"; fi
  done
done
grep -l "FAILED\|SKIPPED" "$RESULTS"/*_benchmarking.tsv 2>/dev/null \
  && echo "!! a method FAILED or was SKIPPED -- see $RESULTS/logs" || true

echo
echo "finished   : $(date -Is)"
echo "results    : $RESULTS"
