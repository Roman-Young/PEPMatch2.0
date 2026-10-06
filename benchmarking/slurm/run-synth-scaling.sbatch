#!/bin/bash
#SBATCH --job-name=synth-scaling
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
# 5 days. The naive baseline is O(proteome) per query, so wall time is set by the largest
# size: 1-indel 100k ~= 43h, 2-indel 100k ~= 74h, plus ~2x observed node-to-node variance.
# 5 days covers 1-indel to 100k and 2-indel to 10k with real margin. For 2-indel 100k, push
# toward the 7-day `normal` QOS ceiling. A generous wall costs nothing -- the job ends when
# it finishes, not at the wall.
#SBATCH --time=5-00:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# THE SCALING FIGURE: PEPMatch 2.0 vs Brute Force across query counts.
#
#   INDELS=1 sbatch benchmarking/slurm/run-synth-scaling.sbatch
#   INDELS=2 sbatch benchmarking/slurm/run-synth-scaling.sbatch
#
# BASELINE CHANGED 2026-08-31. This sweep now races PEPMatch against `brute_force_naive`,
# which sweeps every query against every protein with no prefilter. The previous baseline
# (`brute_force`) narrowed candidates with a pigeonhole prefilter -- PEPMatch's OWN seeding
# strategy -- so the old figure compared Rust+index against Python+the-same-filter rather
# than against brute force. The two are proven to return identical rows
# (scripts/verify-naive-equivalence.py), so ONLY the Brute Force time column moves; recall
# is untouched and every published recall number stands.
#
# Consequence: the naive baseline is ~200x slower, so the old N ceilings are gone. See the
# runtime projection printed below before overriding SIZES.
#
# Three deliberate constraints, all of which change the numbers if broken:
#
# 1. PEPMatch + Brute Force ONLY (--methods). The aligners belong to the separate
#    fixed-size table (run-synth-aligners.sbatch, 5,000 queries). Running BLAST over a
#    1M-query set would cost hundreds of hours and terabytes of RSS to produce a column
#    the figure does not use.
#
# 2. The BASELINE MUST NOT SHARE PEPMatch's candidate generation. That is the whole point
#    of the rebuild: any prefilter that seeds on shared k-mers is PEPMatch's algorithm, and
#    a baseline running PEPMatch's algorithm cannot measure PEPMatch's advantage. Do not
#    "optimise" brute_force_naive by adding a seed/index step to make the job finish sooner
#    -- that silently reverts this fix. Reduce SIZES instead.
#
# 3. ONE job, ONE node, sizes run SEQUENTIALLY. A scaling curve only means something if
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
INDELS="${INDELS:-1}"
# The baseline. `brute_force_naive` (default, 2026-08-31) sweeps every query against every
# protein with NO prefilter -- the honest brute force. `brute_force` is the older
# pigeonhole-prefiltered reference; it shares PEPMatch's own seeding strategy, so its
# speedup UNDERSTATES the real one. Kept runnable for comparison, never for a published
# figure.
BF_METHOD="${BF_METHOD:-brute_force_naive}"
# Sizes. The naive baseline is O(proteome) per query with NO candidate narrowing, so its
# cost is ~flat in query length and the sweep's wall time is set almost entirely by the
# largest N. Both indel counts run the SAME generic oracle (one positional pass,
# early-exit), so cost rises with edit budget as it should. Measured per-query cost
# (dev box, at 32 threads):
#
#   1 indel : ~50 core-s/query  ->  ~1.6 s/query  ->  10k ~= 4h,   100k ~= 43h,  1M ~= 18 days
#   2 indels: ~85 core-s/query  ->  ~2.7 s/query  ->  10k ~= 7h,   100k ~= 74h,  1M ~= 30 days
#
# 100k is reachable for BOTH. 1-indel 100k (~2 days nominal) fits the 5-day wall with margin.
# 2-indel 100k is ~3.4 days nominal, but the ~2x node-to-node variance can push it near 7
# days, so run THAT one under the `long` QOS (14-day ceiling) so a slow node cannot wall-kill
# a multi-day job:
#     INDELS=2 SIZES="100 1k 10k 100k" sbatch --qos=long --time=10-00:00:00 <this script>
# 1M is out of reach for a true brute force at either edit count (~18-30 days) -- that
# infeasibility is itself a reportable result, not a gap to hide. Default is capped at 10k so
# an unattended run always finishes well inside the wall; pass SIZES=... to go higher.
if [[ -z "${SIZES:-}" ]]; then SIZES="100 1k 10k"; fi

echo "=== synthetic ${INDELS}-indel SCALING sweep (PEPMatch vs Brute Force) ==="
echo "job        : ${SLURM_JOB_ID:-manual}"
echo "host       : $(hostname)"
echo "started    : $(date -Is)"

source "$HOME/projects/pepmatch/env.sh"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
echo "threads    : $THREADS"
echo "sizes      : $SIZES"
echo "baseline   : $BF_METHOD"

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

# Cost calibration. The naive baseline is O(proteome) per query with no candidate
# narrowing, so the sweep's wall time is set almost entirely by the largest N. Measure the
# real per-query cost on THIS node against the REAL queries and print the projection
# BEFORE committing hours to it -- the alternative is discovering the sweep could never
# finish 30 hours in. Informational by design: it prints and continues rather than
# guessing a threshold and refusing a run Roman meant to start.
if [[ "$BF_METHOD" == "brute_force_naive" && "${CALIBRATE:-1}" == "1" ]]; then
  echo; echo "=== pre-flight: naive baseline cost calibration (${CAL_N:-24} queries) ==="
  largest="${SIZES##* }"
  python3 "$REPO/benchmarking/scripts/verify-naive-equivalence.py" \
    --indels "$INDELS" \
    --queries "$REPO/benchmarking/queries/synth-${INDELS}indel-${largest}.fasta" \
    --limit "${CAL_N:-24}" --threads "$THREADS" --skip-prefiltered \
    || echo "  (calibration failed -- continuing; it is informational only)"
fi

for size in $SIZES; do
  dataset="synth_${INDELS}indel_${size}"
  echo
  echo "############### SCALING: $dataset  ($(date -Is)) ###############"
  python3 "$REPO/benchmarking/benchmarking.py" \
    -b "$dataset" \
    --methods "PEPMatch,$BF_METHOD" \
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
