#!/bin/bash
# Submit the short-target timing pipeline as two chained SLURM steps:
#
#   st-build (engine, sha256)  --afterok-->  st-timing (every tier, one node, 64 threads)
#
# Logs go to ~/projects/pepmatch/slurm-logs/ and results to
# ~/projects/pepmatch/results/short-targets-<stamp>/, both outside the repo, so `git status`
# stays clean.
#
# Usage (from the repo root, on branch bench/short-targets-timing):
#     SMOKE=1 bash benchmarking/slurm/submit-st-pipeline.sh    # FIRST: minutes, tiny sizes
#     bash benchmarking/slurm/submit-st-pipeline.sh            # the real run
# Env overrides: EXTRA (more script flags, e.g. EXTRA="--tiers A,B").
set -euo pipefail

SLURM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$HOME/projects/pepmatch"
LOGS="$BASE/slurm-logs"
STAMP="$(date +%Y%m%d-%H%M%S)"
if [[ "${SMOKE:-0}" == "1" ]]; then
  export RUN_DIR="$BASE/results/short-targets-smoke-$STAMP"
  export EXTRA="--smoke ${EXTRA:-}"
  TIME=(--time=02:00:00)
else
  export RUN_DIR="$BASE/results/short-targets-$STAMP"
  export EXTRA="${EXTRA:-}"
  TIME=()
fi
mkdir -p "$LOGS" "$RUN_DIR"

build=$(sbatch --parsable --export=ALL --output="$LOGS/%x-%j.out" --error="$LOGS/%x-%j.err" \
        "$SLURM_DIR/run-st-build.sbatch")
timing=$(sbatch --parsable --export=ALL --dependency=afterok:"$build" ${TIME[@]+"${TIME[@]}"} \
         --output="$LOGS/%x-%j.out" --error="$LOGS/%x-%j.err" "$SLURM_DIR/run-st-timing.sbatch")

[[ "${SMOKE:-0}" == "1" ]] && MODE=" (SMOKE)" || MODE=""
echo "short-target timing${MODE}"
echo "  build   = $build"
echo "  timing  = $timing  (held on $build)"
echo "  results = $RUN_DIR"
echo
echo "Queue:  squeue -u \$USER -o '%.12i %.12j %.8T %.10M %.20E'"
echo "Logs:   $LOGS/st-build-$build.out   $LOGS/st-timing-$timing.out"
