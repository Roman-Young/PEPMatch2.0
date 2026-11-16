#!/bin/bash
# Submit the per-length-k comparison as two chained SLURM steps:
#
#   plk-build (base + head, one copy of head per cap)  --afterok-->  plk-timing (one node)
#
# base = one k per batch (master + PRs #39-#42); capN = every query at its own k, at most N.
# Logs go to ~/projects/pepmatch/slurm-logs/ and results to
# ~/projects/pepmatch/results/per-length-k-<stamp>/, both outside the repo.
#
# Usage (from the repo root, on branch bench/per-length-k, after `git fetch origin`):
#     SMOKE=1 bash benchmarking/slurm/submit-plk-pipeline.sh    # FIRST: minutes, tiny sizes
#     bash benchmarking/slurm/submit-plk-pipeline.sh            # the real run
# Env overrides: BASE_COMMIT, HEAD_COMMIT, CAPS (default "4 5 6 7 8 10 16"), THREADS (default "64,8"),
#                EXTRA (more script flags, e.g. EXTRA="--sets mixed,coronavirus --modes mm2").
set -euo pipefail

SLURM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$HOME/projects/pepmatch"
LOGS="$BASE/slurm-logs"
STAMP="$(date +%Y%m%d-%H%M%S)"
export BASE_COMMIT="${BASE_COMMIT:-172207f}"
export HEAD_COMMIT="${HEAD_COMMIT:-95f857b}"
export CAPS="${CAPS:-4 5 6 7 8 10 16}"
if [[ "${SMOKE:-0}" == "1" ]]; then
  export RUN_DIR="$BASE/results/per-length-k-smoke-$STAMP"
  export EXTRA="--smoke ${EXTRA:-}"
  TIME=(--time=02:00:00)
else
  export RUN_DIR="$BASE/results/per-length-k-$STAMP"
  export EXTRA="${EXTRA:-}"
  TIME=()
fi
mkdir -p "$LOGS" "$RUN_DIR"

build=$(sbatch --parsable --export=ALL --output="$LOGS/%x-%j.out" --error="$LOGS/%x-%j.err" \
        "$SLURM_DIR/run-plk-build.sbatch")
timing=$(sbatch --parsable --export=ALL --dependency=afterok:"$build" ${TIME[@]+"${TIME[@]}"} \
         --output="$LOGS/%x-%j.out" --error="$LOGS/%x-%j.err" "$SLURM_DIR/run-plk-timing.sbatch")

[[ "${SMOKE:-0}" == "1" ]] && MODE=" (SMOKE)" || MODE=""
echo "per-length k ceilings${MODE}: base $BASE_COMMIT vs caps [$CAPS] of $HEAD_COMMIT"
echo "  build   = $build"
echo "  timing  = $timing  (held on $build)"
echo "  results = $RUN_DIR"
echo
echo "Queue:  squeue -u \$USER -o '%.12i %.12j %.8T %.10M %.20E'"
echo "Logs:   $LOGS/plk-build-$build.out   $LOGS/plk-timing-$timing.out"
