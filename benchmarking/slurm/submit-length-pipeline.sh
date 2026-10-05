#!/bin/bash
# Submit the per-query-length timing pipeline as three chained SLURM steps:
#
#   length-build (1 job)  --afterok-->  length-generate (array, 1 task per set)
#                         --afterok-->  length-timing (1 job, 1 node)
#
# afterok on an array waits for EVERY task to exit 0, and a task exits 0 only if its set
# self-certified (PEPMatch and Brute Force both 100.0 recall). So timing never starts on
# a bad or missing set, and nothing needs babysitting.
#
# Usage (from the repo root, on branch feature/per-length-timing):
#     SMOKE=1 bash benchmarking/slurm/submit-length-pipeline.sh    # FIRST: ~30 min end to end
#     bash benchmarking/slurm/submit-length-pipeline.sh            # the real run
#
# Env overrides: CONCURRENT (max array tasks at once, default 8).
set -euo pipefail

SLURM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONCURRENT="${CONCURRENT:-8}"

if [[ "${SMOKE:-0}" == "1" ]]; then
  # Tiny sets that exercise every code path: both indel counts, the 2-indel 8-mer k=2
  # partition (the headline point of the figure), a fixed-k/per-length split (L=12), and
  # brute force. Separate prefix + results dir, so the smoke output can never be mistaken
  # for the real run.
  TASKS="2:8 2:12 1:8 1:12"
  export PREFIX=synth-len-smoke NMAX=100 SIZES=50,100 AUDIT=20 CERT_N=50
  export LENGTHS=8,12 TIME_SIZES=50,100 REPEATS=1 BF_N=4
  export RESULTS="$HOME/projects/pepmatch/PEPMatch2.0/benchmarking/results/per-length-smoke"
  GEN_TIME=(--time=02:00:00)
  TIME_TIME=(--time=02:00:00)
else
  # Slowest sets first, so the 2-indel short queries (hours) start while the rest
  # (minutes) fill in around them.
  TASKS=""
  for d in 2 1; do for L in $(seq 8 25); do TASKS+="$d:$L "; done; done
  TASKS="${TASKS% }"
  GEN_TIME=()
  TIME_TIME=()
fi
export TASKS
N=$(wc -w <<< "$TASKS")

build=$(sbatch --parsable "$SLURM_DIR/run-length-build.sbatch")
gen=$(sbatch --parsable --dependency=afterok:"$build" --array=0-$((N - 1))%"$CONCURRENT" \
      ${GEN_TIME[@]+"${GEN_TIME[@]}"} "$SLURM_DIR/run-length-generate.sbatch")
timing=$(sbatch --parsable --dependency=afterok:"$gen" ${TIME_TIME[@]+"${TIME_TIME[@]}"} \
         "$SLURM_DIR/run-length-timing.sbatch")

[[ "${SMOKE:-0}" == "1" ]] && MODE=" (SMOKE)" || MODE=""
echo "per-length pipeline${MODE}: $N sets"
echo "  build    = $build"
echo "  generate = $gen  (array 0-$((N - 1)), max $CONCURRENT at once; held on $build)"
echo "  timing   = $timing  (held on every generate task)"
echo
echo "Queue:  squeue -u \$USER -o '%.14i %.9P %.18j %.8T %.10M %.20E'"
echo "Logs:   length-build-$build.out, length-generate-${gen}_<task>.out, length-timing-$timing.out"
