#!/bin/bash
# Submit the full synthetic-indel scaling pipeline as chained SLURM jobs, for BOTH indel
# counts:
#
#     generate (MINLEN=9, no 8-mers)  --afterok-->  scaling sweep (naive brute force)
#
# The sweep for each indel count is submitted immediately but HELD by a dependency: it runs
# only if that count's generation job self-certified (exit 0, which requires PEPMatch and
# Brute Force to both score exactly 100% recall on the fresh dataset). So a dataset that
# fails its recall gate never reaches a multi-day sweep, and nothing has to be babysat.
#
# Usage (from the repo root):
#     bash benchmarking/slurm/submit-synth-pipeline.sh
#     CORES=32 bash benchmarking/slurm/submit-synth-pipeline.sh          # smaller node
#     SIZES_1="100 1k 10k" SIZES_2="100 1k 10k" bash ...                  # cap lower
#
# Env overrides:
#   CORES    cpus-per-task for every job            (default 64)
#   NMAX     master set size / largest sweep point  (default 100000)
#   MINLEN   min query length; 9 drops the 8-mers   (default 9)
#   SIZES_1  1-indel sweep points (space labels)    (default "100 1k 10k 100k")
#   SIZES_2  2-indel sweep points (space labels)    (default "100 1k 10k 100k")
#
# Notes:
#   - Generation ground truth uses the fast pigeonhole path, so it is cheap: ~<1 h for
#     1-indel 100k, ~5 h for 2-indel 100k at 64 cores. The sweep is the long part.
#   - At CORES<64 the 2-indel 100k sweep can approach the 7-day `normal` QOS wall under
#     ~2x node variance, so this script moves THAT job to `--qos=long` automatically.
#     At 64 cores the in-script 5-day wall is enough and no QOS change is made.
set -euo pipefail

SLURM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CORES="${CORES:-64}"
NMAX="${NMAX:-100000}"
MINLEN="${MINLEN:-9}"
GEN_SIZES="${GEN_SIZES:-100,1000,10000,100000}"   # comma list, for the generator's --sizes
SIZES_1="${SIZES_1:-100 1k 10k 100k}"             # space labels, for the sweep loop
SIZES_2="${SIZES_2:-100 1k 10k 100k}"

submit_chain() {
  local indels="$1" sweep_sizes="$2" extra_sweep="${3:-}"
  local gen sweep
  gen=$(INDELS="$indels" NMAX="$NMAX" MINLEN="$MINLEN" SIZES="$GEN_SIZES" \
        sbatch --parsable --cpus-per-task="$CORES" \
        "$SLURM_DIR/run-synth-generate.sbatch")
  sweep=$(INDELS="$indels" SIZES="$sweep_sizes" \
          sbatch --parsable --cpus-per-task="$CORES" --dependency=afterok:"$gen" \
          $extra_sweep "$SLURM_DIR/run-synth-scaling.sbatch")
  echo "  ${indels}-indel:  generate=${gen}  ->  sweep=${sweep}  (held on afterok:${gen})"
}

# 2-indel gets a longer wall only when it needs one (see header note).
sweep2_extra=""
if [[ "$CORES" -lt 64 ]]; then
  sweep2_extra="--qos=long --time=10-00:00:00"
fi

echo "Synthetic-indel pipeline: MINLEN=${MINLEN} (no 8-mers), NMAX=${NMAX}, ${CORES} cores, naive BF baseline"
submit_chain 1 "$SIZES_1"
submit_chain 2 "$SIZES_2" "$sweep2_extra"
echo
echo "Queue:   squeue -u \$USER -o '%.10i %.9P %.22j %.8T %.10M %.20E'"
echo "The sweeps show State=PENDING Reason=(Dependency) until their generate job finishes."
