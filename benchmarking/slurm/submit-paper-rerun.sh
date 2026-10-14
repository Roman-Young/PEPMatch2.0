#!/bin/bash
# Re-time PEPMatch for every paper figure whose PEPMatch numbers changed with the
# post-profiling speedups. Aligners, brute force and index builds are unchanged and are NOT
# rerun. Run from the repo on herman, on branch paper/runs-with-speedups:
#
#   bash benchmarking/slurm/submit-paper-rerun.sh
#
#   Figure 3 (aligner comparison)  PEPMatch timing + memory, 32 cores -> results/rerun-aligners/
#   Figure 4 (scaling, 1 and 2 indels) PEPMatch timing, 64 cores     -> results/rerun-synth-scaling/
#   Figure 5 (per-length)          engine rebuild, then PEPMatch timing on the existing
#                                  per-length sets, brute force skipped -> results/per-length-rerun/
set -euo pipefail
SLURM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SLURM_DIR/../.." && pwd)"

[[ "$(git -C "$REPO" rev-parse --abbrev-ref HEAD)" == "paper/runs-with-speedups" ]] || {
  echo "FATAL: check out paper/runs-with-speedups first." >&2; exit 1; }
git -C "$REPO" diff --quiet -- pepmatch/ || {
  echo "FATAL: uncommitted changes under pepmatch/ -- the timed engine must be a commit." >&2; exit 1; }

# One build of the engine everyone uses; every timing job waits for it.
build=$(sbatch --parsable "$SLURM_DIR/run-length-build.sbatch")

fig3=$(DATASETS="cosmic_indel cedar_indel synth_1indel_5k synth_2indel_5k" SUBDIR=rerun-aligners MEMORY=1 \
       sbatch --parsable --dependency=afterok:"$build" --cpus-per-task=32 "$SLURM_DIR/run-pepmatch-rerun.sbatch")
fig4_1=$(DATASETS="synth_1indel_100 synth_1indel_1k synth_1indel_10k synth_1indel_100k" SUBDIR=rerun-synth-scaling \
         sbatch --parsable --dependency=afterok:"$build" --cpus-per-task=64 "$SLURM_DIR/run-pepmatch-rerun.sbatch")
fig4_2=$(DATASETS="synth_2indel_100 synth_2indel_1k synth_2indel_10k synth_2indel_100k" SUBDIR=rerun-synth-scaling \
         sbatch --parsable --dependency=afterok:"$build" --cpus-per-task=64 "$SLURM_DIR/run-pepmatch-rerun.sbatch")
fig5=$(BF_N=0 RESULTS="$REPO/benchmarking/results/per-length-rerun" \
       sbatch --parsable --dependency=afterok:"$build" "$SLURM_DIR/run-length-timing.sbatch")

echo "engine build      = $build"
echo "Figure 3          = $fig3   (held on $build)"
echo "Figure 4, 1 indel = $fig4_1"
echo "Figure 4, 2 indel = $fig4_2"
echo "Figure 5          = $fig5"
echo
echo "Queue:  squeue -u \$USER -o '%.14i %.9P %.18j %.8T %.10M %.20E'"
echo "Logs:   pepmatch-rerun-<id>.out, length-build-$build.out, length-timing-$fig5.out"
