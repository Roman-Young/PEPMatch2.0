#!/bin/bash
#SBATCH --job-name=pepmatch-repeat
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --exclusive
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# Repeats PEPMatch's search REPEATS times per dataset on a node of its own and reports median,
# min, max and every run (the first run includes loading the freshly built index from disk).
# Figure 3 datasets at 32 threads, Figure 4 datasets at 64, as in the paper.
#   bash: sbatch benchmarking/slurm/run-pepmatch-repeat.sbatch

set -euo pipefail

REPO="$HOME/projects/pepmatch/PEPMatch2.0"
REPEATS="${REPEATS:-5}"
OUT="$REPO/benchmarking/results/repeat-timing"
SCRATCH="/mnt/hpcscratch/ryoung/pepmatch-repeat-${SLURM_JOB_ID:-manual}"

echo "=== PEPMatch repeated timing ($REPEATS runs per dataset) ==="
echo "job        : ${SLURM_JOB_ID:-manual}   host: $(hostname)   started: $(date -Is)"
echo "commit     : $(git -C "$REPO" rev-parse --short HEAD) ($(git -C "$REPO" rev-parse --abbrev-ref HEAD))"
source "$HOME/projects/pepmatch/env.sh"
python3 -c 'from pepmatch import _rs; import pathlib, hashlib; p = pathlib.Path(_rs.__file__); print("engine     :", p.name, hashlib.sha256(p.read_bytes()).hexdigest()[:12])'

mkdir -p "$SCRATCH/proteomes" "$SCRATCH/work" "$OUT"
cp "$REPO/benchmarking/proteomes/"*.fasta "$SCRATCH/proteomes/"
export PEPMATCH_BENCH_PROTEOME_DIR="$SCRATCH/proteomes"
cd "$SCRATCH/work"

echo; echo "############### Figure 3 datasets, 32 threads ###############"
RAYON_NUM_THREADS=32 python3 "$REPO/benchmarking/scripts/repeat-pepmatch-timing.py" --repeats "$REPEATS" \
  --datasets cosmic_indel,cedar_indel,synth_1indel_5k,synth_2indel_5k --out "$OUT/figure3-32threads.tsv"

echo; echo "############### Figure 4 datasets, 64 threads ###############"
RAYON_NUM_THREADS=64 python3 "$REPO/benchmarking/scripts/repeat-pepmatch-timing.py" --repeats "$REPEATS" \
  --datasets synth_1indel_100,synth_1indel_1k,synth_1indel_10k,synth_1indel_100k,synth_2indel_100,synth_2indel_1k,synth_2indel_10k,synth_2indel_100k \
  --out "$OUT/figure4-64threads.tsv"

echo; echo "=== RESULTS (paste back to Kairo) ==="
for f in "$OUT"/figure3-32threads.tsv "$OUT"/figure4-64threads.tsv; do
  echo; echo "--- $(basename "$f") ---"; cut -f1,3,4,5,6,8 "$f" | column -t -s $'\t'
done
echo; echo "finished   : $(date -Is)"
