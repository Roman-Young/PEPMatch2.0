#!/bin/bash
#SBATCH --job-name=pepmatch-rerun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# Re-times PEPMatch ONLY on the paper's datasets after the post-profiling speedups (#35, #36,
# the indel labeling/annotation fixes). Aligners, brute force and the index build did not
# change, so their existing tables stay as they are; this writes PEPMatch rows to a separate
# results/<SUBDIR>/ so nothing already published gets overwritten.
#
#   DATASETS  space-separated benchmark names (benchmarking_parameters.json)
#   SUBDIR    results subdirectory for the new tables
#   MEMORY    1 = also the memory table (PEPMatch's peak memory changed with #35)

set -euo pipefail

REPO="$HOME/projects/pepmatch/PEPMatch2.0"
DATASETS="${DATASETS:?set DATASETS}"
SUBDIR="${SUBDIR:?set SUBDIR}"
MEMORY="${MEMORY:-0}"
SCRATCH="/mnt/hpcscratch/ryoung/pepmatch-rerun-${SLURM_JOB_ID:-manual}"
RESULTS="$REPO/benchmarking/results"

echo "=== PEPMatch re-time ($SUBDIR) ==="
echo "job        : ${SLURM_JOB_ID:-manual}"
echo "host       : $(hostname)"
echo "started    : $(date -Is)"
echo "commit     : $(git -C "$REPO" rev-parse --short HEAD) ($(git -C "$REPO" rev-parse --abbrev-ref HEAD))"
source "$HOME/projects/pepmatch/env.sh"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
echo "threads    : $THREADS"
echo "datasets   : $DATASETS"

python3 -c 'import pandas, pyarrow, polars' || {
  echo "FATAL: missing python deps (pandas, pyarrow, polars)." >&2; exit 1; }
python3 - "$REPO" <<'EOF' || { echo "FATAL: the engine is stale or not the checkout's build." >&2; exit 1; }
import sys
from pathlib import Path
import pepmatch
from pepmatch import _rs
repo = Path(sys.argv[1]).resolve()
pkg = Path(pepmatch.__file__).resolve().parent
src = pkg / 'rs-engine' / 'src' / 'match.rs'
so = Path(_rs.__file__).resolve()
assert pkg == repo / 'pepmatch', f'imported pepmatch from {pkg}, not {repo}'
assert so.stat().st_mtime > src.stat().st_mtime, 'engine older than match.rs: run maturin develop --release'
assert 'impl Walk' in src.read_text(), 'match.rs lacks the single-walk indel extension (#36): wrong branch'
print(f'  ENGINE OK  {so.name}')
EOF

mkdir -p "$SCRATCH/proteomes" "$SCRATCH/work" "$RESULTS/logs"
cp "$REPO/benchmarking/proteomes/"*.fasta "$SCRATCH/proteomes/"
export PEPMATCH_BENCH_PROTEOME_DIR="$SCRATCH/proteomes"
export PEPMATCH_BENCH_LOGDIR="$RESULTS/logs"
cd "$SCRATCH/work"

echo; echo "=== pre-flight: coordinate convention ==="
python3 "$REPO/benchmarking/scripts/verify_coordinates.py" || {
  echo "FATAL: coordinate check failed -- recall numbers would be meaningless." >&2; exit 1; }

for dataset in $DATASETS; do
  query=$(python3 -c "import json,sys; print(json.load(open('$REPO/benchmarking/benchmarking_parameters.json'))['datasets']['$dataset']['query'])")
  if [[ ! -f "$REPO/benchmarking/$query" ]]; then
    echo "SKIP $dataset: $query not generated on this cluster"; continue
  fi
  echo; echo "############### TIMING: $dataset  ($(date -Is)) ###############"
  python3 "$REPO/benchmarking/benchmarking.py" -b "$dataset" --methods PEPMatch -p "$THREADS" --results-subdir "$SUBDIR"
  if [[ "$MEMORY" == "1" ]]; then
    echo; echo "############### MEMORY: $dataset  ($(date -Is)) ###############"
    python3 "$REPO/benchmarking/benchmarking.py" -b "$dataset" --methods PEPMatch --memory-only -p "$THREADS" --results-subdir "$SUBDIR"
  fi
done

echo; echo "=== RESULTS (paste back to Kairo) ==="
for f in "$RESULTS/$SUBDIR"/*.tsv; do
  echo; echo "--- $(basename "$f") ---"
  column -t -s $'\t' "$f" 2>/dev/null || cat "$f"
  grep -q "FAILED\|SKIPPED" "$f" && echo "  !! FAILED/SKIPPED present -- see $PEPMATCH_BENCH_LOGDIR"
done
echo; echo "finished   : $(date -Is)"
