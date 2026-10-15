#!/bin/bash
#SBATCH --job-name=repeat-timing
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --exclusive
#SBATCH --mem=200G
#SBATCH --time=10:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# Figures 3 and 4 timing, every tool measured the same way: REPEATS x (build, then search),
# median search time reported (scripts/repeat-timing.py). One node of its own for the whole
# job, so every tool and dataset in a figure ran on the same node.
#   Figure 3: PEPMatch, BLAST (short and default), DIAMOND, MMseqs2 at 32 threads
#   Figure 4: PEPMatch at 64 threads (brute force is not rerun: one run per point takes
#             up to 54 h and its scan cost barely varies)
# ~1.5 h, nearly all of it BLAST.
#
# Indexes and databases go on the network scratch (WekaFS), as in every other paper run;
# every tool reads its freshly built index or database from the same file system.
#   bash: sbatch benchmarking/slurm/run-repeat-timing.sbatch

set -euo pipefail

REPO="$HOME/projects/pepmatch/PEPMatch2.0"
REPEATS="${REPEATS:-5}"
OUT="$REPO/benchmarking/results/repeat-all-tools"
ALIGNERS='PEPMatch,blast (short),blast (default),diamond,mmseqs2'

echo "=== repeated timing, all tools ($REPEATS x build+search) ==="
echo "job        : ${SLURM_JOB_ID:-manual}   host: $(hostname)   started: $(date -Is)"
echo "cpu        : $(lscpu | sed -n 's/^Model name: *//p')   (${SLURM_CPUS_ON_NODE:-?} cpus allocated)"
echo "ram        : $(free -g | awk '/^Mem:/ {print $2 " GB"}')"
echo "commit     : $(git -C "$REPO" rev-parse --short HEAD) ($(git -C "$REPO" rev-parse --abbrev-ref HEAD))"
source "$HOME/projects/pepmatch/env.sh"

for tool in blastp makeblastdb diamond mmseqs python3; do
  command -v "$tool" >/dev/null 2>&1 || { echo "FATAL: '$tool' not on PATH after sourcing env.sh." >&2; exit 1; }
done
echo "blast      : $(blastp -version | head -1)"
echo "diamond    : $(diamond version)"
echo "mmseqs     : $(mmseqs version)"
python3 -c 'import pandas, pyarrow, polars' || { echo "FATAL: missing python deps (pandas, pyarrow, polars)." >&2; exit 1; }
python3 - "$REPO" <<'EOF' || { echo "FATAL: the engine is stale or not the checkout's build." >&2; exit 1; }
import hashlib, sys
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
print(f'engine     : {so.name} {hashlib.sha256(so.read_bytes()).hexdigest()[:12]}  OK')
EOF

WORK="/mnt/hpcscratch/ryoung/pepmatch-repeat-${SLURM_JOB_ID:-manual}"
mkdir -p "$WORK/proteomes" "$WORK/work" "$OUT"
echo "workdir    : $WORK"
df -hT "$WORK" | tail -1 | awk '{print "filesystem : " $2 " on " $1 ", " $5 " free"}'
cp "$REPO/benchmarking/proteomes/"*.fasta "$WORK/proteomes/"
export PEPMATCH_BENCH_PROTEOME_DIR="$WORK/proteomes"
export PEPMATCH_BENCH_LOGDIR="$OUT/logs"
mkdir -p "$PEPMATCH_BENCH_LOGDIR"
cd "$WORK/work"

echo; echo "=== pre-flight: coordinate convention ==="
python3 "$REPO/benchmarking/scripts/verify_coordinates.py" || {
  echo "FATAL: coordinate check failed -- recall numbers would be meaningless." >&2; exit 1; }

# Separate processes per thread count: PEPMatch's engine reads RAYON_NUM_THREADS once.
status=0
echo; echo "############### Figure 3: all tools, 32 threads  ($(date -Is)) ###############"
python3 "$REPO/benchmarking/scripts/repeat-timing.py" --repeats "$REPEATS" --threads 32 \
  --datasets cosmic_indel,cedar_indel,synth_2indel_5k --methods "$ALIGNERS" --outdir "$OUT/figure3" || status=1

echo; echo "############### Figure 4: PEPMatch, 64 threads  ($(date -Is)) ###############"
python3 "$REPO/benchmarking/scripts/repeat-timing.py" --repeats "$REPEATS" --threads 64 \
  --datasets synth_1indel_100,synth_1indel_1k,synth_1indel_10k,synth_1indel_100k,synth_2indel_100,synth_2indel_1k,synth_2indel_10k,synth_2indel_100k \
  --methods PEPMatch --outdir "$OUT/figure4" || status=1

echo; echo "=== RESULTS (paste back to Kairo) ==="
for f in "$OUT"/figure3/*_benchmarking.tsv "$OUT"/figure4/*_benchmarking.tsv; do
  echo; echo "--- $(basename "$(dirname "$f")")/$(basename "$f") ---"
  cut -f1,2,4,7,8,9,12 "$f" | column -t -s $'\t'
done
echo; echo "finished   : $(date -Is)"
exit $status
