#!/bin/bash
#SBATCH --job-name=seed-timing
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# New PEPMatch search timing under the rarest n+1-seed strategy (#1+#2), to REPLACE the
# PEPMatch column in the scaling figures. Runs PEPMatch ONLY, both the current full-coverage
# path (default) and rarest, at each size on ONE node, so the default-vs-rarest speedup is
# measured on identical hardware. Brute Force is NOT rerun: its times stand; only PEPMatch
# gets faster, and PEPMatch is already hundreds of x below Brute Force so a node change on
# the (tiny) PEPMatch side does not move the vs-BF story.
#
#   sbatch benchmarking/slurm/run-seed-timing.sbatch
#
# Recall equivalence is validated separately (run-seed-validation.sbatch, STRATEGY=rarest).
# Run that FIRST and confirm PASS before trusting any timing here.

REPO="$HOME/projects/pepmatch/PEPMatch2.0"
cd "$REPO"
export PATH="$HOME/.cargo/bin:$PATH"
source "$HOME/projects/pepmatch/env.sh"
THREADS="${SLURM_CPUS_PER_TASK:-64}"
export SIZES="${SIZES:-100 1k 10k 100k}"

echo "=== host $(hostname)  threads $THREADS  started $(date -Is) ==="
maturin develop --release || { echo "FATAL: engine build failed" >&2; exit 1; }

# Scratch proteome so index (re)builds land on fast local scratch, not the home quota
# (mirrors run-synth-scaling.sbatch). Preprocessing is charged separately from search by the
# harness, so index-build time never contaminates the search numbers reported below.
SCRATCH="/mnt/hpcscratch/ryoung/pepmatch-seed-timing-${SLURM_JOB_ID:-manual}"
mkdir -p "$SCRATCH/proteomes" "$SCRATCH/work"
cp "$REPO/benchmarking/proteomes/human.fasta" "$SCRATCH/proteomes/"
export PEPMATCH_BENCH_PROTEOME_DIR="$SCRATCH/proteomes"
cd "$SCRATCH/work"

for INDELS in 1 2; do
  for size in $SIZES; do
    ds="synth_${INDELS}indel_${size}"
    for strat in default rarest; do
      if [[ "$strat" == "rarest" ]]; then export PEPMATCH_SEED_STRATEGY=rarest; else unset PEPMATCH_SEED_STRATEGY; fi
      echo ">>> timing $ds  ($strat)  $(date -Is)"
      python3 "$REPO/benchmarking/benchmarking.py" -b "$ds" --methods PEPMatch -p "$THREADS" \
        --results-subdir "seed-timing-$strat" > /dev/null 2>&1 \
        || echo "    !! run failed: $ds ($strat)"
    done
  done
done
unset PEPMATCH_SEED_STRATEGY

echo
echo "=== PEPMatch SEARCH TIME (s): default full coverage vs rarest n+1 seeds ==="
python3 - <<'PYEOF'
import csv, os
RES = os.path.expanduser("~/projects/pepmatch/PEPMatch2.0/benchmarking/results")
sizes = os.environ.get("SIZES", "100 1k 10k 100k").split()

def search_s(sub, ds):
    p = f"{RES}/{sub}/{ds}_benchmarking.tsv"
    try:
        for r in csv.DictReader(open(p), delimiter="\t"):
            if r["Method"] == "PEPMatch":
                return float(r["Searching (s)"])
    except FileNotFoundError:
        return None
    return None

for n in (1, 2):
    print(f"\n  {n}-indel      N     default(s)    rarest(s)   speedup")
    for size in sizes:
        ds = f"synth_{n}indel_{size}"
        d = search_s("seed-timing-default", ds)
        r = search_s("seed-timing-rarest", ds)
        if d is not None and r is not None:
            print(f"          {size:>7}   {d:10.2f}   {r:10.2f}   {d/r:6.2f}x")
        else:
            print(f"          {size:>7}   (missing tsv)")
print("\nUse the rarest column as the new PEPMatch search time in the scaling figures;")
print("Brute Force times are unchanged. The speedup column is default-vs-rarest on one node.")
PYEOF
echo
echo "=== finished $(date -Is) ==="
