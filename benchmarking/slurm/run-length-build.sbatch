#!/bin/bash
#SBATCH --job-name=length-build
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# Step 1 of the per-query-length timing pipeline (submit-length-pipeline.sh).
#
# Builds the Rust engine ONCE and records its sha256. Every later job (36 generation
# tasks + the timing job) checks that hash before running. Two reasons this is its own job:
#   * maturin writes the compiled engine into the shared repo, so 36 array tasks each
#     rebuilding would overwrite each other's .so mid-import;
#   * it pins the binary: generation certifies recall with it and timing measures it, so
#     a `git pull` between the two cannot silently swap the engine being timed.
#
# The engine must be the rarest-seed build that produced Figure 4: branch
# feature/per-length-timing (off experiment/seed-reduction), run with
# PEPMATCH_SEED_STRATEGY=rarest. This job refuses any other source tree.

set -euo pipefail

REPO="$HOME/projects/pepmatch/PEPMatch2.0"
OUT="$REPO/benchmarking/results/per-length"
cd "$REPO"
export PATH="$HOME/.cargo/bin:$PATH"
source "$HOME/projects/pepmatch/env.sh"

echo "=== engine build for per-length timing ==="
echo "host       : $(hostname)"
echo "branch     : $(git rev-parse --abbrev-ref HEAD)  ($(git rev-parse --short HEAD))"
git diff --quiet -- pepmatch/ || {
  echo "FATAL: uncommitted changes under pepmatch/ -- the timed engine must be a commit." >&2
  exit 1; }
grep -q 'Ok("rarest") => rarest_pigeonhole_seeds' pepmatch/rs-engine/src/match.rs || {
  echo "FATAL: match.rs has no PEPMATCH_SEED_STRATEGY=rarest switch -- wrong branch." >&2
  exit 1; }

maturin develop --release || { echo "FATAL: engine build failed" >&2; exit 1; }

SO=$(python3 -c 'from pepmatch import _rs; print(_rs.__file__)')
mkdir -p "$OUT"
sha256sum "$SO" | awk '{print $1}' > "$OUT/engine-sha256.txt"
git rev-parse HEAD > "$OUT/engine-commit.txt"
echo "engine     : $SO"
echo "sha256     : $(cat "$OUT/engine-sha256.txt")"
echo "commit     : $(cat "$OUT/engine-commit.txt")"
