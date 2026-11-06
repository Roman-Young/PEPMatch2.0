#!/bin/bash
#SBATCH --job-name=st-build
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#
# Step 1 of the short-target timing pipeline (submit-st-pipeline.sh): build the Rust engine
# ONCE and record its sha256 and commit in the run folder. The timing job refuses any other
# binary, so a `git pull` between the two cannot swap the engine being timed.

set -euo pipefail
REPO="$HOME/projects/pepmatch/PEPMatch2.0"
: "${RUN_DIR:?submit through submit-st-pipeline.sh, which sets RUN_DIR}"
cd "$REPO"
source "$HOME/projects/pepmatch/env.sh"

echo "=== engine build for short-target timing ==="
echo "host   : $(hostname)"
echo "branch : $(git rev-parse --abbrev-ref HEAD)  ($(git rev-parse --short HEAD))"
git diff --quiet -- pepmatch/ || { echo "FATAL: uncommitted changes under pepmatch/; the timed engine must be a commit." >&2; exit 1; }
grep -q 'fn rs_short_proteins' pepmatch/rs-engine/src/lib.rs || { echo "FATAL: no rs_short_proteins in lib.rs; wrong branch." >&2; exit 1; }

maturin develop --release || { echo "FATAL: engine build failed" >&2; exit 1; }
SO=$(python3 -c 'from pepmatch import _rs; print(_rs.__file__)')
python3 -c 'import polars; print("polars", polars.__version__)'
sha256sum "$SO" | awk '{print $1}' > "$RUN_DIR/engine-sha256.txt"
git rev-parse HEAD > "$RUN_DIR/engine-commit.txt"
echo "engine : $SO"
echo "sha256 : $(cat "$RUN_DIR/engine-sha256.txt")"
echo "commit : $(cat "$RUN_DIR/engine-commit.txt")"
