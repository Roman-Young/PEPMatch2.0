#!/bin/bash
#SBATCH --job-name=plk-build
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#
# Step 1 of the per-length-k comparison (submit-plk-pipeline.sh): build BOTH engines once,
# each from a clean export of its commit, and unpack each wheel into its own folder in the
# run directory. The timing job puts one folder at a time on PYTHONPATH, so the two builds
# share the Python and polars of the environment and differ only in the pepmatch package.
# Nothing is installed into the environment and the repo's working tree is not touched.

set -euo pipefail
REPO="$HOME/projects/pepmatch/PEPMatch2.0"
: "${RUN_DIR:?submit through submit-plk-pipeline.sh, which sets RUN_DIR}"
: "${BASE_COMMIT:?}" "${HEAD_COMMIT:?}"
cd "$REPO"
source "$HOME/projects/pepmatch/env.sh"

echo "=== engine builds for the per-length-k comparison ==="
echo "host : $(hostname)"
python3 -c 'import polars, Bio; print("polars", polars.__version__)'

for variant in base head; do
  [[ $variant == base ]] && commit="$BASE_COMMIT" || commit="$HEAD_COMMIT"
  full=$(git rev-parse --verify "$commit^{commit}") || { echo "FATAL: commit $commit not found; run 'git fetch origin' first." >&2; exit 1; }
  src="$RUN_DIR/src-$variant"
  mkdir -p "$src" "$RUN_DIR/site-$variant"
  git archive "$full" | tar -x -C "$src"
  ( cd "$src" && maturin build --release --interpreter "$(command -v python3)" --out "$RUN_DIR/wheel-$variant" ) \
    || { echo "FATAL: $variant build failed" >&2; exit 1; }
  python3 -m zipfile -e "$RUN_DIR"/wheel-$variant/*.whl "$RUN_DIR/site-$variant"
  so=$(ls "$RUN_DIR/site-$variant"/pepmatch/_rs*.so)
  echo "$full" > "$RUN_DIR/$variant-commit.txt"
  sha256sum "$so" | awk '{print $1}' > "$RUN_DIR/$variant-sha256.txt"
  echo "$variant : $full  $(git log -1 --format=%s "$full")"
  echo "       engine sha256 $(cat "$RUN_DIR/$variant-sha256.txt")"
done

# base still takes a k and head must not: catches the two commits being swapped or equal.
grep -q "k=0" "$RUN_DIR/src-base/pepmatch/matcher.py" || { echo "FATAL: base has no k option; wrong BASE_COMMIT." >&2; exit 1; }
grep -q "def k_for" "$RUN_DIR/src-head/pepmatch/matcher.py" || { echo "FATAL: head has no k_for; wrong HEAD_COMMIT." >&2; exit 1; }
echo "builds ready in $RUN_DIR"
