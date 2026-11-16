#!/bin/bash
#SBATCH --job-name=plk-build
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#
# Step 1 of the per-length-k comparison (submit-plk-pipeline.sh): build base and head once,
# each from a clean export of its commit, unpack each wheel into its own folder in the run
# directory, and make one copy of head per k ceiling in CAPS. The timing job puts one folder at a time on PYTHONPATH, so the two builds
# share the Python and polars of the environment and differ only in the pepmatch package.
# Nothing is installed into the environment and the repo's working tree is not touched.

set -euo pipefail
REPO="$HOME/projects/pepmatch/PEPMatch2.0"
: "${RUN_DIR:?submit through submit-plk-pipeline.sh, which sets RUN_DIR}"
: "${BASE_COMMIT:?}" "${HEAD_COMMIT:?}" "${CAPS:?}"
cd "$REPO"
source "$HOME/projects/pepmatch/env.sh"

echo "=== engine builds for the per-length-k comparison ==="
echo "host : $(hostname)"
python3 -c 'import polars, Bio; print("polars", polars.__version__)'

# base: one k per batch, built from BASE_COMMIT.
# capN: HEAD_COMMIT with the k rule's ceiling set to N (k = min(N, length // (edits + 1))).
# The engine is the same for every cap; only that one line of matcher.py differs, so head is
# built once and each cap gets a copy of the build with the line changed.
build() {   # build <commit> <name>
  local full src
  full=$(git rev-parse --verify "$1^{commit}") || { echo "FATAL: commit $1 not found; run 'git fetch origin' first." >&2; exit 1; }
  src="$RUN_DIR/src-$2"
  mkdir -p "$src" "$RUN_DIR/site-$2"
  git archive "$full" | tar -x -C "$src"
  ( cd "$src" && maturin build --release --interpreter "$(command -v python3)" --out "$RUN_DIR/wheel-$2" ) \
    || { echo "FATAL: $2 build failed" >&2; exit 1; }
  python3 -m zipfile -e "$RUN_DIR"/wheel-$2/*.whl "$RUN_DIR/site-$2"
  echo "$full" > "$RUN_DIR/$2-commit.txt"
  echo "$2 : $full  $(git log -1 --format=%s "$full")"
  echo "       engine sha256 $(sha256sum "$RUN_DIR/site-$2"/pepmatch/_rs*.so | awk '{print $1}')"
}
build "$BASE_COMMIT" base
build "$HEAD_COMMIT" head
grep -q "k=0" "$RUN_DIR/src-base/pepmatch/matcher.py" || { echo "FATAL: base has no k option; wrong BASE_COMMIT." >&2; exit 1; }
RULE='  return min(MAX_INDEX_K, max(2, length // (edits + 1)))'
grep -qxF "$RULE" "$RUN_DIR/site-head/pepmatch/matcher.py" || { echo "FATAL: head's k rule line not found; wrong HEAD_COMMIT." >&2; exit 1; }

for cap in $CAPS; do
  site="$RUN_DIR/site-cap$cap"
  cp -r "$RUN_DIR/site-head" "$site"
  python3 - "$site/pepmatch/matcher.py" "$cap" <<'PY'
import sys
path, cap = sys.argv[1], int(sys.argv[2])
rule = '  return min(MAX_INDEX_K, max(2, length // (edits + 1)))'
text = open(path).read()
assert text.count(rule) == 1 and 2 <= cap <= 16
open(path, 'w').write(text.replace(rule, f'  return min({cap}, max(2, length // (edits + 1)))'))
PY
  # Run from inside the copy: python -c looks in the current folder first, and the repo
  # root (where this job starts) has its own pepmatch/.
  (cd "$site" && PYTHONPATH="$site" python3 -c "
import pepmatch.matcher as m, sys
assert m.__file__.startswith('$site'), m.__file__
got = [m.k_for(40, 0), m.k_for(9, 0), m.k_for(40, 3)]
want = [$cap, min($cap, 9), min($cap, 10)]
assert got == want, (got, want)
print('cap$cap : k for a 40-mer exact / 9-mer exact / 40-mer at 3 mm =', got)")
done
echo "builds ready in $RUN_DIR"
