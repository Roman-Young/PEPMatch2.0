#!/bin/bash
#SBATCH --job-name=annot-validate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
# 12h. Projected well under 2h; the slowest point is the old code on the 2-indel 100k set
# (3 repeats + one full-table run). A generous wall is free: the job ends when it finishes.
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# PHASE 1: prove the indel-annotation speedup changes nothing but time, at full scale.
#
#   git checkout bench/annotation-speedup
#   sbatch benchmarking/slurm/run-annotation-validation.sbatch
#
# For every dataset behind Figures 3 and 4, runs the ORIGINAL matcher.py ("old", from
# OLD_REF) and the new one against the SAME engine build, on ONE node, and requires a
# byte-identical full output table (Indel Positions included) and recall 100.0 for both.
# All logic and the pass/fail verdict: scripts/validate-annotation-speedup.py.
#
# Thread counts are the figures' own: 32 for the aligner sets (Figure 3), 64 for the
# scaling sets (Figure 4). Old/new order alternates per dataset so slow drift on the
# node cannot favour one version.

set -euo pipefail

REPO="$HOME/projects/pepmatch/PEPMatch2.0"
OLD_REF="${OLD_REF:-352ce92}"     # IEDB master before the speedup
DATASETS="${DATASETS:-cosmic_indel:32 cedar_indel:32 synth_2indel_5k:32 \
synth_1indel_100:64 synth_1indel_1k:64 synth_1indel_10k:64 synth_1indel_100k:64 \
synth_2indel_100:64 synth_2indel_1k:64 synth_2indel_10k:64 synth_2indel_100k:64}"
REPEATS="${REPEATS:-3}"
OUT="$REPO/benchmarking/results/annotation-validation-${SLURM_JOB_ID:-manual}"
VALIDATE="$REPO/benchmarking/scripts/validate-annotation-speedup.py"

echo "=== Phase 1: annotation speedup validation ==="
echo "job        : ${SLURM_JOB_ID:-manual}"
echo "host       : $(hostname)"
echo "started    : $(date -Is)"

cd "$REPO"
export PATH="$HOME/.cargo/bin:$PATH"
source "$HOME/projects/pepmatch/env.sh"

branch="$(git rev-parse --abbrev-ref HEAD)"
[[ "$branch" == "bench/annotation-speedup" ]] || {
  echo "FATAL: on branch '$branch', expected bench/annotation-speedup." >&2; exit 1; }
git diff --quiet -- pepmatch/ || {
  echo "FATAL: uncommitted changes under pepmatch/ -- the tested code must be a commit." >&2
  exit 1; }
git cat-file -e "$OLD_REF:pepmatch/matcher.py" || {
  echo "FATAL: $OLD_REF:pepmatch/matcher.py not found -- cannot build the old version." >&2
  exit 1; }
echo "commit     : $(git rev-parse --short HEAD)   old matcher from $OLD_REF"

maturin develop --release || { echo "FATAL: engine build failed" >&2; exit 1; }
SO="$(python3 -c 'from pepmatch import _rs; print(_rs.__file__)')"

SCRATCH="/mnt/hpcscratch/ryoung/pepmatch-annot-validate-${SLURM_JOB_ID:-manual}"
mkdir -p "$SCRATCH/proteomes" "$SCRATCH/work" "$SCRATCH/old_pkg" "$OUT"
cp "$REPO/benchmarking/proteomes/human.fasta" "$SCRATCH/proteomes/"
export PEPMATCH_BENCH_PROTEOME_DIR="$SCRATCH/proteomes"

# The old version = this exact package (same compiled engine) with the original
# matcher.py swapped in. Imported via PYTHONPATH, which precedes the editable install.
cp -r "$REPO/pepmatch" "$SCRATCH/old_pkg/"
git show "$OLD_REF:pepmatch/matcher.py" > "$SCRATCH/old_pkg/pepmatch/matcher.py"
OLD_SO="$(ls "$SCRATCH"/old_pkg/pepmatch/_rs*.so)"
[[ "$(sha256sum "$SO" | cut -c1-64)" == "$(sha256sum "$OLD_SO" | cut -c1-64)" ]] || {
  echo "FATAL: old and new packages do not share one engine build." >&2; exit 1; }

# Source fingerprint of the code under test. After the PR merges, the same hashes on the
# merged files mean every number from this run still describes the released code.
{
  echo "commit	$(git rev-parse HEAD)"
  echo "old_ref	$OLD_REF"
  for f in pepmatch/matcher.py pepmatch/benchmarker.py pepmatch/preprocessor.py \
           pepmatch/rs-engine/src/*.rs; do
    echo "$f	$(sha256sum "$f" | cut -c1-64)"
  done
  echo "engine_so	$(sha256sum "$SO" | cut -c1-64)"
  echo "host	$(hostname)"
} > "$OUT/source-fingerprint.txt"

echo; echo "=== pre-flight: every dataset's files exist ==="
python3 - "$DATASETS" <<'PYEOF'
import json, sys
from pathlib import Path
bench = Path.home() / 'projects/pepmatch/PEPMatch2.0/benchmarking'
cfg = json.load(open(bench / 'benchmarking_parameters.json'))['datasets']
missing = []
for item in sys.argv[1].split():
  name = item.split(':')[0]
  if name not in cfg:
    missing.append(f'{name}: not registered')
    print(f'  MISSING: {name} (not registered)')
    continue
  gone = [cfg[name][k] for k in ('query', 'expected') if not (bench / cfg[name][k]).exists()]
  missing += [f'{name}: {g}' for g in gone]
  print(f'  MISSING: {name} {gone}' if gone else f'  ok: {name}')
if missing:
  sys.exit('FATAL: ' + '; '.join(missing))
PYEOF

cd "$SCRATCH/work"
i=0
for item in $DATASETS; do
  ds="${item%%:*}"; threads="${item##*:}"
  if (( i % 2 == 0 )); then order="old new"; else order="new old"; fi
  i=$((i + 1))
  echo; echo "############### $ds ($threads threads)  $(date -Is) ###############"
  for version in $order; do
    if [[ "$version" == "old" ]]; then pp="$SCRATCH/old_pkg"; else pp="$REPO"; fi
    PYTHONPATH="$pp" python3 -u "$VALIDATE" run --version "$version" --dataset "$ds" \
      --threads "$threads" --repeats "$REPEATS" --out "$OUT" --tables "$SCRATCH/tables"
  done
done

echo; echo "=== VERDICT ==="
python3 "$VALIDATE" compare --out "$OUT"

echo
echo "=== paste these back to Kairo ==="
echo "----- source-fingerprint.txt"; cat "$OUT/source-fingerprint.txt"
echo "----- validation-summary.tsv"; cat "$OUT/validation-summary.tsv"
echo
echo "finished   : $(date -Is)"
echo "results    : $OUT"
echo "tables     : $SCRATCH/tables  (full outputs, old and new; scratch auto-purges after 5 days)"
