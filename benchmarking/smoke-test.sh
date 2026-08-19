#!/bin/bash
# Pre-submission smoke test. Run this ONCE before sbatch'ing the real benchmark.
#
#   bash benchmarking/smoke-test.sh
#
# On the LJI cluster, run it inside an interactive allocation so it doesn't hog a
# login node:
#   srun --cpus-per-task=8 --mem=32G --time=1:00:00 --pty bash benchmarking/smoke-test.sh
#
# It checks, in order:
#   1. every binary and python dep is present, and the Rust engine imports
#   2. all harness files parse
#   3. DIAMOND/MMseqs2 report 1-based coordinates (the one unverified correctness risk)
#   4. a full CEDAR run (50 queries) actually completes for every method
# Anything that fails here would otherwise fail at 2am with nobody watching.

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(dirname "$HERE")"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
FAILED=0

step() { echo; echo "=== $* ==="; }
ok()   { echo "  PASS: $*"; }
bad()  { echo "  FAIL: $*"; FAILED=1; }

step "1. environment"
for tool in blastp makeblastdb diamond mmseqs python3; do
  if command -v "$tool" >/dev/null 2>&1; then
    ok "$tool -> $(command -v "$tool")"
  else
    bad "$tool NOT on PATH (source env.sh first?)"
  fi
done

if python3 -c 'import pandas, pyarrow, polars' 2>/dev/null; then
  ok "python deps (pandas, pyarrow, polars)"
else
  bad "missing python deps -- pip install pandas pyarrow polars"
fi

if python3 -c 'import pepmatch; from pepmatch import _rs' 2>/dev/null; then
  ok "PEPMatch Rust engine imports (ENGINE OK)"
else
  bad "PEPMatch engine not importable -- run: maturin develop --release"
fi

step "2. harness files parse"
if python3 -m py_compile \
    "$HERE/benchmarking.py" "$HERE/verify_coordinates.py" \
    "$HERE/methods/_shell.py" "$HERE/methods/blast.py" "$HERE/methods/diamond.py" \
    "$HERE/methods/mmseqs2.py" "$HERE/methods/brute_force.py" \
    "$REPO/pepmatch/benchmarker.py" 2>&1; then
  ok "all harness files compile"
else
  bad "syntax error in the harness"
fi

if python3 -c "import json; json.load(open('$HERE/benchmarking_parameters.json'))" 2>/dev/null; then
  ok "benchmarking_parameters.json is valid"
else
  bad "benchmarking_parameters.json is malformed"
fi

step "3. coordinate convention (DIAMOND / MMseqs2 unverified until now)"
if python3 "$HERE/verify_coordinates.py"; then
  ok "all tools agree on 1-based Index start"
else
  bad "coordinate mismatch -- a tool's recall would be a FALSE 0%"
fi

step "4. full CEDAR run (50 queries, every method)"
if [[ $FAILED -eq 1 ]]; then
  echo "  SKIPPED -- fix the failures above first."
else
  cd "$HERE"
  if python3 benchmarking.py -b cedar_indel -p "$THREADS"; then
    echo
    if grep -q "FAILED\|SKIPPED" cedar_indel_benchmarking.tsv 2>/dev/null; then
      bad "a method reported FAILED/SKIPPED -- see the table + logs/ directory"
    else
      ok "every method produced results"
    fi
    echo
    echo "  --- table ---"
    column -t -s $'\t' cedar_indel_benchmarking.tsv 2>/dev/null || cat cedar_indel_benchmarking.tsv
    echo
    echo "  Sanity check: PEPMatch and Brute Force should both be 100% recall."
    echo "  Note whether MMseqs2 moved off its previous 55.6% -- masking is now OFF,"
    echo "  and if it moved materially the old numbers understated it."
  else
    bad "CEDAR run crashed"
  fi
fi

echo
if [[ $FAILED -eq 0 ]]; then
  echo "ALL CHECKS PASSED -- safe to submit:"
  echo "  sbatch benchmarking/run-indel-benchmark.sbatch"
  exit 0
fi
echo "SOME CHECKS FAILED -- do not submit the overnight job yet."
exit 1
