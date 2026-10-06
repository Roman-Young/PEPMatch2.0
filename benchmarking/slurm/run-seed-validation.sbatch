#!/bin/bash
#SBATCH --job-name=seed-validate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# Large-scale recall validation of the reduced n+1-seed indel search (PEPMATCH_SEED_STRATEGY
# =first) vs full coverage and the proteome-wide oracle, at 100k for BOTH 1- and 2-indel.
#
# Runs on a compute node so nothing has to be babysat on the login node. Submit from the repo
# root while on the experiment/seed-reduction branch:
#
#   sbatch benchmarking/slurm/run-seed-validation.sbatch
#
# It (1) builds the gated engine (flag defaults OFF, so behavior is identical until a script
# sets the env var), then (2) runs both validators and both classifiers. Reads the 100k
# query + expected sets already on the cluster from the regen. PEPMatch-only against the
# precomputed oracle TSVs (no brute force), so it is fast: ~seconds for 1-indel, a few
# minutes for 2-indel. A generous 4h wall covers index (re)builds with margin.
#
# PASS looks like, per dataset:
#   lost 0 / spurious 0 / oracle missed 0 / oracle extra 0  ->  100% recall, result-identical.

REPO="$HOME/projects/pepmatch/PEPMatch2.0"
cd "$REPO"
export PATH="$HOME/.cargo/bin:$PATH"     # cargo, for the maturin build on the compute node
source "$HOME/projects/pepmatch/env.sh"  # activates the project venv (python3.10)

echo "=== host $(hostname)   branch $(git branch --show-current)   started $(date -Is) ==="
echo "=== building the gated engine (fast: target/ is warm) ==="
maturin develop --release || { echo "FATAL: engine build failed" >&2; exit 1; }

# STRATEGY defaults to the #1+#2 candidate (rarest). Override with STRATEGY=first to check
# the plain n+1 reduction instead.
STRATEGY="${STRATEGY:-rarest}"
status=0
for spec in "synth-1indel 1" "synth-2indel 2"; do
  set -- $spec
  prefix=$1; n=$2
  echo
  echo "################### ${n}-indel 100k  (strategy=${STRATEGY}) ###################"
  python3 benchmarking/scripts/validate-seed-reduction.py --prefix "$prefix" --size 100k --indels "$n" --strategy "$STRATEGY" \
    || { echo ">>> ${n}-indel VALIDATE returned nonzero (see rows above)"; status=1; }
  python3 benchmarking/scripts/classify-seed-survivors.py --prefix "$prefix" --size 100k --indels "$n" || true
done

echo
if [[ "$status" -eq 0 ]]; then
  echo "=== ALL VALIDATIONS PASSED: reduced n+1-seed search is 100% recall at 100k. ==="
else
  echo "=== ONE OR MORE VALIDATIONS FAILED -- do NOT PR; see the per-dataset rows above. ==="
fi
echo "=== finished $(date -Is) ==="
exit "$status"
