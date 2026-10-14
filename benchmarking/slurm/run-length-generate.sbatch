#!/bin/bash
#SBATCH --job-name=length-generate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
# 36h. The slowest task, 2-indel 8-mers at 10k queries, projects to ~9.3h at 64 cores
# (--cost-estimate, incl. its 3.3x heavy-tail factor); x2 node variance plus the audit
# still fits. Every other task finishes in minutes to ~2h.
#SBATCH --time=36:00:00
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err
#
# Step 2 of the per-query-length timing pipeline. ONE ARRAY TASK = ONE (indels, length)
# query set. Submitted by submit-length-pipeline.sh, which passes TASKS and the --array
# range; do not submit by hand without them.
#
# Parallel across nodes on purpose: generation is NOT timed, so node-to-node variance does
# not matter here. (The timing job is the opposite: one node, sequential.) Each task
# writes only its own prefixed files and verifies from its own scratch dir, so tasks never
# touch each other's outputs or .pepidx files.
#
# Each set: --length-dist uniform --min-len L --max-len L (single length; the default
# `cedar` distribution has no mass at 19/21/24 and would exit). Sizes 100/1k/10k are
# nested prefixes: 1k + 10k are timed, 100 is the self-certification subset.

set -euo pipefail

REPO="$HOME/projects/pepmatch/PEPMatch2.0"
OUT="$REPO/benchmarking/results/per-length"
PREFIX="${PREFIX:-synth-len}"
NMAX="${NMAX:-10000}"
SIZES="${SIZES:-100,1000,10000}"
AUDIT="${AUDIT:-200}"
CERT_N="${CERT_N:-100}"
read -r -a TASK_LIST <<< "${TASKS:?TASKS is unset -- submit through submit-length-pipeline.sh}"
TASK="${TASK_LIST[${SLURM_ARRAY_TASK_ID:?not an array task}]}"
INDELS="${TASK%%:*}"
LEN="${TASK##*:}"
LL=$(printf '%02d' "$LEN")
SET="${PREFIX}-${INDELS}indel-L${LL}"

echo "=== per-length generation: $SET ==="
echo "job        : ${SLURM_ARRAY_JOB_ID:-?}_${SLURM_ARRAY_TASK_ID}"
echo "host       : $(hostname)"
echo "started    : $(date -Is)"

source "$HOME/projects/pepmatch/env.sh"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
export PEPMATCH_SEED_STRATEGY=rarest

SO=$(python3 -c 'from pepmatch import _rs; print(_rs.__file__)')
if [[ "$(sha256sum "$SO" | awk '{print $1}')" != "$(cat "$OUT/engine-sha256.txt")" ]]; then
  echo "FATAL: engine is not the binary length-build recorded -- rerun the pipeline." >&2
  exit 1
fi

GEN="$REPO/benchmarking/scripts/build-synthetic-indel-benchmark.py"
COMMON=(--indels "$INDELS" --min-len "$LEN" --max-len "$LEN" --length-dist uniform
        --n-max "$NMAX" --sizes "$SIZES" --threads "$THREADS"
        --out-prefix "$SET" --name-prefix "${SET//-/_}")

echo; echo "=== projected cost ==="
python3 "$GEN" "${COMMON[@]}" --cost-estimate

echo; echo "=== generating ==="
python3 "$GEN" "${COMMON[@]}" --audit-sample "$AUDIT"

echo; echo "=== self-certifying (PEPMatch + Brute Force must both score 100.0) ==="
# The recall gate builds a .pepidx in the CWD; a per-task scratch dir keeps parallel
# tasks from reading each other's half-written index.
SCRATCH="/mnt/hpcscratch/ryoung/pepmatch-length-gen-${SLURM_ARRAY_JOB_ID:-x}-${SLURM_ARRAY_TASK_ID}"
mkdir -p "$SCRATCH/proteomes" "$SCRATCH/work"
cp "$REPO/benchmarking/proteomes/human.fasta" "$SCRATCH/proteomes/"
export PEPMATCH_BENCH_PROTEOME_DIR="$SCRATCH/proteomes"
(
  cd "$SCRATCH/work"
  python3 "$REPO/benchmarking/scripts/verify-synthetic-dataset.py" \
    --prefix "$SET" --run-tools-on "$CERT_N" --threads "$THREADS"
) || { echo "FATAL: $SET failed self-certification." >&2; exit 1; }
rm -rf "$SCRATCH"

echo
echo "finished   : $(date -Is)"
