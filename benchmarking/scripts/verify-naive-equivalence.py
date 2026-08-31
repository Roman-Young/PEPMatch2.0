#!/usr/bin/env python3
"""Prove `brute_force_naive` returns EXACTLY the rows `brute_force` returns, and measure
what the naive sweep costs.

Two jobs, one run:

1. EQUIVALENCE (the correctness gate). The prefilter in `brute_force.py` is sound (a
   <=n-indel edit cannot touch all n+1 disjoint pieces, so at least one survives verbatim),
   which means dropping it must change runtime and NOTHING else. If these two ever disagree
   on a single row, one of them is wrong and no number from either is publishable. Compared
   as exact row sets, not counts -- equal counts with different rows is precisely the bug
   this is here to catch.

2. COST (the planning number). The naive sweep is O(proteome) per query with no candidate
   narrowing, so it is orders of magnitude slower. Its measured per-query cost is what
   decides the largest N the scaling sweeps can actually reach, and this prints the
   projection instead of leaving it to be discovered 14 hours into a cluster job.

Usage:
  verify-naive-equivalence.py --indels 1 --limit 20
  verify-naive-equivalence.py --indels 2 --limit 10 --queries queries/synth-2indel-100.fasta
"""
import argparse
import importlib.util
import random
import sys
import tempfile
import time
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH / 'methods'))


def _load_method(name):
  spec = importlib.util.spec_from_file_location(name, BENCH / 'methods' / f'{name}.py')
  mod = importlib.util.module_from_spec(spec)
  sys.modules[name] = mod          # Pool workers re-resolve the module by name after fork
  spec.loader.exec_module(mod)
  return mod


def _rows(df):
  """Row SET, so ordering and any duplicate rows are both normalised away."""
  return {tuple(r) for r in df.itertuples(index=False, name=None)}


def _check_oracle_equivalence(args):
  """Prove: generic-aligner(max_indels=1) rows == committed-1-indel-oracle rows, per
  (query, protein). If they ever differ, the unified baseline would report a recall other
  than the committed ground truth's, so this is the gate that licenses the unification."""
  import csv
  import random as _random

  def _load_oracle(name):
    spec = importlib.util.spec_from_file_location(name, BENCH.parent / 'pepmatch' / 'tests' / f'{name}.py')
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m.brute_force_search

  o_specific = _load_oracle('indel_brute_force')    # what the 1-indel ground truth used
  o_generic = _load_oracle('indel2_brute_force')    # what the naive baseline now uses

  seqs, acc, chunks = {}, None, []
  for line in open(args.proteome):
    line = line.rstrip('\n')
    if line.startswith('>'):
      if acc is not None:
        seqs[acc] = ''.join(chunks)
      parts = line[1:].split('|')
      acc = parts[1] if len(parts) >= 2 else line[1:].split()[0]
      chunks = []
    else:
      chunks.append(line)
  seqs[acc] = ''.join(chunks)

  qfile = BENCH / 'queries' / 'synth-1indel-100.fasta'
  efile = BENCH / 'expected' / 'synth-1indel-100-expected.tsv'
  qs = sorted({ln.strip() for ln in open(qfile) if ln.strip() and not ln.startswith('>')})
  exp = {}
  for r in csv.DictReader(open(efile), delimiter='\t'):
    exp.setdefault(r['Query Sequence'], set()).add(r['Protein ID'])

  rng = _random.Random(args.seed)
  sample = rng.sample(qs, min(args.limit, len(qs)))
  others = list(seqs)

  def rows(oracle, q, accs, generic):
    out = set()
    for a in accs:
      res = oracle(q, seqs[a], 1) if generic else oracle(q, seqs[a])
      for start, matched in res:
        if len(matched) != len(q):
          out.add((a, start, matched))
    return out

  print(f'oracle-equivalence check: {len(sample)} real 1-indel queries, '
        f'each vs its matches + 150 random proteins')
  mism = 0
  for q in sample:
    accs = set(exp.get(q, set())) | set(rng.sample(others, 150))
    r_spec, r_gen = rows(o_specific, q, accs, False), rows(o_generic, q, accs, True)
    if r_spec != r_gen:
      mism += 1
      if mism <= 5:
        print(f'  MISMATCH {q}: only committed={list(r_spec - r_gen)[:2]} '
              f'only generic={list(r_gen - r_spec)[:2]}')
  if mism:
    print(f'\n*** {mism} queries disagreed -- the unified baseline is NOT safe. ***')
    raise SystemExit(1)
  print('\nALL AGREE: the generic aligner (max_indels=1) reproduces the committed 1-indel '
        'oracle exactly.\nThe unified naive baseline preserves 100% recall against the '
        'committed ground truth.')


def _project(per_q, threads):
  """Print what the naive sweep costs at each sweep size on THIS machine.

  Sampling is stratified by length, but per-query cost still varies between queries, so
  treat these as an order-of-magnitude guide for choosing SIZES, not a promised wall time.
  """
  print()
  print(f'projected naive Brute Force search time (wall, this machine, {threads} threads):')
  print(f'  per query   {per_q:8.3f} s   ({per_q * threads:.1f} core-s)')
  for n in (100, 1000, 10000, 100000, 1000000):
    hrs = per_q * n / 3600
    note = ''
    if hrs > 168:
      note = '   <-- exceeds the 7-day QOS wall'
    elif hrs > 48:
      note = '   <-- exceeds the 48h sbatch wall (raise #SBATCH --time)'
    print(f'  N={n:>9,}  {hrs:10.2f} h{note}')
  print()
  print('Cumulative wall for a sweep is the SUM of its sizes, and the largest size')
  print('dominates. Divide by (cluster threads / these threads) for a cluster estimate.')


def main():
  p = argparse.ArgumentParser()
  p.add_argument('--indels', type=int, choices=[1, 2], required=True)
  p.add_argument('--queries', default=None,
                 help='FASTA of queries. Default: queries/synth-<n>indel-100.fasta')
  p.add_argument('--proteome', default=str(BENCH / 'proteomes' / 'human.fasta'))
  p.add_argument('--limit', type=int, default=20,
                 help='Sample this many queries (naive BF is O(proteome) per query).')
  p.add_argument('--seed', type=int, default=20260831)
  p.add_argument('--threads', type=int, default=4)
  p.add_argument('--skip-prefiltered', action='store_true',
                 help='Only time the naive sweep (for large --limit cost probes).')
  p.add_argument('--check-oracle-equivalence', action='store_true',
                 help='For --indels 1 only: prove the generic aligner (max_indels=1) that '
                      'the naive baseline uses returns the SAME rows as the committed '
                      '1-indel oracle the ground truth was built with. Runs against each '
                      "query's real matches plus a random protein sample, then exits.")
  args = p.parse_args()

  if args.check_oracle_equivalence:
    _check_oracle_equivalence(args)
    return

  qpath = Path(args.queries) if args.queries else (
    BENCH / 'queries' / f'synth-{args.indels}indel-100.fasta')
  if not qpath.is_absolute():
    qpath = BENCH / qpath
  if not qpath.exists():
    raise SystemExit(f'FATAL: query file not found: {qpath}')
  if not Path(args.proteome).exists():
    raise SystemExit(f'FATAL: proteome not found: {args.proteome}')

  seqs = [ln.strip() for ln in open(qpath) if ln.strip() and not ln.startswith('>')]
  seqs = sorted(set(seqs))
  if args.limit and args.limit < len(seqs):
    # Stratify by length: per-query cost spans orders of magnitude with length, so a
    # uniform sample would be dominated by the cheap long queries and under-report cost.
    by_len = {}
    for q in seqs:
      by_len.setdefault(len(q), []).append(q)
    rng = random.Random(args.seed)
    pick, lengths = [], sorted(by_len)
    while len(pick) < args.limit and any(by_len[L] for L in lengths):
      for L in lengths:
        if by_len[L] and len(pick) < args.limit:
          pick.append(by_len[L].pop(rng.randrange(len(by_len[L]))))
    seqs = sorted(pick)

  with tempfile.NamedTemporaryFile('w', suffix='.fasta', delete=False) as f:
    for i, q in enumerate(seqs, start=1):
      f.write(f'>{i}\n{q}\n')
    sample_path = f.name

  lengths = sorted({len(q) for q in seqs})
  print(f'queries   : {len(seqs)} sampled from {qpath.name} (lengths {lengths[0]}-{lengths[-1]})')
  print(f'proteome  : {args.proteome}')
  print(f'indels    : {args.indels}   threads: {args.threads}')
  print()

  kwargs = dict(
    benchmark='equivalence-check', query=sample_path, proteome=args.proteome,
    lengths=[], max_mismatches=0, method_parameters={'threads': args.threads},
    indels=args.indels,
  )

  naive = _load_method('brute_force_naive')
  t0 = time.perf_counter()
  naive_df = naive.Benchmarker(**kwargs).search()
  t_naive = time.perf_counter() - t0
  print(f'naive (no prefilter)   : {t_naive:8.2f}s   {len(naive_df):>8,} rows')

  if args.skip_prefiltered:
    _project(t_naive / len(seqs), args.threads)
    return

  pre = _load_method('brute_force')
  t0 = time.perf_counter()
  pre_df = pre.Benchmarker(**kwargs).search()
  t_pre = time.perf_counter() - t0
  print(f'prefiltered (current)  : {t_pre:8.2f}s   {len(pre_df):>8,} rows')

  a, b = _rows(naive_df), _rows(pre_df)
  print()
  if a == b:
    ratio = (t_naive / t_pre) if t_pre else float('nan')
    print(f'EQUIVALENT: both methods returned the same {len(a):,} rows.')
    print(f'naive is {ratio:.1f}x slower than the prefiltered baseline on this sample.')
  else:
    only_naive, only_pre = a - b, b - a
    print('*** NOT EQUIVALENT -- do not publish numbers from either method. ***')
    print(f'  rows only in naive       : {len(only_naive):,}')
    print(f'  rows only in prefiltered : {len(only_pre):,}')
    for label, rows in (('naive', only_naive), ('prefiltered', only_pre)):
      for r in list(sorted(rows))[:5]:
        print(f'    only in {label}: {r}')
    # A row found ONLY by the naive sweep means the prefilter dropped a real match, which
    # would also mean the committed ground truth is incomplete. Loudest possible failure.
    if only_naive:
      print('\n  The naive sweep found matches the prefilter MISSED: the prefilter is')
      print('  unsound and the synthetic ground truth built with it is incomplete.')
    raise SystemExit(1)

  _project(t_naive / len(seqs), args.threads)


if __name__ == '__main__':
  main()
