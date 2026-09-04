#!/usr/bin/env python3
"""Companion to validate-seed-reduction.py: over the committed oracle matches, classify which
of the n+1 'first' seeds survive verbatim in each matched window, and count the matches that
FORCE extension from a non-leftmost / purely-interior anchor.

This is what turns 'the reduction passed' into 'the reduction was actually stress-tested':
a match whose only surviving first-seed is interior could ONLY have been found by extending
from an interior anchor (both directions). A large count here, combined with 0 oracle-misses
in validate-seed-reduction.py, means interior-anchor extension was exercised against
independent ground truth that many times.

Pure string analysis over the expected TSV -- no engine, so it's cheap even at 100k.

Usage:
  python benchmarking/scripts/classify-seed-survivors.py --prefix synth-2indel --size 100k --indels 2
"""
import argparse
import csv
from collections import Counter
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]


def partition_k(lengths, n):
  """Mirror Matcher.indel_search grouping: queries too short for k>=3 use k=2, the rest use
  max(2, min_len // (n+1)). Returns {length: k}."""
  short = [L for L in lengths if L // (n + 1) < 3]
  rest = [L for L in lengths if L // (n + 1) >= 3]
  g = {L: 2 for L in short}
  if rest:
    kk = max(2, min(rest) // (n + 1))
    g.update({L: kk for L in rest})
  return g


def first_seeds(query, k, n):
  return [(query[i * k:(i + 1) * k], i * k) for i in range(n + 1) if (i + 1) * k <= len(query)]


def main():
  ap = argparse.ArgumentParser()
  ap.add_argument('--prefix', required=True)
  ap.add_argument('--size', required=True)
  ap.add_argument('--indels', type=int, required=True)
  args = ap.parse_args()
  n = args.indels

  epath = BENCH / 'expected' / f'{args.prefix}-{args.size}-expected.tsv'
  if not epath.exists():
    raise SystemExit(f'FATAL: missing {epath}')

  lengths, rows = set(), []
  with open(epath) as f:
    for r in csv.DictReader(f, delimiter='\t'):
      rows.append((r['Query Sequence'], r['Matched Sequence']))
      lengths.add(len(r['Query Sequence']))
  k_for = partition_k(sorted(lengths), n)

  cls = Counter()
  for q, matched in rows:
    k = k_for[len(q)]
    seeds = first_seeds(q, k, n)
    if len(seeds) < n + 1:
      cls['too_short_fallback'] += 1
      continue
    survivors = [off for (s, off) in seeds if s in matched]
    left, right = seeds[0][1], seeds[-1][1]
    if not survivors:
      cls['NO_SURVIVOR (recall bug!)'] += 1
    elif survivors == [left]:
      cls['only_leftmost'] += 1
    elif left not in survivors and right not in survivors:
      cls['INTERIOR_ONLY (needs L+R extension)'] += 1
    elif left not in survivors:
      cls['no_leftmost (needs LEFT extension)'] += 1
    else:
      cls['leftmost_plus_others'] += 1

  total = sum(cls.values())
  print(f'{args.prefix}-{args.size}  (n={n}, {total:,} oracle matches)')
  for key, c in cls.most_common():
    print(f'  {c:>10,}  {100*c/total:5.1f}%  {key}')
  forced = (cls['INTERIOR_ONLY (needs L+R extension)'] + cls['no_leftmost (needs LEFT extension)'])
  print(f'\n  matches forcing non-leftmost-anchor extension: {forced:,} ({100*forced/total:.1f}%)')
  if cls['NO_SURVIVOR (recall bug!)']:
    print('  *** NO_SURVIVOR > 0: the n+1 seed set fails the pigeonhole guarantee here -- investigate. ***')


if __name__ == '__main__':
  main()
