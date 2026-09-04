#!/usr/bin/env python3
"""Large-scale recall validation of the reduced n+1-seed indel search (PEPMATCH_SEED_STRATEGY
=first) against the FULL-coverage engine and the committed proteome-wide oracle.

Runs the real Rust engine twice on the same query set -- strategy=all (current shipped
behavior) and strategy=first (n+1 disjoint seeds) -- and checks three things:

  1. first == all           no match is lost, none is spuriously added by the reduction
  2. oracle subset of first every committed ground-truth match is still found
  3. first subset of oracle  the reduction introduces no match the oracle doesn't have

Because it runs against the synthetic sets' PROTEOME-WIDE ground truth (every placement in
every protein, from the committed brute-force oracle, spot-checked by an exhaustive
no-prefilter scan at generation time), a clean result at 100k is a genuine recall guarantee,
not a sample.

Usage (from the repo root, in the pepmatch venv):
  python benchmarking/scripts/validate-seed-reduction.py --prefix synth-1indel --size 100k --indels 1
  python benchmarking/scripts/validate-seed-reduction.py --prefix synth-2indel --size 100k --indels 2

Paths resolve relative to this script, so it runs the same on the dev box and the cluster.
"""
import argparse
import csv
import subprocess
import sys
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
PROTEOME = BENCH / 'proteomes' / 'human.fasta'


def run_engine(qfasta, n, strategy, out_path):
  """Run one seed strategy in a subprocess (env var read fresh by the Rust each time),
  writing indel matches to out_path. The Matcher prints progress to stdout, so results go
  to a file, never the stream."""
  code = f'''
import os
os.environ["PEPMATCH_SEED_STRATEGY"] = {strategy!r}
from pepmatch import Matcher
df = Matcher(query={str(qfasta)!r}, proteome_file={str(PROTEOME)!r}, max_indels={n},
             output_format="dataframe", sequence_version=False).match()
df = df.select(["Query Sequence","Matched Sequence","Protein ID","Index start"])
import csv
with open({str(out_path)!r}, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\\t")
    for q, m, p, s in df.iter_rows():
        if m is None or len(str(m)) == len(str(q)):   # indel matches only, mirror oracle TSVs
            continue
        w.writerow([q, m, p, s])
'''
  r = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, cwd=str(BENCH.parent))
  if r.returncode != 0:
    sys.stderr.write(r.stderr[-3000:])
    raise SystemExit(f'engine run failed (strategy={strategy})')


def load(path, cols4=True):
  rows = set()
  with open(path, newline='') as fh:
    if cols4:
      for row in csv.reader(fh, delimiter='\t'):
        if row:
          rows.add((row[0], row[1], row[2], str(row[3])))
    else:
      for r in csv.DictReader(fh, delimiter='\t'):
        rows.add((r['Query Sequence'], r['Matched Sequence'], r['Protein ID'], str(r['Index start'])))
  return rows


def main():
  ap = argparse.ArgumentParser()
  ap.add_argument('--prefix', required=True)
  ap.add_argument('--size', required=True)
  ap.add_argument('--indels', type=int, required=True)
  args = ap.parse_args()

  qfasta = BENCH / 'queries' / f'{args.prefix}-{args.size}.fasta'
  expected = BENCH / 'expected' / f'{args.prefix}-{args.size}-expected.tsv'
  for p in (qfasta, expected, PROTEOME):
    if not p.exists():
      raise SystemExit(f'FATAL: missing {p}')

  tmp = BENCH / 'results' / '_seedcheck'
  tmp.mkdir(parents=True, exist_ok=True)
  all_out, first_out = tmp / 'all.tsv', tmp / 'first.tsv'

  print(f'query set : {qfasta.name}   (n={args.indels})')
  print('running strategy=all  ...', flush=True)
  run_engine(qfasta, args.indels, 'all', all_out)
  print('running strategy=first...', flush=True)
  run_engine(qfasta, args.indels, 'first', first_out)

  all_rows, first_rows = load(all_out), load(first_out)
  oracle = load(expected, cols4=False)

  lost = all_rows - first_rows          # in full coverage, dropped by reduction
  spurious = first_rows - all_rows      # produced only by reduction
  oracle_missed = oracle - first_rows   # ground-truth match reduction fails to find
  oracle_extra = first_rows - oracle    # reduction finds something not in ground truth

  print(f'\n  all    : {len(all_rows):,} matches')
  print(f'  first  : {len(first_rows):,} matches')
  print(f'  oracle : {len(oracle):,} matches')
  print(f'\n  lost by reduction (all\\first)      : {len(lost):,}')
  print(f'  spurious from reduction (first\\all): {len(spurious):,}')
  print(f'  oracle matches missed by first     : {len(oracle_missed):,}')
  print(f'  first matches not in oracle        : {len(oracle_extra):,}')

  ok = not (lost or spurious or oracle_missed or oracle_extra)
  for tag, s in (('lost', lost), ('spurious', spurious),
                 ('oracle_missed', oracle_missed), ('oracle_extra', oracle_extra)):
    for row in list(s)[:5]:
      print(f'    {tag}: {row}')
  print('\n' + ('PASS: reduced n+1-seed search is recall- and result-identical to full '
                'coverage and the oracle.' if ok else 'FAIL: see rows above.'))
  return 0 if ok else 1


if __name__ == '__main__':
  sys.exit(main())
