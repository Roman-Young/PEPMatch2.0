#!/usr/bin/env python3
"""Repeat PEPMatch's search on benchmark datasets and report the median, so a figure does not
hinge on one run on a shared node.

Uses the same Benchmarker the timing tables use (pepmatch/benchmarker.py): the index is built
once per dataset (not timed, as in the paper), then Benchmarker.search() is timed REPEATS
times. Recall is checked on the first repeat with benchmarking.py's own recall().

  repeat-pepmatch-timing.py --datasets cosmic_indel,cedar_indel --repeats 5 --out results/x.tsv
"""
import argparse
import json
import statistics
import sys
import time
from pathlib import Path

import pandas as pd

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH))
from benchmarking import recall, resolve_proteome   # noqa: E402  the exact functions the tables use


def main():
  p = argparse.ArgumentParser()
  p.add_argument('--datasets', required=True)
  p.add_argument('--repeats', type=int, default=5)
  p.add_argument('--out', required=True)
  args = p.parse_args()

  from pepmatch.benchmarker import Benchmarker
  config = json.load(open(BENCH / 'benchmarking_parameters.json'))['datasets']
  rows = []
  for name in args.datasets.split(','):
    ds = config[name]
    if not (BENCH / ds['query']).exists():
      print(f'SKIP {name}: {ds["query"]} not generated here', flush=True)
      continue
    b = Benchmarker(name, BENCH / ds['query'], resolve_proteome(ds), ds['lengths'], ds['mismatches'], {},
                    indels=ds.get('indels', 0))
    t = time.perf_counter()
    b.preprocess_proteome()
    build = time.perf_counter() - t
    times, rec = [], None
    for i in range(args.repeats):
      t = time.perf_counter()
      df = b.search()
      times.append(time.perf_counter() - t)
      if i == 0:
        rec = recall(df, pd.read_csv(BENCH / ds['expected'], sep='\t'))
    row = {'Dataset': name, 'Repeats': args.repeats, 'Search median (s)': round(statistics.median(times), 3),
           'Search min (s)': round(min(times), 3), 'Search max (s)': round(max(times), 3),
           'All runs (s)': ' '.join(f'{x:.3f}' for x in times), 'Index build (s)': round(build, 3),
           'Recall (%)': round(rec, 1)}
    rows.append(row)
    print('  '.join(f'{k}={v}' for k, v in row.items()), flush=True)
  Path(args.out).parent.mkdir(parents=True, exist_ok=True)
  pd.DataFrame(rows).to_csv(args.out, sep='\t', index=False)
  print(f'wrote {args.out}')


if __name__ == '__main__':
  main()
