#!/usr/bin/env python3
"""Phase 1 validation of the indel-annotation speedup: old vs new code, full scale.

WHAT IT PROVES
  The speedup (annotate each distinct (query, matched) pair once; find placements from
  string agreement instead of every combination) must change NOTHING but time. For each
  dataset behind Figures 3 and 4 this runs the original matcher.py ("old") and the new
  one ("new") against the SAME compiled engine and requires:
    * the full output table, Indel Positions column included, byte-identical;
    * recall 100.0 against the generator's ground truth, for both;
  and reports the search-time speedup. Any failure exits non-zero.

HOW THE TWO VERSIONS RUN
  Only pepmatch/matcher.py differs, so the job builds the engine once and makes a copy
  of the package with the original matcher.py swapped in; "old" runs import that copy via
  PYTHONPATH. Each (dataset, version) is its own process, because the Rust thread pool is
  sized once per process and Figures 3 and 4 used different thread counts (32 and 64).

TIMING = THE FIGURES' DEFINITION
  pepmatch.benchmarker.Benchmarker.search(), exactly what benchmarking.py timed for the
  PEPMatch "Searching (s)" column, after preprocess_proteome() (index build, excluded).
  Repeated --repeats times; the first repeat is the published protocol (no warm-up), the
  median is the robust number.

  validate-annotation-speedup.py run --version old|new --dataset cosmic_indel --threads 32 --out DIR
  validate-annotation-speedup.py compare --out DIR
"""
import argparse
import csv
import hashlib
import inspect
import os
import platform
import statistics
import sys
import time
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
RUN_FIELDS = ['dataset', 'version', 'threads', 'repeat', 'search_s', 'preprocess_s', 'rows',
              'hit_rows', 'recall', 'output_sha256', 'pepmatch_file', 'host']


def sha256_file(path, buf=1 << 20):
  h = hashlib.sha256()
  with open(path, 'rb') as f:
    while chunk := f.read(buf):
      h.update(chunk)
  return h.hexdigest()


def run(args):
  # Rayon sizes its pool at first use, so this must precede anything touching pepmatch.
  os.environ['RAYON_NUM_THREADS'] = str(args.threads)
  os.environ['OMP_NUM_THREADS'] = str(args.threads)
  import pandas as pd
  import pepmatch
  from pepmatch import Matcher
  from pepmatch import matcher as M
  sys.path.insert(0, str(BENCH))
  import benchmarking as H

  # Refuse to report a number for the wrong code: the copy on PYTHONPATH must be the
  # one actually imported, and it must be the version we were told it is.
  is_new = 'runs[s][x]' in inspect.getsource(M._indel_placements)
  if is_new != (args.version == 'new'):
    sys.exit(f'FATAL: asked for {args.version!r} but imported '
             f'{"new" if is_new else "old"} matcher.py from {M.__file__}')

  config = H.load_config()
  dataset = config['datasets'][args.dataset]
  expected = pd.read_csv(BENCH / dataset['expected'], sep='\t')
  tool = H.load_method('PEPMatch', args.dataset, dataset, {'threads': args.threads})

  t = time.perf_counter()
  tool.preprocess_proteome()
  pre = time.perf_counter() - t

  out = Path(args.out)
  out.mkdir(parents=True, exist_ok=True)
  rows_out = []
  for rep in range(1, args.repeats + 1):
    t = time.perf_counter()
    found = tool.search()
    secs = time.perf_counter() - t
    rec = H.recall(found, expected)
    rows_out.append((rep, secs, rec))
    print(f'  {args.dataset:18s} {args.version} rep {rep}: search {secs:9.3f} s  '
          f'recall {rec:.1f}', flush=True)

  # The full table, Indel Positions included: the search() frame above keeps only the
  # four recall columns, and the annotation is the part this change touched.
  full = Matcher(query=tool.query, proteome_file=tool.proteome, max_indels=tool.indels,
                 output_format='dataframe', sequence_version=False).match()
  tables = Path(args.tables) if args.tables else out / 'tables'
  table = tables / f'{args.dataset}.{args.version}.tsv'
  table.parent.mkdir(parents=True, exist_ok=True)
  full.write_csv(table, separator='\t')
  digest = sha256_file(table)
  hits = int(full['Matched Sequence'].is_not_null().sum())

  runs_path = out / 'validation-runs.tsv'
  new_file = not runs_path.exists()
  with open(runs_path, 'a', newline='') as f:
    w = csv.DictWriter(f, fieldnames=RUN_FIELDS, delimiter='\t')
    if new_file:
      w.writeheader()
    for rep, secs, rec in rows_out:
      w.writerow(dict(dataset=args.dataset, version=args.version, threads=args.threads,
                      repeat=rep, search_s=f'{secs:.4f}', preprocess_s=f'{pre:.3f}',
                      rows=full.height, hit_rows=hits, recall=f'{rec:.1f}',
                      output_sha256=digest, pepmatch_file=pepmatch.__file__,
                      host=platform.node()))
  print(f'  {args.dataset:18s} {args.version} table {full.height:,} rows  sha {digest[:16]}',
        flush=True)


def compare(args):
  out = Path(args.out)
  rows = list(csv.DictReader(open(out / 'validation-runs.tsv'), delimiter='\t'))
  by = {}
  for r in rows:
    by.setdefault((r['dataset'], r['version']), []).append(r)
  datasets = list(dict.fromkeys(r['dataset'] for r in rows))
  failures = []
  summary = []
  print(f'\n{"dataset":18s} {"thr":>3s} {"rows":>9s}  {"identical":>9s}  {"recall":>11s}  '
        f'{"old s (1st/med)":>17s}  {"new s (1st/med)":>17s}  {"speedup":>7s}')
  for ds in datasets:
    old, new = by.get((ds, 'old')), by.get((ds, 'new'))
    if not old or not new:
      failures.append(f'{ds}: missing the {"old" if not old else "new"} run')
      continue
    same = old[0]['output_sha256'] == new[0]['output_sha256'] and old[0]['rows'] == new[0]['rows']
    recalls = {float(r['recall']) for r in old + new}
    o = [float(r['search_s']) for r in old]
    n = [float(r['search_s']) for r in new]
    speed = statistics.median(o) / statistics.median(n)
    if not same:
      failures.append(f'{ds}: OUTPUT DIFFERS (old {old[0]["output_sha256"][:12]} '
                      f'vs new {new[0]["output_sha256"][:12]})')
    if recalls != {100.0}:
      failures.append(f'{ds}: recall {sorted(recalls)}, not 100.0')
    print(f'{ds:18s} {old[0]["threads"]:>3s} {int(old[0]["rows"]):>9,d}  '
          f'{"YES" if same else "NO":>9s}  {"/".join(f"{x:.1f}" for x in sorted(recalls)):>11s}  '
          f'{o[0]:8.3f}/{statistics.median(o):8.3f}  {n[0]:8.3f}/{statistics.median(n):8.3f}  '
          f'{speed:6.1f}x')
    summary.append(dict(dataset=ds, threads=old[0]['threads'], rows=old[0]['rows'],
                        hit_rows=old[0]['hit_rows'], identical='yes' if same else 'NO',
                        recall='/'.join(f'{x:.1f}' for x in sorted(recalls)),
                        old_first_s=f'{o[0]:.4f}', old_median_s=f'{statistics.median(o):.4f}',
                        new_first_s=f'{n[0]:.4f}', new_median_s=f'{statistics.median(n):.4f}',
                        speedup=f'{speed:.2f}', output_sha256=old[0]['output_sha256']))
  with open(out / 'validation-summary.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(summary[0]) if summary else ['dataset'],
                       delimiter='\t')
    w.writeheader()
    w.writerows(summary)
  if failures:
    print('\nVALIDATION FAILED:')
    for x in failures:
      print(f'  - {x}')
    sys.exit(1)
  print(f'\nVALIDATION PASSED: {len(summary)} datasets, every output byte-identical, '
        f'recall 100.0 for old and new.')


def main():
  p = argparse.ArgumentParser(description=__doc__,
                              formatter_class=argparse.RawDescriptionHelpFormatter)
  sub = p.add_subparsers(dest='cmd', required=True)
  r = sub.add_parser('run')
  r.add_argument('--version', choices=['old', 'new'], required=True)
  r.add_argument('--dataset', required=True)
  r.add_argument('--threads', type=int, required=True)
  r.add_argument('--repeats', type=int, default=3)
  r.add_argument('--out', required=True)
  r.add_argument('--tables', default=None,
                 help='where the full output tables go (default OUT/tables); large at 100k')
  c = sub.add_parser('compare')
  c.add_argument('--out', required=True)
  args = p.parse_args()
  run(args) if args.cmd == 'run' else compare(args)


if __name__ == '__main__':
  main()
