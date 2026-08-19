import argparse
import importlib
import inspect
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

METHODS_DIR = str(Path(__file__).parent / 'methods')
if METHODS_DIR not in sys.path:
  sys.path.insert(0, METHODS_DIR)

RESULT_COLUMNS = ['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start']


def load_config():
  with open(Path(__file__).parent / 'benchmarking_parameters.json') as f:
    return json.load(f)


def load_method(name, benchmark, dataset, method_params):
  try:
    if name == 'PEPMatch':
      module = importlib.import_module('pepmatch.benchmarker')
    else:
      module = importlib.import_module(name)

    kwargs = dict(
      benchmark=benchmark,
      query=Path(__file__).parent / dataset['query'],
      proteome=Path(__file__).parent / dataset['proteome'],
      lengths=dataset['lengths'],
      max_mismatches=dataset['mismatches'],
      method_parameters=method_params,
    )
    # Pass `indels` only to wrappers that declare it, so the existing mismatch
    # datasets and their wrappers are completely unaffected.
    if 'indels' in inspect.signature(module.Benchmarker.__init__).parameters:
      kwargs['indels'] = dataset.get('indels', 0)
    return module.Benchmarker(**kwargs)
  except ValueError as e:
    print(f'  Skipping {name}: {e}')
    return None


def pin_threads(threads):
  """Pin every tool to the same thread count for a fair timing comparison.

  PEPMatch's Rust engine parallelises over the global rayon pool, which reads
  RAYON_NUM_THREADS at first use; the aligners and brute force take a per-method
  `threads` parameter (see each wrapper). Setting the env var here covers PEPMatch
  in-process and is inherited by the memory-measurement subprocesses.
  """
  os.environ['RAYON_NUM_THREADS'] = str(threads)
  os.environ['OMP_NUM_THREADS'] = str(threads)


def with_threads(method_params, threads):
  """Inject the pinned thread count into a method's parameters without mutating
  the shared config dict."""
  return dict(method_params, threads=threads)


def time_step(fn):
  try:
    start = time.perf_counter()
    result = fn()
    return time.perf_counter() - start, result
  except TypeError:
    return None, None


def recall(results_df, expected_df):
  results = results_df[RESULT_COLUMNS].drop_duplicates()
  expected = expected_df[RESULT_COLUMNS].drop_duplicates()

  results['Index start'] = results['Index start'].fillna(0).astype(int)
  expected['Index start'] = expected['Index start'].fillna(0).astype(int)

  matched = pd.merge(results, expected, how='inner', on=RESULT_COLUMNS).drop_duplicates()
  return min((len(matched) / len(expected)) * 100, 100)


def measure_peak_rss_mb(benchmark, name, threads):
  """Peak resident memory (MB) for one method, measured in a FRESH subprocess.

  A fresh process gives a clean per-method high-water mark. We take the max of the
  process's own peak RSS (in-process tools -- PEPMatch's Rust engine, brute force)
  and its largest child's peak RSS (subprocess tools -- BLAST/DIAMOND/MMseqs2).
  This is the fair cross-tool measure the old tracemalloc path could not give: it
  only saw Python heap allocations, so it undercounted the Rust engine and missed
  the aligner subprocesses entirely. Runs in a subprocess so a crashing tool cannot
  take down the driver.
  """
  proc = subprocess.run(
    [sys.executable, __file__, '-b', benchmark,
     '--mem-method', name, '--threads', str(threads)],
    capture_output=True, text=True, env=os.environ,
  )
  for line in proc.stdout.splitlines():
    if line.startswith('PEAK_RSS_KB='):
      return int(line.split('=', 1)[1]) / 1024.0  # ru_maxrss is KB on Linux
  return None


def run_single_method_memory(benchmark, name, threads):
  """Internal entry point (--mem-method): run exactly one method end to end, then
  print its peak RSS as `PEAK_RSS_KB=<n>` for measure_peak_rss_mb to parse."""
  config = load_config()
  dataset = config['datasets'][benchmark]
  method = next(m for m in config['methods'] if m['name'] == name)

  pin_threads(threads)
  tool = load_method(name, benchmark, dataset, with_threads(method['method_parameters'], threads))
  if tool is None:
    print('PEAK_RSS_KB=0')
    return

  for step in (tool.preprocess_proteome, tool.preprocess_query, tool.search):
    try:
      step()
    except TypeError:
      pass  # tool legitimately skips this step (e.g. no query preprocessing)
  if hasattr(tool, 'cleanup'):
    tool.cleanup()

  peak = max(
    resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
  )
  print(f'PEAK_RSS_KB={peak}')


def run_benchmark(benchmark, include_memory=False, include_text_shifting=False, threads=1):
  config = load_config()
  dataset = config['datasets'][benchmark]
  methods = config['methods']

  if not include_text_shifting:
    methods = [m for m in methods if not m['text_shifting']]

  pin_threads(threads)

  expected_df = pd.read_csv(
    Path(__file__).parent / dataset['expected'], sep='\t'
  )

  rows = []

  for method in methods:
    name = method['name']
    print(f'\n{"=" * 60}')
    print(f'  {name}')
    print(f'{"=" * 60}')

    # One failing method must never take down the whole run: the real timed run
    # happens on the cluster where the driver cannot be debugged live, so a tool
    # that crashes (missing binary, segfault, bad output) is skipped, not fatal.
    try:
      tool = load_method(name, benchmark, dataset, with_threads(method['method_parameters'], threads))
      if tool is None:
        continue

      # preprocess proteome
      print('  Preprocessing proteome...')
      preprocess_proteome_time, _ = time_step(tool.preprocess_proteome)
      if preprocess_proteome_time is not None:
        print(f'  -> {preprocess_proteome_time:.3f}s')
      else:
        print('  -> N/A')

      # preprocess query
      print('  Preprocessing query...')
      preprocess_query_time, _ = time_step(tool.preprocess_query)
      if preprocess_query_time is not None:
        print(f'  -> {preprocess_query_time:.3f}s')
      else:
        print('  -> N/A')

      # search
      print('  Searching...')
      search_time, results_df = time_step(tool.search)
      print(f'  -> {search_time:.3f}s')

      # total
      total = sum(t for t in [preprocess_proteome_time, preprocess_query_time, search_time] if t is not None)

      # memory (fresh subprocess; fair across Python/Rust/subprocess tools)
      memory = None
      if include_memory:
        print('  Measuring memory...')
        memory = measure_peak_rss_mb(benchmark, name, threads)
        if memory is not None:
          print(f'  -> {memory:.1f} MB')
        else:
          print('  -> N/A')

      # recall
      recall_pct = recall(results_df, expected_df)
      print(f'  Recall: {recall_pct:.1f}%')

      # cleanup
      if hasattr(tool, 'cleanup'):
        tool.cleanup()

      rows.append({
        'Method': str(tool),
        'Proteome Preprocessing (s)': f'{preprocess_proteome_time:.3f}' if preprocess_proteome_time is not None else 'N/A',
        'Query Preprocessing (s)': f'{preprocess_query_time:.3f}' if preprocess_query_time is not None else 'N/A',
        'Searching (s)': f'{search_time:.3f}',
        'Total (s)': f'{total:.3f}',
        'Memory (MB)': f'{memory:.1f}' if memory is not None else 'N/A',
        'Recall (%)': f'{recall_pct:.1f}',
      })
    except Exception as e:  # noqa: BLE001 -- deliberately broad; keep the run alive
      print(f'  !! {name} failed, skipping: {type(e).__name__}: {e}')
      continue

  results = pd.DataFrame(rows)

  print(f'\n\n{"=" * 80}')
  print(f'  {benchmark.upper()} RESULTS')
  print(f'{"=" * 80}\n')
  print(results.to_string(index=False))

  output_path = f'{benchmark}_benchmarking.tsv'
  results.to_csv(output_path, sep='\t', index=False)
  print(f'\nSaved to {output_path}')

  return results


def main():
  parser = argparse.ArgumentParser(description='PEPMatch Benchmarking Framework')
  parser.add_argument(
    '-b', '--benchmark',
    choices=['mhc_ligands', 'milk', 'coronavirus', 'neoepitopes', 'cosmic_indel', 'cedar_indel'],
    required=True,
  )
  parser.add_argument('-m', '--memory', action='store_true', default=False)
  parser.add_argument('-t', '--text_shifting', action='store_true', default=False)
  parser.add_argument(
    '-p', '--threads', type=int, default=1,
    help='Threads every tool is pinned to, for a fair timing comparison (default 1). '
         'On the cluster, set this to match --cpus-per-task.',
  )
  # Internal: measure one method's peak RSS in a fresh process (see measure_peak_rss_mb).
  parser.add_argument('--mem-method', default=None, help=argparse.SUPPRESS)
  args = parser.parse_args()

  if args.mem_method:
    run_single_method_memory(args.benchmark, args.mem_method, args.threads)
    return

  run_benchmark(args.benchmark, args.memory, args.text_shifting, args.threads)


if __name__ == '__main__':
  main()
