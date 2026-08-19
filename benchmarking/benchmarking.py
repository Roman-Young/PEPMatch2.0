import argparse
import importlib
import inspect
import json
import os
import resource
import subprocess
import sys
import time
import traceback
from pathlib import Path

import pandas as pd

METHODS_DIR = str(Path(__file__).parent / 'methods')
if METHODS_DIR not in sys.path:
  sys.path.insert(0, METHODS_DIR)

RESULT_COLUMNS = ['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start']

RESULT_TABLE_COLUMNS = [
  'Method', 'Proteome Preprocessing (s)', 'Query Preprocessing (s)',
  'Searching (s)', 'Total (s)', 'Memory (MB)', 'Recall (%)', 'Status',
]


class StepNotApplicable(Exception):
  """Raised by a wrapper for a step it legitimately does not perform (e.g. tools
  that do not preprocess queries).

  This exists because "step not applicable" and "step crashed" used to be the same
  signal: time_step caught bare TypeError, so a genuine TypeError bug inside search()
  was swallowed, the method was dropped from the table, and the only trace was a log
  line that looked exactly like the normal N/A prints. A missing row in a published
  table must never be that easy to produce.
  """


def load_config():
  with open(Path(__file__).parent / 'benchmarking_parameters.json') as f:
    return json.load(f)


def resolve_proteome(dataset):
  """Where to read the proteome from.

  BLAST/DIAMOND/MMseqs2 all write their databases NEXT TO the proteome file, and an
  MMseqs2 index over the human proteome runs to multiple GB. On a cluster that would
  land on the home quota. Setting PEPMATCH_BENCH_PROTEOME_DIR to a scratch copy moves
  every generated database there, leaving only the (tiny) results tables in home.
  """
  override = os.environ.get('PEPMATCH_BENCH_PROTEOME_DIR')
  if override:
    return Path(override) / Path(dataset['proteome']).name
  return Path(__file__).parent / dataset['proteome']


def load_method(name, benchmark, dataset, method_params):
  try:
    if name == 'PEPMatch':
      module = importlib.import_module('pepmatch.benchmarker')
    else:
      module = importlib.import_module(name)

    kwargs = dict(
      benchmark=benchmark,
      query=Path(__file__).parent / dataset['query'],
      proteome=resolve_proteome(dataset),
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


def failed_row(name, status, detail):
  """A table row recording that a method did not produce results, and why."""
  row = {c: 'N/A' for c in RESULT_TABLE_COLUMNS}
  row['Method'] = name
  row['Status'] = f'{status}: {detail}'
  return row


def append_row(rows, output_path, row):
  """Append one method's row to the in-memory table AND flush it to disk immediately.

  Writing only at the end meant a wall-clock kill produced no output at all; now a
  truncated run still yields a partial table of everything that finished.
  """
  rows.append(row)
  write_header = not Path(output_path).exists()
  pd.DataFrame([row], columns=RESULT_TABLE_COLUMNS).to_csv(
    output_path, sep='\t', index=False, mode='a', header=write_header
  )


def time_optional_step(fn):
  """Time a step a tool may legitimately not implement (the preprocessing steps).

  StepNotApplicable is the intended signal; bare TypeError is still accepted because
  the existing wrappers raise it for 'does not preprocess queries'.
  """
  try:
    start = time.perf_counter()
    result = fn()
    return time.perf_counter() - start, result
  except (StepNotApplicable, TypeError):
    return None, None


def time_required_step(fn):
  """Time a step every tool must implement (search).

  Deliberately does NOT swallow exceptions: a crash here is a real failure and must
  reach the per-method handler so it is recorded as FAILED rather than vanishing.
  """
  start = time.perf_counter()
  result = fn()
  return time.perf_counter() - start, result


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
  output_path = str(Path(__file__).parent / f'{benchmark}_benchmarking.tsv')
  # Start clean, then append after every method. A SLURM wall-clock kill mid-run used
  # to lose everything, including methods that had already finished hours earlier.
  Path(output_path).unlink(missing_ok=True)

  for method in methods:
    name = method['name']
    # A method may appear more than once with different parameters (e.g. BLAST is run
    # in both blastp-short and default blastp modes), so `label` is what identifies the
    # row; it falls back to the module name for single-configuration methods.
    label = method.get('label', name)
    print(f'\n{"=" * 60}')
    print(f'  {label}')
    print(f'{"=" * 60}')
    sys.stdout.flush()

    # One failing method must never take down the whole run: the real timed run
    # happens on the cluster where the driver cannot be debugged live, so a tool
    # that crashes (missing binary, segfault, bad output) is recorded as FAILED
    # and the run continues.
    try:
      tool = load_method(name, benchmark, dataset, with_threads(method['method_parameters'], threads))
      if tool is None:
        append_row(rows, output_path, failed_row(label, 'SKIPPED', 'not applicable to this dataset'))
        continue

      # preprocess proteome
      print('  Preprocessing proteome...')
      preprocess_proteome_time, _ = time_optional_step(tool.preprocess_proteome)
      if preprocess_proteome_time is not None:
        print(f'  -> {preprocess_proteome_time:.3f}s')
      else:
        print('  -> N/A')

      # preprocess query
      print('  Preprocessing query...')
      preprocess_query_time, _ = time_optional_step(tool.preprocess_query)
      if preprocess_query_time is not None:
        print(f'  -> {preprocess_query_time:.3f}s')
      else:
        print('  -> N/A')

      # search -- required; a crash here propagates and is recorded as FAILED
      print('  Searching...')
      search_time, results_df = time_required_step(tool.search)
      print(f'  -> {search_time:.3f}s')

      # total
      total = sum(t for t in [preprocess_proteome_time, preprocess_query_time, search_time] if t is not None)

      # memory (fresh subprocess; off by default -- see measure_peak_rss_mb)
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

      append_row(rows, output_path, {
        'Method': str(tool),
        'Proteome Preprocessing (s)': f'{preprocess_proteome_time:.3f}' if preprocess_proteome_time is not None else 'N/A',
        'Query Preprocessing (s)': f'{preprocess_query_time:.3f}' if preprocess_query_time is not None else 'N/A',
        'Searching (s)': f'{search_time:.3f}',
        'Total (s)': f'{total:.3f}',
        'Memory (MB)': f'{memory:.1f}' if memory is not None else 'N/A',
        'Recall (%)': f'{recall_pct:.1f}',
        'Status': 'OK',
      })
    except Exception as e:  # noqa: BLE001 -- deliberately broad; keep the run alive
      detail = f'{type(e).__name__}: {e}'
      print(f'  !! {label} FAILED: {detail}')
      traceback.print_exc()
      # Record the failure IN THE TABLE. A silently missing row is indistinguishable
      # from a method that was never run, which is exactly how a bad table gets published.
      append_row(rows, output_path, failed_row(label, 'FAILED', detail))
      continue

  results = pd.DataFrame(rows, columns=RESULT_TABLE_COLUMNS)

  print(f'\n\n{"=" * 80}')
  print(f'  {benchmark.upper()} RESULTS')
  print(f'{"=" * 80}\n')
  print(results.to_string(index=False))

  # The table is already on disk (written incrementally); rewrite it once so the
  # final file is clean, then state plainly whether every method actually reported.
  results.to_csv(output_path, sep='\t', index=False)
  print(f'\nSaved to {output_path}')

  failures = [r['Method'] for r in rows if r.get('Status', 'OK') != 'OK']
  if len(rows) != len(methods):
    print(f'\n!! INCOMPLETE: {len(rows)} rows for {len(methods)} methods -- rows were lost.')
  if failures:
    print(f'!! {len(failures)} method(s) did not produce results: {", ".join(failures)}')
  else:
    print(f'All {len(rows)} methods reported successfully.')

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
