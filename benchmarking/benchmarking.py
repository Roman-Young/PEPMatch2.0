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

# The memory table reports the headline figure PLUS every raw component it was derived
# from, so the baseline subtraction can be audited instead of taken on trust.
MEMORY_TABLE_COLUMNS = [
  'Method', 'Memory (MB)', 'Peak RSS self (MB)', 'Peak RSS child (MB)',
  'Python baseline (MB)', 'Import cost (MB)', 'External binary', 'Status',
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


def failed_row(name, status, detail, columns=None):
  """A table row recording that a method did not produce results, and why."""
  columns = columns or RESULT_TABLE_COLUMNS
  row = {c: 'N/A' for c in columns}
  row['Method'] = name
  row['Status'] = f'{status}: {detail}'
  return row


def append_row(rows, output_path, row, columns=None):
  """Append one method's row to the in-memory table AND flush it to disk immediately.

  Writing only at the end meant a wall-clock kill produced no output at all; now a
  truncated run still yields a partial table of everything that finished.
  """
  columns = columns or RESULT_TABLE_COLUMNS
  rows.append(row)
  write_header = not Path(output_path).exists()
  pd.DataFrame([row], columns=columns).to_csv(
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


KB_PER_MB = 1024.0   # ru_maxrss is in kilobytes on Linux


def external_tool_runs():
  """How many external binaries the current process exec'd (see methods/_shell.py)."""
  try:
    import _shell
    return _shell.EXTERNAL_RUNS
  except Exception:  # noqa: BLE001 -- absence just means "no external tools ran"
    return 0


def measure_memory(benchmark, method_index, threads):
  """Peak memory for ONE method, measured in a fresh subprocess. Returns a dict.

  Takes the method's INDEX rather than its name: a method may be registered more than
  once with different parameters (BLAST runs in both blastp-short and default modes),
  and selecting by name would silently measure whichever entry came first.

  Raises RuntimeError if the measurement did not produce a result, so a failure lands
  in the table as FAILED instead of a bare 'N/A' that reads like "not applicable".
  """
  proc = subprocess.run(
    [sys.executable, __file__, '-b', benchmark,
     '--mem-index', str(method_index), '--threads', str(threads)],
    capture_output=True, text=True, env=os.environ,
  )

  payload = None
  for line in proc.stdout.splitlines():
    if line.startswith('MEM_JSON='):
      payload = json.loads(line.split('=', 1)[1])

  if payload is None:
    tail = (proc.stderr or '').strip().splitlines()[-5:]
    raise RuntimeError(
      f'memory subprocess exited {proc.returncode} without reporting a result. '
      f'stderr tail: {" | ".join(tail) if tail else "(empty)"}'
    )
  if payload.get('error'):
    raise RuntimeError(payload['error'])
  return payload


def run_single_method_memory(benchmark, method_index, threads):
  """Internal entry point (--mem-index): run one method end to end in a clean process
  and report its memory as a single `MEM_JSON={...}` line.

  WHY THIS IS NOT JUST max(self, children):

  These tools split into two architectures that a raw peak-RSS number cannot compare.
  PEPMatch and brute force run INSIDE this Python process, so their peak includes the
  interpreter, pandas and polars -- a fixed floor of several hundred MB that has nothing
  to do with the algorithm. BLAST/DIAMOND/MMseqs2 exec standalone C binaries whose peak
  carries no Python at all. Reporting both raw made the two in-process methods land on
  an identical number (both dominated by the shared floor) while the aligners looked
  artificially lean -- a measurement artifact, not a result.

  So the reported figure is peak memory ABOVE the interpreter floor:

    * exec'd external binaries  -> max(child peak, self peak - baseline)
      The child is already baseline-free; the parent term covers wrapper work (the
      aligner wrappers build a full proteome dict in-process to resolve match strings).
      Both terms are floor-free, so taking the max is consistent.

    * in-process / fork-parallel -> max(self peak, child peak) - baseline
      A forked worker INHERITS this process's pages, so its RSS contains the same
      baseline; subtracting once is correct for either term.

  Caveat recorded honestly: for fork-parallel methods (brute force) a worker's private
  growth is not added to the parent's, so its figure is the dominant shared footprint
  rather than the sum across workers. The raw components are all reported alongside so
  the reduction can be audited rather than trusted.
  """
  try:
    config = load_config()
    dataset = config['datasets'][benchmark]
    method = config['methods'][method_index]
    pin_threads(threads)

    # The floor: everything already resident before the tool's module is imported.
    baseline_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

    tool = load_method(
      method['name'], benchmark, dataset,
      with_threads(method['method_parameters'], threads),
    )
    if tool is None:
      print('MEM_JSON=' + json.dumps({'skipped': True}))
      return

    after_import_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    label = str(tool)

    # Same step semantics as the timed run: preprocessing may be inapplicable, but a
    # crash in search is a real failure and must propagate.
    for step in (tool.preprocess_proteome, tool.preprocess_query):
      try:
        step()
      except (StepNotApplicable, TypeError):
        pass
    tool.search()
    if hasattr(tool, 'cleanup'):
      tool.cleanup()

    self_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    child_kb = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    external = external_tool_runs() > 0

    if external:
      reported_kb = max(child_kb, self_kb - baseline_kb)
    else:
      reported_kb = max(self_kb, child_kb) - baseline_kb

    print('MEM_JSON=' + json.dumps({
      'label': label,
      'memory_mb': round(max(reported_kb, 0) / KB_PER_MB, 1),
      'self_peak_mb': round(self_kb / KB_PER_MB, 1),
      'child_peak_mb': round(child_kb / KB_PER_MB, 1),
      'baseline_mb': round(baseline_kb / KB_PER_MB, 1),
      'import_cost_mb': round((after_import_kb - baseline_kb) / KB_PER_MB, 1),
      'external_tool': bool(external),
    }))
  except Exception as e:  # noqa: BLE001 -- report the failure through the protocol
    traceback.print_exc()
    print('MEM_JSON=' + json.dumps({'error': f'{type(e).__name__}: {e}'}))


def run_benchmark(benchmark, include_memory=False, include_text_shifting=False, threads=1):
  config = load_config()
  dataset = config['datasets'][benchmark]
  # Carry each method's index in the config, because the memory subprocess selects by
  # index -- a method can appear twice under one name (BLAST short vs default).
  methods = list(enumerate(config['methods']))

  if not include_text_shifting:
    methods = [(i, m) for i, m in methods if not m['text_shifting']]

  pin_threads(threads)

  expected_df = pd.read_csv(
    Path(__file__).parent / dataset['expected'], sep='\t'
  )

  rows = []
  output_path = str(Path(__file__).parent / f'{benchmark}_benchmarking.tsv')
  # Start clean, then append after every method. A SLURM wall-clock kill mid-run used
  # to lose everything, including methods that had already finished hours earlier.
  Path(output_path).unlink(missing_ok=True)

  for method_index, method in methods:
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

      # memory (fresh subprocess; off by default -- see run_single_method_memory)
      memory = None
      if include_memory:
        print('  Measuring memory...')
        payload = measure_memory(benchmark, method_index, threads)
        memory = None if payload.get('skipped') else payload['memory_mb']
        print(f'  -> {memory:.1f} MB' if memory is not None else '  -> N/A')

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


def run_memory_benchmark(benchmark, include_text_shifting=False, threads=1):
  """Produce the MEMORY table only, without redoing the timing run.

  Each method is measured in its own fresh subprocess, so peaks cannot bleed between
  methods. Alongside the reported figure this writes every raw component (self peak,
  child peak, interpreter baseline, import cost, whether external binaries ran) so the
  number in the paper can be audited rather than taken on faith.
  """
  config = load_config()
  methods = list(enumerate(config['methods']))
  if not include_text_shifting:
    methods = [(i, m) for i, m in methods if not m['text_shifting']]

  pin_threads(threads)

  rows = []
  output_path = str(Path(__file__).parent / f'{benchmark}_memory.tsv')
  Path(output_path).unlink(missing_ok=True)

  for method_index, method in methods:
    label = method.get('label', method['name'])
    print(f'\n{"=" * 60}')
    print(f'  {label}  (memory)')
    print(f'{"=" * 60}')
    sys.stdout.flush()

    try:
      payload = measure_memory(benchmark, method_index, threads)
      if payload.get('skipped'):
        append_row(rows, output_path,
                   failed_row(label, 'SKIPPED', 'not applicable to this dataset',
                              MEMORY_TABLE_COLUMNS),
                   MEMORY_TABLE_COLUMNS)
        print('  -> SKIPPED')
        continue

      print(f'  -> {payload["memory_mb"]:.1f} MB '
            f'(self {payload["self_peak_mb"]:.1f} / child {payload["child_peak_mb"]:.1f} / '
            f'baseline {payload["baseline_mb"]:.1f}, '
            f'external={payload["external_tool"]})')

      append_row(rows, output_path, {
        'Method': payload['label'],
        'Memory (MB)': f'{payload["memory_mb"]:.1f}',
        'Peak RSS self (MB)': f'{payload["self_peak_mb"]:.1f}',
        'Peak RSS child (MB)': f'{payload["child_peak_mb"]:.1f}',
        'Python baseline (MB)': f'{payload["baseline_mb"]:.1f}',
        'Import cost (MB)': f'{payload["import_cost_mb"]:.1f}',
        'External binary': 'yes' if payload['external_tool'] else 'no',
        'Status': 'OK',
      }, MEMORY_TABLE_COLUMNS)
    except Exception as e:  # noqa: BLE001 -- one method must not end the run
      detail = f'{type(e).__name__}: {e}'
      print(f'  !! {label} FAILED: {detail}')
      append_row(rows, output_path,
                 failed_row(label, 'FAILED', detail, MEMORY_TABLE_COLUMNS),
                 MEMORY_TABLE_COLUMNS)
      continue

  results = pd.DataFrame(rows, columns=MEMORY_TABLE_COLUMNS)

  print(f'\n\n{"=" * 80}')
  print(f'  {benchmark.upper()} MEMORY')
  print(f'{"=" * 80}\n')
  print(results.to_string(index=False))

  results.to_csv(output_path, sep='\t', index=False)
  print(f'\nSaved to {output_path}')

  failures = [r['Method'] for r in rows if r.get('Status', 'OK') != 'OK']
  if failures:
    print(f'!! {len(failures)} method(s) did not report memory: {", ".join(failures)}')
  else:
    print(f'All {len(rows)} methods reported memory successfully.')
  print('\nReported memory = peak RSS attributable to the method, EXCLUDING the shared '
        'Python interpreter baseline.\nRaw components are in the table so the reduction '
        'can be checked.')

  return results


def main():
  parser = argparse.ArgumentParser(description='PEPMatch Benchmarking Framework')
  parser.add_argument(
    '-b', '--benchmark',
    choices=['mhc_ligands', 'milk', 'coronavirus', 'neoepitopes', 'cosmic_indel', 'cedar_indel'],
    required=True,
  )
  parser.add_argument(
    '-m', '--memory', action='store_true', default=False,
    help='Also measure memory during the timed run. Re-runs every method, so it '
         'roughly doubles wall time; prefer --memory-only as a separate job.',
  )
  parser.add_argument(
    '--memory-only', action='store_true', default=False,
    help='Produce ONLY the memory table (no timing, no recall).',
  )
  parser.add_argument('-t', '--text_shifting', action='store_true', default=False)
  parser.add_argument(
    '-p', '--threads', type=int, default=1,
    help='Threads every tool is pinned to, for a fair timing comparison (default 1). '
         'On the cluster, set this to match --cpus-per-task.',
  )
  # Internal: measure one method's memory in a fresh process. Selected by INDEX, not
  # name, because a method can be registered twice (BLAST short vs default).
  parser.add_argument('--mem-index', type=int, default=None, help=argparse.SUPPRESS)
  args = parser.parse_args()

  if args.mem_index is not None:
    run_single_method_memory(args.benchmark, args.mem_index, args.threads)
    return

  if args.memory_only:
    run_memory_benchmark(args.benchmark, args.text_shifting, args.threads)
    return

  run_benchmark(args.benchmark, args.memory, args.text_shifting, args.threads)


if __name__ == '__main__':
  main()
