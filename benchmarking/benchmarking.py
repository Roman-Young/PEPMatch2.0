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

# Memory is split the same way the timing table is split, because preprocessing is a
# one-time amortised cost while search is paid per query. 'Peak total' is measured
# end-to-end rather than derived, since the true peak can exceed max(phases).
MEMORY_TABLE_COLUMNS = [
  'Method', 'Preprocessing (MB)', 'Search (MB)', 'Peak total (MB)',
  'External binary', 'Status',
]

# Full per-phase components, written beside the table so the baseline subtraction can
# be audited rather than taken on trust.
MEMORY_DETAIL_COLUMNS = [
  'Method', 'Phase', 'Memory (MB)', 'Peak RSS self (MB)', 'Peak RSS child (MB)',
  'Python baseline (MB)', 'Import cost (MB)', 'External binary',
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


def results_dir():
  """Directory for generated result tables. Kept separate from the framework/source
  files so runs don't clutter the top level; created on demand."""
  d = Path(__file__).parent / 'results'
  d.mkdir(exist_ok=True)
  return d


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


def measure_memory(benchmark, method_index, threads, phase='all'):
  """Peak memory for ONE method in ONE phase, measured in a fresh subprocess.

  Takes the method's INDEX rather than its name: a method may be registered more than
  once with different parameters (BLAST runs in both blastp-short and default modes),
  and selecting by name would silently measure whichever entry came first.

  Raises RuntimeError if the measurement did not produce a result, so a failure lands
  in the table as FAILED instead of a bare 'N/A' that reads like "not applicable".
  """
  proc = subprocess.run(
    [sys.executable, __file__, '-b', benchmark,
     '--mem-index', str(method_index), '--mem-phase', phase,
     '--threads', str(threads)],
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


def run_single_method_memory(benchmark, method_index, threads, phase='all'):
  """Internal entry point (--mem-index): run one method in a clean process and print its
  memory as a single `MEM_JSON={...}` line.

  Phases (each in a separate process, since peak RSS cannot be reset within one):
    'preprocess' -- build the index/database only; skips cleanup so the artifact survives.
    'search'     -- search against the artifact the preprocess phase left behind.
    'all'        -- authoritative end-to-end peak (measured, not derived: the true peak can
                    exceed max(preprocess, search) if preprocess memory isn't freed first).
  Methods declaring `search_requires_preprocess` (brute force -- proteome lives in-process
  with no on-disk handoff) report 'inseparable' for the split phases.

  The reported figure is peak RSS ABOVE the Python interpreter floor, so in-process tools
  (PEPMatch, brute force -- peak includes interpreter + pandas + polars) stay comparable to
  the aligners (standalone C binaries, no Python). Without it they weren't: PEPMatch and
  brute force previously reported an identical number, both dominated by the shared floor.
    * external binary  -> max(child peak, self peak - baseline)   both terms floor-free
    * in-process       -> max(self peak, child peak) - baseline   a fork inherits the floor
  Every raw component is reported alongside so the reduction can be audited, not trusted.
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

    if getattr(tool, 'search_requires_preprocess', False) and phase != 'all':
      print('MEM_JSON=' + json.dumps({'inseparable': True, 'label': label}))
      return

    # Same step semantics as the timed run: preprocessing may be inapplicable, but a
    # crash in search is a real failure and must propagate.
    if phase in ('all', 'preprocess'):
      for step in (tool.preprocess_proteome, tool.preprocess_query):
        try:
          step()
        except (StepNotApplicable, TypeError):
          pass

    if phase in ('all', 'search'):
      tool.search()
      # Cleanup only ever runs in a phase that searched, so the 'preprocess' phase
      # leaves its index/database in place for the 'search' phase to consume.
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
      'phase': phase,
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
  output_path = str(results_dir() / f'{benchmark}_benchmarking.tsv')
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


def phase_cell(phase_results, phase):
  """Render one phase's figure for the memory table, distinguishing 'this method cannot
  be split' from 'this phase produced no number'."""
  p = phase_results.get(phase, {})
  if p.get('inseparable'):
    return 'inseparable'
  return f'{p["memory_mb"]:.1f}' if 'memory_mb' in p else 'N/A'


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
  detail_rows = []
  output_path = str(results_dir() / f'{benchmark}_memory.tsv')
  detail_path = str(results_dir() / f'{benchmark}_memory_detail.tsv')
  Path(output_path).unlink(missing_ok=True)
  Path(detail_path).unlink(missing_ok=True)

  for method_index, method in methods:
    label = method.get('label', method['name'])
    print(f'\n{"=" * 60}')
    print(f'  {label}  (memory)')
    print(f'{"=" * 60}')
    sys.stdout.flush()

    try:
      # Order matters: 'preprocess' leaves its index/database behind for 'search' to
      # consume, then 'all' rebuilds from scratch for the authoritative end-to-end peak.
      phase_results = {}
      for phase in ('preprocess', 'search', 'all'):
        payload = measure_memory(benchmark, method_index, threads, phase)
        phase_results[phase] = payload

        if payload.get('skipped'):
          print(f'  {phase:<10} -> SKIPPED (not applicable to this dataset)')
          break
        if payload.get('inseparable'):
          print(f'  {phase:<10} -> inseparable (search needs preprocessing in-process)')
          continue

        print(f'  {phase:<10} -> {payload["memory_mb"]:.1f} MB '
              f'(self {payload["self_peak_mb"]:.1f} / child {payload["child_peak_mb"]:.1f} / '
              f'baseline {payload["baseline_mb"]:.1f})')
        append_row(detail_rows, detail_path, {
          'Method': payload['label'],
          'Phase': phase,
          'Memory (MB)': f'{payload["memory_mb"]:.1f}',
          'Peak RSS self (MB)': f'{payload["self_peak_mb"]:.1f}',
          'Peak RSS child (MB)': f'{payload["child_peak_mb"]:.1f}',
          'Python baseline (MB)': f'{payload["baseline_mb"]:.1f}',
          'Import cost (MB)': f'{payload["import_cost_mb"]:.1f}',
          'External binary': 'yes' if payload['external_tool'] else 'no',
        }, MEMORY_DETAIL_COLUMNS)
        sys.stdout.flush()

      if phase_results.get('preprocess', {}).get('skipped'):
        append_row(rows, output_path,
                   failed_row(label, 'SKIPPED', 'not applicable to this dataset',
                              MEMORY_TABLE_COLUMNS),
                   MEMORY_TABLE_COLUMNS)
        continue

      total = phase_results['all']
      append_row(rows, output_path, {
        'Method': total['label'],
        'Preprocessing (MB)': phase_cell(phase_results, 'preprocess'),
        'Search (MB)': phase_cell(phase_results, 'search'),
        'Peak total (MB)': f'{total["memory_mb"]:.1f}',
        'External binary': 'yes' if total['external_tool'] else 'no',
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
  pd.DataFrame(detail_rows, columns=MEMORY_DETAIL_COLUMNS).to_csv(
    detail_path, sep='\t', index=False
  )
  print(f'\nSaved to {output_path}')
  print(f'Per-phase components in {detail_path}')

  failures = [r['Method'] for r in rows if r.get('Status', 'OK') != 'OK']
  if failures:
    print(f'!! {len(failures)} method(s) did not report memory: {", ".join(failures)}')
  else:
    print(f'All {len(rows)} methods reported memory successfully.')
  print('\nReported memory = peak RSS attributable to the method, EXCLUDING the shared '
        'Python interpreter baseline.')
  print('Preprocessing and Search are measured in separate processes (peak RSS cannot '
        'be reset within one).')
  print('Peak total is measured end-to-end, not derived, because the true peak can '
        'exceed max(phases).')

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
  parser.add_argument('--mem-phase', choices=['all', 'preprocess', 'search'],
                      default='all', help=argparse.SUPPRESS)
  args = parser.parse_args()

  if args.mem_index is not None:
    run_single_method_memory(args.benchmark, args.mem_index, args.threads, args.mem_phase)
    return

  if args.memory_only:
    run_memory_benchmark(args.benchmark, args.text_shifting, args.threads)
    return

  run_benchmark(args.benchmark, args.memory, args.text_shifting, args.threads)


if __name__ == '__main__':
  main()
