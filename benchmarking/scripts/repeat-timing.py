#!/usr/bin/env python3
"""Time every tool the same way, REPEATS times, and report the median search time.

One repeat of one tool is exactly one benchmarking.py run: build the tool's index or
database (timed, but not part of the search time), then search once. The aligner wrappers
delete their database at the end of search(), so every repeat has to rebuild it anyway;
PEPMatch follows the same build-then-search protocol, so all tools are measured the same
way. Repeats run round-robin (repeat 1 of every tool, then repeat 2, ...), so any drift on
the node during the job is spread across the tools instead of landing on one of them.

Uses benchmarking.py's own load_method / with_threads / pin_threads / recall, so a repeat
is configured identically to the single-run tables. Writes, per dataset, into --outdir:

  <dataset>_benchmarking.tsv   benchmarking.py's table layout with the MEDIAN of every
                               timing column (the figure scripts read it unchanged), plus
                               search min/max and every run
  <dataset>_runs.tsv           one row per run: build time, search time, recall, status

  repeat-timing.py --datasets cosmic_indel,cedar_indel --threads 32 \
    --methods 'PEPMatch,blast (short),blast (default),diamond,mmseqs2' --outdir results/x
"""
import argparse
import statistics
import sys
import traceback
from pathlib import Path

import pandas as pd

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH))
import benchmarking as bm   # noqa: E402  the harness's own helpers, not copies


def fmt(x):
  return 'N/A' if x is None else f'{x:.3f}'


def median_or_none(values):
  values = [v for v in values if v is not None]
  return statistics.median(values) if values else None


def summarize(runs, repeats):
  """One benchmarking.py-style row per method from its repeats."""
  rows = []
  for label in dict.fromkeys(r['Label'] for r in runs):
    ok = [r for r in runs if r['Label'] == label and r['Status'] == 'OK']
    mine = [r for r in runs if r['Label'] == label]
    if not ok:
      rows.append({'Method': label, 'Status': mine[-1]['Status'], 'Runs OK': f'0/{len(mine)}'})
      continue
    build = median_or_none([r['Build (s)'] for r in ok])
    query = median_or_none([r['Query preprocessing (s)'] for r in ok])
    search = statistics.median(r['Search (s)'] for r in ok)
    recalls = sorted({r['Recall (%)'] for r in ok})
    rows.append({
      'Method': ok[0]['Method'],
      'Proteome Preprocessing (s)': fmt(build),
      'Query Preprocessing (s)': fmt(query),
      'Searching (s)': fmt(search),
      'Total (s)': fmt(sum(t for t in (build, query, search) if t is not None)),
      'Memory (MB)': 'N/A',
      # Every tool here is deterministic, so a recall that moves between runs means
      # something is wrong; say so in the table instead of quietly picking one.
      'Recall (%)': f'{recalls[0]:.1f}' if len(recalls) == 1 else f'{recalls[0]:.1f}-{recalls[-1]:.1f}',
      'Status': 'OK' if len(ok) == repeats and len(recalls) == 1 else
                ('RECALL VARIES' if len(recalls) > 1 else 'INCOMPLETE'),
      'Runs OK': f'{len(ok)}/{len(mine)}',
      'Search min (s)': fmt(min(r['Search (s)'] for r in ok)),
      'Search max (s)': fmt(max(r['Search (s)'] for r in ok)),
      'All searches (s)': ' '.join(fmt(r['Search (s)']) for r in ok),
    })
  return pd.DataFrame(rows, columns=bm.RESULT_TABLE_COLUMNS +
                      ['Runs OK', 'Search min (s)', 'Search max (s)', 'All searches (s)'])


def main():
  p = argparse.ArgumentParser()
  p.add_argument('--datasets', required=True)
  p.add_argument('--methods', required=True,
                 help='comma-separated labels or names from benchmarking_parameters.json')
  p.add_argument('--repeats', type=int, default=5)
  p.add_argument('--threads', type=int, required=True)
  p.add_argument('--outdir', required=True)
  args = p.parse_args()

  bm.pin_threads(args.threads)   # before PEPMatch's engine first runs: rayon reads it once
  config = bm.load_config()
  wanted = [w.strip().lower() for w in args.methods.split(',') if w.strip()]
  methods = [m for m in config['methods'] if m.get('label', m['name']).lower() in wanted]
  missing = set(wanted) - {m.get('label', m['name']).lower() for m in methods}
  if missing:
    raise SystemExit(f'FATAL: no method labelled {sorted(missing)} in benchmarking_parameters.json')

  outdir = Path(args.outdir)
  outdir.mkdir(parents=True, exist_ok=True)
  failures = 0
  for name in args.datasets.split(','):
    ds = config['datasets'][name]
    if not (BENCH / ds['query']).exists():
      print(f'SKIP {name}: {ds["query"]} not generated here', flush=True)
      continue
    expected = pd.read_csv(BENCH / ds['expected'], sep='\t')
    runs = []
    for rep in range(1, args.repeats + 1):
      for m in methods:
        label = m.get('label', m['name'])
        run = {'Dataset': name, 'Label': label, 'Method': label, 'Repeat': rep, 'Build (s)': None,
               'Query preprocessing (s)': None, 'Search (s)': None, 'Recall (%)': None,
               'Status': 'OK', 'Detail': ''}
        try:
          tool = bm.load_method(m['name'], name, ds, bm.with_threads(m['method_parameters'], args.threads))
          if tool is None:
            run['Status'] = 'SKIPPED'
          else:
            run['Method'] = str(tool)
            run['Build (s)'], _ = bm.time_optional_step(tool.preprocess_proteome)
            run['Query preprocessing (s)'], _ = bm.time_optional_step(tool.preprocess_query)
            run['Search (s)'], df = bm.time_required_step(tool.search)
            run['Recall (%)'] = round(bm.recall(df, expected), 1)
            if hasattr(tool, 'cleanup'):
              tool.cleanup()
        except Exception as e:  # noqa: BLE001 -- one tool failing must not end the job
          run['Status'], run['Detail'] = 'FAILED', f'{type(e).__name__}: {e}'
          traceback.print_exc()
          failures += 1
        runs.append(run)
        print(f'{name}  rep {rep}/{args.repeats}  {label:16s}  build {fmt(run["Build (s)"])}  '
              f'search {fmt(run["Search (s)"])}  recall {run["Recall (%)"]}  {run["Status"]}', flush=True)
        # Written after every run, so a wall-clock kill keeps everything finished so far.
        pd.DataFrame(runs).to_csv(outdir / f'{name}_runs.tsv', sep='\t', index=False)
    table = summarize(runs, args.repeats)
    table.to_csv(outdir / f'{name}_benchmarking.tsv', sep='\t', index=False)
    print(f'\n--- {name} (median of {args.repeats}) ---\n{table.to_string(index=False)}\n', flush=True)
  if failures:
    print(f'!! {failures} run(s) FAILED -- see the traceback(s) above and *_runs.tsv')
    sys.exit(1)


if __name__ == '__main__':
  main()
