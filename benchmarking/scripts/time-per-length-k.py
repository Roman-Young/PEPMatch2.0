"""Per-length k against one k per batch: speed, first-use cost, memory and identical hits.

Two builds of PEPMatch are compared on one node, one config at a time:
  base  one k for the whole batch, taken from its shortest query (master + PRs #39-#42)
  head  every query at its own k, one index per k (branch test/per-length-k)

Each build is an unpacked wheel in its own folder (--base-site / --head-site); a config runs
in a fresh process with that folder first on PYTHONPATH, so both builds use the same Python
and the same polars and differ only in the pepmatch package. A worker refuses to run if
`pepmatch` was imported from anywhere else.

Query sets:
  the four benchmark sets on their own proteomes (mhc-ligands, neoepitopes, milk-peptides on
  human; coronavirus on betacoronaviruses)
  synthetic sets cut from the human proteome (fixed seed; half carry 1-2 substitutions):
    mixed-N    lengths uniform over 8-25, for N in --sizes
    outlier-N  15-mers with 0.1% 8-mers (the short query that drags a one-k batch down)
    only9-N, only15-N  single-length controls, where the two builds should tie
Modes: exact, 1-3 mismatches, 1-2 indels, best match (best match only up to --best-max queries).

Per config and build, two processes:
  cold  empty index folder, one search: total time, time inside the index builder, which
        indexes were built and their size on disk. This is what a first use costs.
  warm  the same folder: one untimed search, then --repeats timed ones (median reported),
        peak memory, and a checksum of the result rows.
The builds alternate which goes first. A process that runs past --cap-min is stopped and
reported as capped. Rows are compared between builds as a sorted multiset of (query, matched
sequence, protein, start, edits); best match is compared on each query's edit count, since a
tie may name a different protein.

  time-per-length-k.py --workdir W --outdir R --base-site B --head-site H [--smoke]
"""
import argparse
import contextlib
import csv
import hashlib
import io
import json
import os
import platform
import random
import re
import resource
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BENCH = REPO / 'benchmarking'
MODES = {
  'exact': {'max_mismatches': 0},
  'mm1': {'max_mismatches': 1},
  'mm2': {'max_mismatches': 2},
  'mm3': {'max_mismatches': 3},
  'indel1': {'max_indels': 1},
  'indel2': {'max_indels': 2},
  'best': {'best_match': True},
}
PAPER_SETS = {
  'mhc-ligands': ('queries/mhc-ligands-test.fasta', 'human'),
  'neoepitopes': ('queries/neoepitopes-test.fasta', 'human'),
  'milk-peptides': ('queries/milk-peptides-test.fasta', 'human'),
  'coronavirus': ('queries/coronavirus-test.fasta', 'betacoronaviruses'),
}
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
SEED = 20261116


# ---------------------------------------------------------------- worker (one build, one config)

def worker(args):
  site = Path(args.site).resolve()
  import pepmatch
  from pepmatch import Matcher, _rs
  import pepmatch.preprocessor as preprocessor
  if site not in Path(pepmatch.__file__).resolve().parents:
    sys.exit(f'FATAL: pepmatch imported from {pepmatch.__file__}, not from {site}')

  build = {'seconds': 0.0, 'count': 0}
  real_preprocess = preprocessor.rs_preprocess

  def timed_preprocess(*a, **kw):
    start = time.perf_counter()
    result = real_preprocess(*a, **kw)
    build['seconds'] += time.perf_counter() - start
    build['count'] += 1
    return result

  preprocessor.rs_preprocess = timed_preprocess

  def search():
    matcher = Matcher(
      query=args.query, proteome_file=args.proteome, preprocessed_files_path=args.index_dir,
      output_format='dataframe', **MODES[args.mode],
    )
    with contextlib.redirect_stdout(io.StringIO()):
      start = time.perf_counter()
      frame = matcher.match()
      return time.perf_counter() - start, frame

  out = {'so_sha256': hashlib.sha256(Path(_rs.__file__).read_bytes()).hexdigest()}
  if args.phase == 'cold':
    seconds, _ = search()
    indexes = sorted(Path(args.index_dir).glob('*.pepidx'), key=lambda p: int(re.search(r'_(\d+)mers', p.name)[1]))
    out.update(
      cold_s=seconds, build_s=build['seconds'], builds=build['count'],
      ks=[int(re.search(r'_(\d+)mers', p.name)[1]) for p in indexes],
      index_bytes=sum(p.stat().st_size for p in indexes),
    )
  else:
    search()   # untimed: pages the indexes in
    times = []
    for _ in range(args.repeats):
      seconds, frame = search()
      times.append(seconds)
    if build['count']:
      sys.exit('FATAL: a warm run built an index; the cold run did not leave them all behind')
    edits = 'Indels' if 'Indels' in frame.columns else 'Mismatches'
    hits = frame.filter(frame['Matched Sequence'].is_not_null())
    columns = ['Query Sequence', edits] if args.mode == 'best' else \
              ['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start', edits]
    rows = frame.select(columns).unique() if args.mode == 'best' else frame.select(columns)
    out.update(
      warm_s=statistics.median(times), warm_min_s=min(times), warm_all=times,
      rows=frame.height, hits=hits.height, queries_hit=hits['Query Sequence'].n_unique(),
      rows_sha256=hashlib.sha256(rows.sort(columns).write_csv().encode()).hexdigest(),
      peak_rss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
    )
  print('RESULT ' + json.dumps(out))


# ---------------------------------------------------------------- driver

def read_proteins(fasta, limit=None):
  proteins, current = [], []
  with open(fasta) as f:
    for line in f:
      if line.startswith('>'):
        if current:
          proteins.append(''.join(current))
          if limit and len(proteins) >= limit:
            return proteins
        current = []
      else:
        current.append(line.strip().upper())
  if current:
    proteins.append(''.join(current))
  return proteins


def cut_peptides(proteins, lengths, rng):
  """One peptide per entry of `lengths`, cut from a random protein; half get 1-2 substitutions."""
  peptides = []
  for length in lengths:
    while True:
      protein = rng.choice(proteins)
      if len(protein) < length:
        continue
      start = rng.randrange(len(protein) - length + 1)
      peptide = protein[start:start + length]
      if set(peptide) <= set(AMINO_ACIDS):
        break
    if rng.random() < 0.5:
      residues = list(peptide)
      for position in rng.sample(range(length), rng.choice((1, 2))):
        residues[position] = rng.choice(AMINO_ACIDS.replace(residues[position], ''))
      peptide = ''.join(residues)
    peptides.append(peptide)
  return peptides


def synthetic_sets(proteins, sizes, workdir):
  """{name: query file}. Written as .txt, one peptide per line."""
  folder = workdir / 'queries'
  folder.mkdir(parents=True, exist_ok=True)
  sets = {}

  def write(name, lengths):
    rng = random.Random(f'{SEED}-{name}')
    path = folder / f'{name}.txt'
    path.write_text('\n'.join(cut_peptides(proteins, lengths, rng)) + '\n')
    sets[name] = path

  layout = random.Random(SEED)
  for n in sizes:
    write(f'mixed-{n}', [layout.randint(8, 25) for _ in range(n)])
  n = sizes[min(1, len(sizes) - 1)]
  short = max(1, n // 1000)
  outlier = [15] * (n - short) + [8] * short
  layout.shuffle(outlier)
  write(f'outlier-{n}', outlier)
  write(f'only9-{n}', [9] * n)
  write(f'only15-{n}', [15] * n)
  return sets


def count_queries(path):
  with open(path) as f:
    if str(path).endswith('.txt'):
      return sum(1 for line in f if line.strip())
    return sum(1 for line in f if line.startswith('>'))


def run_worker(args, site, phase, config, index_dir, threads):
  env = dict(os.environ, PYTHONPATH=str(site), RAYON_NUM_THREADS=str(threads), POLARS_MAX_THREADS=str(threads))
  command = [
    sys.executable, str(Path(__file__).resolve()), '--worker', '--site', str(site), '--phase', phase,
    '--query', str(config['query']), '--proteome', str(config['proteome']), '--mode', config['mode'],
    '--index-dir', str(index_dir), '--repeats', str(args.repeats),
  ]
  try:
    done = subprocess.run(command, env=env, capture_output=True, text=True, timeout=args.cap_min * 60, cwd=str(args.workdir))
  except subprocess.TimeoutExpired:
    return {'capped': True}
  for line in done.stdout.splitlines():
    if line.startswith('RESULT '):
      return json.loads(line[7:])
  return {'error': (done.stderr or done.stdout).strip().splitlines()[-1:] or ['no output']}


def fmt(value, digits=3):
  return '' if value is None else f'{value:.{digits}f}'


def write_summary(results, path):
  header = ['set', 'queries', 'mode', 'threads',
            'base_warm_s', 'head_warm_s', 'warm_speedup',
            'base_cold_s', 'head_cold_s', 'base_build_s', 'head_build_s',
            'base_ks', 'head_ks', 'base_index_gb', 'head_index_gb',
            'base_rss_mb', 'head_rss_mb', 'hits', 'identical', 'note']
  with open(path, 'w', newline='') as f:
    out = csv.writer(f, delimiter='\t')
    out.writerow(header)
    for r in results:
      base, head = r['base'], r['head']
      notes = [f'{v} {k}' for v in ('base', 'head') for k in ('capped', 'error') if r[v].get(k)]
      both = 'rows_sha256' in base and 'rows_sha256' in head
      speedup = base['warm_s'] / head['warm_s'] if both and head['warm_s'] else None
      out.writerow([
        r['set'], r['queries'], r['mode'], r['threads'],
        fmt(base.get('warm_s')), fmt(head.get('warm_s')), fmt(speedup, 2),
        fmt(base.get('cold_s'), 2), fmt(head.get('cold_s'), 2), fmt(base.get('build_s'), 2), fmt(head.get('build_s'), 2),
        ','.join(map(str, base.get('ks', []))), ','.join(map(str, head.get('ks', []))),
        fmt(base['index_bytes'] / 1e9, 2) if 'index_bytes' in base else '',
        fmt(head['index_bytes'] / 1e9, 2) if 'index_bytes' in head else '',
        fmt(base.get('peak_rss_mb'), 0), fmt(head.get('peak_rss_mb'), 0),
        head.get('hits', base.get('hits', '')),
        ('yes' if base['rows_sha256'] == head['rows_sha256'] else 'NO') if both else '',
        '; '.join(notes),
      ])


def driver(args):
  args.workdir, args.outdir = Path(args.workdir).resolve(), Path(args.outdir).resolve()
  args.workdir.mkdir(parents=True, exist_ok=True)
  args.outdir.mkdir(parents=True, exist_ok=True)
  sites = {'base': Path(args.base_site).resolve(), 'head': Path(args.head_site).resolve()}
  threads_list = [int(t) for t in args.threads.split(',')]
  sizes = [int(s) for s in args.sizes.split(',')]

  proteomes = {'human': Path(args.proteome).resolve(), 'betacoronaviruses': BENCH / 'proteomes' / 'betacoronaviruses.fasta'}
  if args.smoke:
    # A small proteome and tiny sets: every code path in minutes, numbers meaningless.
    proteins = read_proteins(proteomes['human'], limit=3000)
    small = args.workdir / 'human-smoke.fasta'
    small.write_text(''.join(f'>sp|S{i:05d}|S{i:05d}_HUMAN Smoke protein OS=Homo sapiens OX=9606 GN=G{i} PE=1 SV=1\n{p}\n'
                             for i, p in enumerate(proteins)))
    proteomes['human'] = small
    sizes, args.repeats = [200], 1
    query_sets = {}
  else:
    proteins = read_proteins(proteomes['human'])
    query_sets = {name: (BENCH / query, proteomes[proteome]) for name, (query, proteome) in PAPER_SETS.items()}
  for name, path in synthetic_sets(proteins, sizes, args.workdir).items():
    query_sets[name] = (path, proteomes['human'])
  if args.sets:
    wanted = args.sets.split(',')
    query_sets = {name: v for name, v in query_sets.items() if any(name == w or name.startswith(w + '-') for w in wanted)}
  modes = args.modes.split(',') if args.modes else list(MODES)

  configs = []
  for name, (query, proteome) in query_sets.items():
    n = count_queries(query)
    for mode in modes:
      if mode == 'best' and n > args.best_max:
        continue
      for threads in threads_list:
        configs.append({'set': name, 'queries': n, 'mode': mode, 'threads': threads, 'query': query, 'proteome': proteome})

  info = [
    f'host        : {platform.node()}',
    f'python      : {sys.version.split()[0]}',
    f'cpus        : {os.cpu_count()}',
    f'threads     : {args.threads}',
    f'repeats     : {args.repeats}   cap per process: {args.cap_min} min',
    f'base        : {args.base_label}',
    f'head        : {args.head_label}',
    f'configs     : {len(configs)}',
    f'started     : {time.strftime("%Y-%m-%dT%H:%M:%S%z")}',
  ]
  (args.outdir / 'run-info.txt').write_text('\n'.join(info) + '\n')
  print('\n'.join(info), flush=True)

  results = []
  shas = {'base': set(), 'head': set()}
  with open(args.outdir / 'per-length-k-results.jsonl', 'w') as log:
    for number, config in enumerate(configs):
      record = {k: config[k] for k in ('set', 'queries', 'mode', 'threads')}
      order = ['base', 'head'] if number % 2 == 0 else ['head', 'base']
      for variant in order:
        index_dir = args.workdir / 'idx' / variant
        if index_dir.exists():
          shutil.rmtree(index_dir)
        index_dir.mkdir(parents=True)
        result = run_worker(args, sites[variant], 'cold', config, index_dir, config['threads'])
        if 'cold_s' in result:
          warm = run_worker(args, sites[variant], 'warm', config, index_dir, config['threads'])
          shas[variant].add(result.pop('so_sha256'))
          warm.pop('so_sha256', None)
          result.update(warm)
        record[variant] = result
        shutil.rmtree(index_dir)
      results.append(record)
      log.write(json.dumps(record) + '\n')
      log.flush()
      write_summary(results, args.outdir / 'per-length-k-summary.tsv')
      base, head = record['base'], record['head']
      same = ''
      if 'rows_sha256' in base and 'rows_sha256' in head:
        same = 'identical' if base['rows_sha256'] == head['rows_sha256'] else 'ROWS DIFFER'
      print(f"[{number + 1}/{len(configs)}] {config['set']} {config['mode']} t={config['threads']}: "
            f"base {fmt(base.get('warm_s')) or base}  head {fmt(head.get('warm_s')) or head}  {same}", flush=True)

  different = [r for r in results if 'rows_sha256' in r['base'] and 'rows_sha256' in r['head']
               and r['base']['rows_sha256'] != r['head']['rows_sha256']]
  unfinished = [r for r in results if 'rows_sha256' not in r['base'] or 'rows_sha256' not in r['head']]
  with open(args.outdir / 'run-info.txt', 'a') as f:
    f.write(f'finished    : {time.strftime("%Y-%m-%dT%H:%M:%S%z")}\n')
    f.write(f'base engine : {sorted(shas["base"])}\n')
    f.write(f'head engine : {sorted(shas["head"])}\n')
    f.write(f'compared    : {len(results) - len(unfinished)} of {len(results)} configs\n')
    f.write(f'rows differ : {len(different)}\n')
    f.write(f'capped/error: {len(unfinished)}\n')
  print(f'\ncompared {len(results) - len(unfinished)} of {len(results)}; rows differ in {len(different)}; capped or failed {len(unfinished)}')


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--worker', action='store_true')
  parser.add_argument('--site')
  parser.add_argument('--phase', choices=['cold', 'warm'])
  parser.add_argument('--query')
  parser.add_argument('--mode', choices=list(MODES))
  parser.add_argument('--index-dir')
  parser.add_argument('--workdir')
  parser.add_argument('--outdir')
  parser.add_argument('--base-site')
  parser.add_argument('--head-site')
  parser.add_argument('--base-label', default='')
  parser.add_argument('--head-label', default='')
  parser.add_argument('--proteome', default=str(BENCH / 'proteomes' / 'human.fasta'))
  parser.add_argument('--threads', default=str(os.cpu_count()))
  parser.add_argument('--sizes', default='1000,10000,100000')
  parser.add_argument('--sets', default='', help='only these sets, by name or prefix (mixed, outlier, coronavirus, ...)')
  parser.add_argument('--modes', default='')
  parser.add_argument('--repeats', type=int, default=3)
  parser.add_argument('--cap-min', type=float, default=30)
  parser.add_argument('--best-max', type=int, default=10000)
  parser.add_argument('--smoke', action='store_true')
  args = parser.parse_args()
  if args.worker:
    worker(args)
  else:
    for required in ('workdir', 'outdir', 'base_site', 'head_site'):
      if not getattr(args, required):
        parser.error(f'--{required.replace("_", "-")} is required')
    driver(args)


if __name__ == '__main__':
  main()
