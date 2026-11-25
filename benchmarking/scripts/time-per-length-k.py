"""Where should per-length k stop? One k per batch against per-length k with a ceiling.

Variants are builds of PEPMatch, each an unpacked wheel in its own folder (--variant
NAME=FOLDER, repeated; the first is the reference the others are compared with):
  base    one k for the whole batch, taken from its shortest query (master + PRs #39-#42)
  capN    every query at its own k, k = min(N, length // (edits + 1)); cap16 is uncapped
A config runs in a fresh process with that folder first on PYTHONPATH, so all variants use
the same Python and the same polars and differ only in the pepmatch package. A worker
refuses to run if `pepmatch` was imported from anywhere else.

Query sets:
  the four benchmark sets on their own proteomes (mhc-ligands, neoepitopes, milk-peptides on
  human; coronavirus on betacoronaviruses)
  synthetic sets cut from the human proteome (fixed seed; half carry 1-2 substitutions):
    mixed-N    lengths uniform over 8-25, for N in --sizes (100 queries to 100,000: a small
               batch pays per-group overhead, a large one shows throughput)
  and at --fixed-size queries each:
    short-N    lengths uniform over 5-12 (the floor of the k rule; some lengths fall below
               the recall guarantee at 3 mismatches, in every variant alike)
    outlier-N  15-mers with 0.1% 8-mers (the short query that drags a one-k batch down)
    long-N     lengths uniform over 30-50 (where a ceiling on k matters most)
    only9-N, only15-N  single-length controls
Modes: exact, 1-3 mismatches, 1-2 indels, best match (best match only up to --best-max queries).

Per config and variant, two processes:
  cold  empty index folder, one search: total time, time inside the index builder, which
        indexes were built and their size on disk. This is what a first use costs.
  warm  the same folder: one untimed search, then timed ones: at least --repeats, and more
        (up to --max-repeats) until --min-timed-s seconds have been timed, so a 10 ms search
        is not judged on three samples. The median is reported, with the fastest, the
        slowest and the number of runs; plus peak memory and a checksum of the result rows.
The variants rotate which goes first. A process that runs past --cap-min is stopped and
reported as capped. Rows are compared with the first variant as a sorted multiset of (query,
matched sequence, protein, start, edits); best match is compared on each query's edit count,
since a tie may name a different protein.

Output: per-length-k-results.jsonl (everything), per-length-k-long.tsv (one row per config
and variant) and per-length-k-summary.tsv (one row per config, variants side by side).

  time-per-length-k.py --workdir W --outdir R --variant base=B --variant cap6=C ... [--smoke]
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
    while len(times) < args.repeats or (sum(times) < args.min_timed_s and len(times) < args.max_repeats):
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
      warm_s=statistics.median(times), warm_min_s=min(times), warm_max_s=max(times), warm_runs=len(times), warm_all=times,
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


def synthetic_sets(proteins, sizes, fixed, workdir):
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
  n = fixed
  write(f'short-{n}', [layout.randint(5, 12) for _ in range(n)])
  few = max(1, n // 1000)
  outlier = [15] * (n - few) + [8] * few
  layout.shuffle(outlier)
  write(f'outlier-{n}', outlier)
  write(f'long-{n}', [layout.randint(30, 50) for _ in range(n)])
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
    '--min-timed-s', str(args.min_timed_s), '--max-repeats', str(args.max_repeats),
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


def write_tables(results, names, outdir):
  reference = names[0]
  long_header = ['set', 'queries', 'mode', 'threads', 'variant', 'warm_s', 'warm_min_s', 'warm_max_s', 'warm_runs', 'speedup_vs_' + reference,
                 'cold_s', 'build_s', 'ks', 'indexes', 'index_gb', 'peak_rss_mb', 'hits', 'identical', 'note']
  with open(outdir / 'per-length-k-long.tsv', 'w', newline='') as f:
    out = csv.writer(f, delimiter='\t')
    out.writerow(long_header)
    for r in results:
      ref = r['variants'][reference]
      for name in names:
        v = r['variants'][name]
        speedup = ref['warm_s'] / v['warm_s'] if 'warm_s' in ref and v.get('warm_s') else None
        same = ('yes' if v['rows_sha256'] == ref['rows_sha256'] else 'NO') if 'rows_sha256' in v and 'rows_sha256' in ref else ''
        out.writerow([
          r['set'], r['queries'], r['mode'], r['threads'], name, fmt(v.get('warm_s')), fmt(v.get('warm_min_s')),
          fmt(v.get('warm_max_s')), v.get('warm_runs', ''), fmt(speedup, 2),
          fmt(v.get('cold_s'), 2), fmt(v.get('build_s'), 2), ','.join(map(str, v.get('ks', []))), len(v.get('ks', [])),
          fmt(v['index_bytes'] / 1e9, 2) if 'index_bytes' in v else '', fmt(v.get('peak_rss_mb'), 0),
          v.get('hits', ''), same, 'capped' if v.get('capped') else ('error' if v.get('error') else ''),
        ])
  with open(outdir / 'per-length-k-summary.tsv', 'w', newline='') as f:
    out = csv.writer(f, delimiter='\t')
    out.writerow(['set', 'queries', 'mode', 'threads']
                 + [f'{n}_warm_s' for n in names] + [f'{n}_cold_s' for n in names]
                 + [f'{n}_indexes' for n in names] + ['identical'])
    for r in results:
      vs = [r['variants'][n] for n in names]
      done = all('rows_sha256' in v for v in vs)
      differ = [n for n, v in zip(names, vs) if done and v['rows_sha256'] != vs[0]['rows_sha256']]
      unfinished = [n for n, v in zip(names, vs) if 'rows_sha256' not in v]
      out.writerow([r['set'], r['queries'], r['mode'], r['threads']]
                   + [fmt(v.get('warm_s')) for v in vs] + [fmt(v.get('cold_s'), 2) for v in vs]
                   + [len(v.get('ks', [])) for v in vs]
                   + ['yes' if done and not differ else 'NO ' + ','.join(differ + unfinished)])


def driver(args):
  args.workdir, args.outdir = Path(args.workdir).resolve(), Path(args.outdir).resolve()
  args.workdir.mkdir(parents=True, exist_ok=True)
  args.outdir.mkdir(parents=True, exist_ok=True)
  sites = {}
  for spec in args.variant:
    name, _, folder = spec.partition('=')
    sites[name] = Path(folder).resolve()
  names = list(sites)
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
    sizes, args.fixed_size, args.repeats, args.min_timed_s = [200], 200, 1, 0
    query_sets = {}
  else:
    proteins = read_proteins(proteomes['human'])
    query_sets = {name: (BENCH / query, proteomes[proteome]) for name, (query, proteome) in PAPER_SETS.items()}
  for name, path in synthetic_sets(proteins, sizes, args.fixed_size, args.workdir).items():
    query_sets[name] = (path, proteomes['human'])
  if args.sets:
    wanted = args.sets.split(',')
    query_sets = {name: v for name, v in query_sets.items() if any(name == w or name.startswith(w + '-') for w in wanted)}
  modes = args.modes.split(',') if args.modes else list(MODES)
  skip = {tuple(item.split(':')) for item in args.skip.split(',') if item}

  configs = []
  for name, (query, proteome) in query_sets.items():
    n = count_queries(query)
    for mode in modes:
      if mode == 'best' and n > args.best_max:
        continue
      if (name.split('-')[0], mode) in skip:
        continue
      for threads in threads_list:
        configs.append({'set': name, 'queries': n, 'mode': mode, 'threads': threads, 'query': query, 'proteome': proteome})

  info = [
    f'host        : {platform.node()}',
    f'python      : {sys.version.split()[0]}',
    f'cpus        : {os.cpu_count()}',
    f'threads     : {args.threads}',
    f'repeats     : at least {args.repeats}, until {args.min_timed_s} s timed, at most {args.max_repeats}; cap per process: {args.cap_min} min',
    f'human       : {len(proteins)} proteins, {sum(map(len, proteins))} residues ({proteomes["human"]})',
    f'variants    : {" ".join(names)}  (reference: {names[0]})',
    f'labels      : {args.labels}',
    f'configs     : {len(configs)}',
    f'started     : {time.strftime("%Y-%m-%dT%H:%M:%S%z")}',
  ]
  (args.outdir / 'run-info.txt').write_text('\n'.join(info) + '\n')
  print('\n'.join(info), flush=True)

  results = []
  shas = {name: set() for name in names}
  with open(args.outdir / 'per-length-k-results.jsonl', 'w') as log:
    for number, config in enumerate(configs):
      record = {k: config[k] for k in ('set', 'queries', 'mode', 'threads')}
      record['variants'] = {}
      shift = number % len(names)
      for name in names[shift:] + names[:shift]:
        index_dir = args.workdir / 'idx' / name
        if index_dir.exists():
          shutil.rmtree(index_dir)
        index_dir.mkdir(parents=True)
        result = run_worker(args, sites[name], 'cold', config, index_dir, config['threads'])
        if 'cold_s' in result:
          warm = run_worker(args, sites[name], 'warm', config, index_dir, config['threads'])
          shas[name].add(result.pop('so_sha256'))
          warm.pop('so_sha256', None)
          result.update(warm)
        record['variants'][name] = result
        shutil.rmtree(index_dir)
      results.append(record)
      log.write(json.dumps(record) + '\n')
      log.flush()
      write_tables(results, names, args.outdir)
      ref = record['variants'][names[0]]
      line = '  '.join(f"{n} {fmt(record['variants'][n].get('warm_s')) or 'FAIL'}" for n in names)
      same = all(record['variants'][n].get('rows_sha256') == ref.get('rows_sha256') for n in names) and 'rows_sha256' in ref
      print(f"[{number + 1}/{len(configs)}] {config['set']} {config['mode']} t={config['threads']}: {line}  "
            f"{'identical' if same else 'ROWS DIFFER OR FAILED'}", flush=True)

  bad = [r for r in results
         if not all(r['variants'][n].get('rows_sha256') == r['variants'][names[0]].get('rows_sha256')
                    and 'rows_sha256' in r['variants'][n] for n in names)]
  with open(args.outdir / 'run-info.txt', 'a') as f:
    f.write(f'finished    : {time.strftime("%Y-%m-%dT%H:%M:%S%z")}\n')
    for name in names:
      f.write(f'engine {name:5}: {sorted(shas[name])}\n')
    f.write(f'identical   : {len(results) - len(bad)} of {len(results)} configs (all variants)\n')
    f.write(f'not identical or capped/failed: {len(bad)}\n')
  print(f'\nidentical in {len(results) - len(bad)} of {len(results)} configs; not identical or failed: {len(bad)}')


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
  parser.add_argument('--variant', action='append', default=[], help='NAME=FOLDER; the first is the reference')
  parser.add_argument('--labels', default='', help='free text for run-info, e.g. the commits')
  parser.add_argument('--proteome', default=str(BENCH / 'proteomes' / 'human.fasta'))
  parser.add_argument('--threads', default=str(os.cpu_count()))
  parser.add_argument('--sizes', default='100,1000,10000,100000', help='query counts for the mixed sets')
  parser.add_argument('--fixed-size', type=int, default=10000, help='query count for the other synthetic sets')
  parser.add_argument('--min-timed-s', type=float, default=1.0)
  parser.add_argument('--max-repeats', type=int, default=15)
  parser.add_argument('--sets', default='', help='only these sets, by name or prefix (mixed, outlier, coronavirus, ...)')
  parser.add_argument('--modes', default='')
  parser.add_argument('--skip', default='short:mm3', help='set-prefix:mode pairs to leave out. short:mm3 is '
                      'left out by default: 5-12-mers at 3 mismatches match a large share of the proteome '
                      'and ran out of memory in every build (2026-11-16)')
  parser.add_argument('--repeats', type=int, default=3)
  parser.add_argument('--cap-min', type=float, default=30)
  parser.add_argument('--best-max', type=int, default=10000)
  parser.add_argument('--smoke', action='store_true')
  args = parser.parse_args()
  if args.worker:
    worker(args)
  else:
    for required in ('workdir', 'outdir'):
      if not getattr(args, required):
        parser.error(f'--{required} is required')
    if len(args.variant) < 2 or any('=' not in v for v in args.variant):
      parser.error('give at least two --variant NAME=FOLDER')
    driver(args)


if __name__ == '__main__':
  main()
