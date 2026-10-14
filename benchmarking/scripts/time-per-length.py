#!/usr/bin/env python3
"""Per-query-length timing: how PEPMatch 2.0's search cost changes with query length.

WHY THIS EXISTS
  The 2-indel slowdown relative to 1 indel is explained in the text by an
  average-case seed-count argument: k = floor(L / (d + 1)), and every residue of seed
  lost makes a seed ~20x more common in the proteome. This script replaces the argument
  with a measurement. For each indel count d in {1, 2} and each length L in 8..25 it
  times PEPMatch on a query set where EVERY query has length L, so the seed length that
  set runs at is that length's own k. 2-indel 8-mers fall to k=2 (floor(8/3) = 2 < 3) and
  should be the expensive point.

TWO ARMS, SAME QUERIES
  per-length  k = the pigeonhole optimum for L (what a user searching only L-mers gets).
  fixed-k     k = 4 for 1 indel, 3 for 2 indels: the k the 9-25 aa scaling sets in
              Figure 4 actually ran at. Seed length is held constant, so length is the
              only variable. 2-indel 8-mers cannot run at k=3 (only two 3-mer tiles fit,
              and pigeonhole needs three), so that one point runs at k=2 in both arms.

TWO TIMERS, BOTH LOGGED
  engine   rs_indel_match() alone: the Rust search, nothing else. This is the figure.
  harness  Matcher(...).match() -> select -> to_pandas, the exact call benchmarking.py
           times for the PEPMatch "Searching (s)" column (pepmatch/benchmarker.py). It
           adds FASTA parsing, the metadata join and the per-hit indel annotation, which
           is single-threaded Python. Kept as the cross-check against Figure 4, never
           plotted as the algorithmic curve: at small N it is mostly overhead.

  Per-query cost is the MARGINAL cost (T_10k - T_1k) / 9,000, which cancels any fixed
  per-call cost. The 1k set is a byte-identical prefix of the 10k set (the generator
  nests subsets), so the difference is 9,000 more queries of the same population.
  Every (arm, d, L, N, timer) is measured REPEATS times; the summary takes the median.

  Index builds are timed once per k, before any search, and never charged to search.
  One untimed warm-up search per batch pages the index in, so the first timed repeat is
  not a cold-cache outlier.

RECALL
  Scored on every timed set through benchmarking.recall() (the harness's own function)
  against the generator's proteome-wide ground truth. Engine and harness run the same
  Rust call on the same queries, so engine hits are also checked to equal the harness's
  hit rows: a mismatch means the two timers did not measure the same work.

BRUTE FORCE CONTRAST
  methods/brute_force_naive.py (no seeding, no prefilter) on the first BF_N queries of
  each 1k set. Its Pool runs one query per worker, so each query is timed inside its
  worker: core-seconds per query is independent of how queries packed onto cores. The
  64-thread-equivalent per-query time is mean core-s / threads (the sweep is
  embarrassingly parallel). Never add a seed or prefilter step to make this faster.

PROVENANCE GUARDS (any failure exits non-zero before a single number is written)
  * PEPMATCH_SEED_STRATEGY=rarest is set, and the match.rs next to the IMPORTED pepmatch
    package routes indel seeding to rarest_pigeonhole_seeds.
  * The compiled engine (.so) is newer than match.rs -- Rust edits are invisible until
    `maturin develop --release`, the project's #1 footgun.
  * Optional --expect-so-sha256: the engine must be the exact binary the build job made.

Usage
  time-per-length.py --threads 64 --outdir results/per-length
  time-per-length.py --indels 2 --lengths 8,9 --sizes 100 --repeats 1 --bf-n 4   # smoke
"""
import argparse
import csv
import hashlib
import importlib.util
import os
import platform
import statistics
import subprocess
import sys
import time
from multiprocessing import Pool
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
REPO = BENCH.parent
LENGTHS = list(range(8, 26))
FIXED_K = {1: 4, 2: 3}

RUN_FIELDS = ['arm', 'indels', 'length', 'k', 'n', 'timer', 'repeat', 'seconds',
              'hit_rows', 'recall']
SUMMARY_FIELDS = [
  'arm', 'indels', 'length', 'k',
  'engine_marginal_ms_per_query', 'engine_ms_per_query_at_nmax',
  'harness_marginal_ms_per_query', 'harness_ms_per_query_at_nmax',
  'hit_rows_per_query', 'recall_min',
]
BF_FIELDS = ['indels', 'length', 'n', 'threads', 'proteome_load_s', 'wall_s',
             'core_s_mean', 'core_s_median', 'core_s_max', 'ms_per_query_at_threads',
             'recall']


def label_for(size):
  return f'{size // 1000}k' if size >= 1000 and size % 1000 == 0 else str(size)


def set_paths(prefix, d, L, size):
  stem = f'{prefix}-{d}indel-L{L:02d}-{label_for(size)}'
  return BENCH / 'queries' / f'{stem}.fasta', BENCH / 'expected' / f'{stem}-expected.tsv'


def own_k(L, d):
  # Mirrors Matcher.indel_search's partition: too short for 3-mer seeds -> its own k=2
  # table, otherwise the pigeonhole optimum.
  return 2 if L // (d + 1) < 3 else L // (d + 1)


def arm_k(arm, L, d):
  k = own_k(L, d)
  return k if arm == 'per-length' else min(FIXED_K[d], k)


def sha256(path, buf=1 << 20):
  h = hashlib.sha256()
  with open(path, 'rb') as f:
    while chunk := f.read(buf):
      h.update(chunk)
  return h.hexdigest()


def read_fasta(path):
  records, qid = [], None
  with open(path) as f:
    for line in f:
      line = line.strip()
      if line.startswith('>'):
        qid = line[1:].split()[0]
      elif line:
        records.append((qid, line))
  return records


def check_provenance(args):
  """Refuse to time an engine we cannot name. Returns facts for the run header."""
  import pepmatch
  from pepmatch import _rs
  pkg = Path(pepmatch.__file__).resolve().parent
  match_rs = pkg / 'rs-engine' / 'src' / 'match.rs'
  so = Path(_rs.__file__).resolve()
  problems = []

  strategy = os.environ.get('PEPMATCH_SEED_STRATEGY', '')
  src = match_rs.read_text() if match_rs.exists() else ''
  gated = 'Ok("rarest") => rarest_pigeonhole_seeds' in src
  unconditional = 'let seeds = rarest_pigeonhole_seeds(' in src
  if not (gated or unconditional):
    problems.append(f'{match_rs} does not route indel seeding to rarest_pigeonhole_seeds')
  if gated and strategy != 'rarest':
    problems.append('this engine needs PEPMATCH_SEED_STRATEGY=rarest (it is '
                    f'{strategy!r}); without it the full-tiling seeds are timed')
  if match_rs.exists() and so.stat().st_mtime < match_rs.stat().st_mtime:
    problems.append(f'{so.name} is older than match.rs: the engine is stale. Run '
                    f'`maturin develop --release` first.')
  so_sha = sha256(so)
  if args.expect_so_sha256 and so_sha != args.expect_so_sha256:
    problems.append(f'engine sha256 {so_sha[:12]} != the build job\'s '
                    f'{args.expect_so_sha256[:12]}: not the binary that was built')
  if problems:
    for p in problems:
      print(f'FATAL: {p}', file=sys.stderr)
    sys.exit(2)

  try:
    head = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=pkg,
                          capture_output=True, text=True, timeout=15).stdout.strip()
  except Exception:
    head = 'unknown'
  return {'pepmatch_pkg': str(pkg), 'git_head': head or 'unknown', 'engine_so': str(so),
          'engine_sha256': so_sha,
          'seed_strategy': 'rarest (gated)' if gated else 'rarest (unconditional)'}


def load_harness():
  sys.path.insert(0, str(BENCH))
  import benchmarking
  return benchmarking


def build_indexes(ks, proteome, workdir, out):
  """One index per k, timed once, before any search touches it."""
  from pepmatch.preprocessor import Preprocessor
  rows = []
  for k in sorted(ks):
    path = Path(workdir) / f'{Path(proteome).stem}_{k}mers.pepidx'
    t = time.perf_counter()
    Preprocessor(str(proteome), preprocessed_files_path=str(workdir)).preprocess(k=k)
    secs = time.perf_counter() - t
    rows.append({'k': k, 'seconds': f'{secs:.3f}', 'bytes': path.stat().st_size})
    print(f'  index k={k:<2}  {secs:7.2f} s  {path.stat().st_size / 1e6:8.1f} MB', flush=True)
  with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['k', 'seconds', 'bytes'], delimiter='\t')
    w.writeheader()
    w.writerows(rows)


def time_pepmatch(args, harness, proteome, workdir, run_writer, run_file):
  import pandas as pd
  from pepmatch import Matcher
  from pepmatch._rs import rs_indel_match

  for d in args.indels:
    for L in args.lengths:
      expected = {}
      for size in args.sizes:
        expected[size] = pd.read_csv(set_paths(args.prefix, d, L, size)[1], sep='\t')
      for arm in args.arms:
        k = arm_k(arm, L, d)
        pepidx = str(Path(workdir) / f'{Path(proteome).stem}_{k}mers.pepidx')
        warm = read_fasta(set_paths(args.prefix, d, L, args.sizes[0])[0])[:64]
        rs_indel_match(pepidx, warm, d)                     # untimed warm-up
        for size in args.sizes:
          fasta = set_paths(args.prefix, d, L, size)[0]
          queries = read_fasta(fasta)
          if len(queries) != size:
            sys.exit(f'FATAL: {fasta.name} has {len(queries)} queries, expected {size}')

          engine_hits = None
          for rep in range(1, args.repeats + 1):
            t = time.perf_counter()
            cols = rs_indel_match(pepidx, queries, d)
            secs = time.perf_counter() - t
            engine_hits = sum(m is not None for m in cols[2])
            run_writer.writerow(dict(arm=arm, indels=d, length=L, k=k, n=size,
                                     timer='engine', repeat=rep, seconds=f'{secs:.6f}',
                                     hit_rows=engine_hits, recall=''))
            run_file.flush()

          for rep in range(1, args.repeats + 1):
            t = time.perf_counter()
            df = Matcher(query=str(fasta), proteome_file=str(proteome), max_indels=d, k=k,
                         preprocessed_files_path=str(workdir), output_format='dataframe',
                         sequence_version=False).match()
            out = df.select(harness.RESULT_COLUMNS).to_pandas()
            secs = time.perf_counter() - t
            hit_rows = int(out['Matched Sequence'].notna().sum())
            rec = harness.recall(out, expected[size]) if rep == 1 else ''
            if rep == 1 and hit_rows != engine_hits:
              sys.exit(f'FATAL: d={d} L={L} {arm} N={size}: engine returned '
                       f'{engine_hits} hits, harness {hit_rows}. The timers measured '
                       f'different work.')
            if rep == 1 and rec != 100.0:
              sys.exit(f'FATAL: d={d} L={L} {arm} N={size}: recall {rec}, not 100.0')
            run_writer.writerow(dict(arm=arm, indels=d, length=L, k=k, n=size,
                                     timer='harness', repeat=rep, seconds=f'{secs:.6f}',
                                     hit_rows=hit_rows, recall=rec))
            run_file.flush()
          print(f'  d={d} L={L:2d} {arm:<10} k={k:<2} N={size:>6,}  '
                f'{engine_hits:>8,} hits  recall 100.0', flush=True)


def summarise(run_path, out_path, sizes):
  rows = list(csv.DictReader(open(run_path), delimiter='\t'))
  groups = {}
  for r in rows:
    key = (r['arm'], int(r['indels']), int(r['length']), int(r['k']))
    groups.setdefault(key, []).append(r)
  lo, hi = min(sizes), max(sizes)

  def med(rs, timer, n):
    vals = [float(r['seconds']) for r in rs if r['timer'] == timer and int(r['n']) == n]
    return statistics.median(vals) if vals else None

  out = []
  for (arm, d, L, k), rs in sorted(groups.items()):
    row = {'arm': arm, 'indels': d, 'length': L, 'k': k}
    for timer in ('engine', 'harness'):
      t_hi, t_lo = med(rs, timer, hi), med(rs, timer, lo)
      marginal = ((t_hi - t_lo) / (hi - lo) * 1e3
                  if t_hi is not None and t_lo is not None and hi > lo else None)
      row[f'{timer}_marginal_ms_per_query'] = '' if marginal is None else f'{marginal:.4f}'
      row[f'{timer}_ms_per_query_at_nmax'] = '' if t_hi is None else f'{t_hi / hi * 1e3:.4f}'
    hits = [int(r['hit_rows']) for r in rs if r['timer'] == 'engine' and int(r['n']) == hi]
    row['hit_rows_per_query'] = f'{hits[0] / hi:.3f}' if hits else ''
    recalls = [float(r['recall']) for r in rs if r['recall']]
    row['recall_min'] = f'{min(recalls):.1f}' if recalls else ''
    out.append(row)
  with open(out_path, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=SUMMARY_FIELDS, delimiter='\t')
    w.writeheader()
    w.writerows(out)
  return out


_BF = None


def _bf_timed(query):
  t = time.process_time()
  rows = _BF._search_one(query)
  return query, rows, time.process_time() - t


def time_brute_force(args, harness, proteome, out_path):
  """Naive brute force on the first bf_n queries of each smallest set, one node."""
  import pandas as pd
  global _BF
  spec = importlib.util.spec_from_file_location('bf_naive_timing',
                                                BENCH / 'methods' / 'brute_force_naive.py')
  _BF = importlib.util.module_from_spec(spec)
  sys.modules['bf_naive_timing'] = _BF       # Pool workers resolve the module by name
  spec.loader.exec_module(_BF)

  with open(out_path, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=BF_FIELDS, delimiter='\t')
    w.writeheader()
    for d in args.indels:
      tool = _BF.Benchmarker(benchmark='per-length', query='/dev/null',
                             proteome=str(proteome), lengths=[], max_mismatches=0,
                             method_parameters={'threads': args.threads}, indels=d)
      t = time.perf_counter()
      tool._load_proteome()                 # sets the module globals, incl. _N = d
      load_s = time.perf_counter() - t
      for L in args.lengths:
        fasta, exp_path = set_paths(args.prefix, d, L, min(args.sizes))
        queries = sorted({s for _, s in read_fasta(fasta)[:args.bf_n]})
        t = time.perf_counter()
        with Pool(min(args.threads, len(queries))) as pool:   # fork: shares the proteome
          results = pool.map(_bf_timed, queries, chunksize=1)
        wall = time.perf_counter() - t
        rows = [r for _, rs, _ in results for r in rs]
        core = [c for _, _, c in results]
        found = pd.DataFrame(rows, columns=harness.RESULT_COLUMNS)
        expected = pd.read_csv(exp_path, sep='\t')
        expected = expected[expected['Query Sequence'].isin(set(queries))]
        rec = harness.recall(found, expected)
        if rec != 100.0:
          sys.exit(f'FATAL: brute force d={d} L={L}: recall {rec}, not 100.0')
        mean = statistics.mean(core)
        w.writerow(dict(indels=d, length=L, n=len(queries), threads=args.threads,
                        proteome_load_s=f'{load_s:.2f}', wall_s=f'{wall:.2f}',
                        core_s_mean=f'{mean:.3f}',
                        core_s_median=f'{statistics.median(core):.3f}',
                        core_s_max=f'{max(core):.3f}',
                        ms_per_query_at_threads=f'{mean / args.threads * 1e3:.3f}',
                        recall=rec))
        f.flush()
        print(f'  BF d={d} L={L:2d}  {len(queries)} q  {mean:7.2f} core-s/q  '
              f'wall {wall:7.1f} s  recall 100.0', flush=True)


def main():
  p = argparse.ArgumentParser(description=__doc__,
                              formatter_class=argparse.RawDescriptionHelpFormatter)
  p.add_argument('--indels', default='1,2')
  p.add_argument('--lengths', default='8-25', help="e.g. '8-25' or '8,9,20'")
  p.add_argument('--sizes', default='1000,10000',
                 help='query counts per set; the marginal cost uses the smallest and largest')
  p.add_argument('--arms', default='per-length,fixed-k')
  p.add_argument('--repeats', type=int, default=3)
  p.add_argument('--threads', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK',
                                                                   os.cpu_count() or 1)))
  p.add_argument('--bf-n', type=int, default=64,
                 help='brute-force queries per length (0 skips brute force)')
  p.add_argument('--prefix', default='synth-len', help='dataset name prefix')
  p.add_argument('--proteome', default=None,
                 help='default: $PEPMATCH_BENCH_PROTEOME_DIR/human.fasta, else the repo copy')
  p.add_argument('--workdir', default='.', help='where the .pepidx files are built')
  p.add_argument('--outdir', default=str(BENCH / 'results' / 'per-length'))
  p.add_argument('--expect-so-sha256', default=None)
  args = p.parse_args()

  args.indels = [int(x) for x in args.indels.split(',')]
  if '-' in args.lengths:
    a, b = args.lengths.split('-')
    args.lengths = list(range(int(a), int(b) + 1))
  else:
    args.lengths = [int(x) for x in args.lengths.split(',')]
  args.sizes = sorted(int(x) for x in args.sizes.split(','))
  args.arms = [a.strip() for a in args.arms.split(',')]
  if not set(args.arms) <= {'per-length', 'fixed-k'}:
    sys.exit(f'unknown arm in {args.arms}')

  # Rayon reads this at first use, so it must be set before the engine runs anything.
  os.environ['RAYON_NUM_THREADS'] = str(args.threads)
  os.environ['OMP_NUM_THREADS'] = str(args.threads)

  facts = check_provenance(args)
  harness = load_harness()
  proteome_dir = os.environ.get('PEPMATCH_BENCH_PROTEOME_DIR')
  proteome = Path(args.proteome or (Path(proteome_dir) / 'human.fasta' if proteome_dir
                                    else BENCH / 'proteomes' / 'human.fasta'))
  outdir = Path(args.outdir)
  outdir.mkdir(parents=True, exist_ok=True)

  missing = [str(path) for d in args.indels for L in args.lengths for s in args.sizes
             for path in set_paths(args.prefix, d, L, s) if not path.exists()]
  if missing:
    print('FATAL: query sets not generated yet:', *missing[:10], sep='\n  ', file=sys.stderr)
    sys.exit(1)

  header = {**facts, 'host': platform.node(), 'threads': args.threads,
            'proteome': str(proteome), 'proteome_sha256': sha256(proteome),
            'indels': args.indels, 'lengths': args.lengths, 'sizes': args.sizes,
            'arms': args.arms, 'repeats': args.repeats, 'bf_n': args.bf_n,
            'started': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
  with open(outdir / 'run-info.txt', 'w') as f:
    for key, val in header.items():
      f.write(f'{key}\t{val}\n')
      print(f'  {key:<16} {val}')

  ks = {arm_k(arm, L, d) for arm in args.arms for d in args.indels for L in args.lengths}
  print('\n=== index builds (timed once, excluded from search) ===', flush=True)
  build_indexes(ks, proteome, args.workdir, outdir / 'index-build.tsv')

  print('\n=== PEPMatch 2.0 ===', flush=True)
  run_path = outdir / 'per-length-runs.tsv'
  with open(run_path, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=RUN_FIELDS, delimiter='\t')
    w.writeheader()
    time_pepmatch(args, harness, proteome, args.workdir, w, f)
  summary = summarise(run_path, outdir / 'per-length-summary.tsv', args.sizes)

  if args.bf_n > 0:
    print('\n=== naive brute force (contrast) ===', flush=True)
    time_brute_force(args, harness, proteome, outdir / 'per-length-bf.tsv')

  print('\n=== SUMMARY: engine ms/query (marginal), per arm ===')
  print(f'  {"arm":<10} {"d":>1} {"L":>3} {"k":>3} {"engine ms/q":>12} '
        f'{"harness ms/q":>13} {"hits/q":>8}')
  for r in summary:
    print(f'  {r["arm"]:<10} {r["indels"]:>1} {r["length"]:>3} {r["k"]:>3} '
          f'{r["engine_marginal_ms_per_query"]:>12} {r["harness_marginal_ms_per_query"]:>13} '
          f'{r["hit_rows_per_query"]:>8}')
  print(f'\nwrote {outdir}/  (run-info, index-build, per-length-runs, per-length-summary'
        f'{", per-length-bf" if args.bf_n > 0 else ""})')


if __name__ == '__main__':
  main()
