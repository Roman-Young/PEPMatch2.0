"""Timing for the short-target search: Matcher(short_targets=True).

With short_targets on, after the normal (forward) search PEPMatch also reports proteome
entries shorter than the query that lie inside it, by searching those entries against a
small index of the queries. This script measures what that costs and checks that it finds
what it should, on one node, one config at a time.

Tiers (--tiers, default all):
  A  human proteome x the 4 benchmark query sets x 5 modes x flag off/on
  B  synthetic epitope DB (N entries of 8-25 aa cut from human.fasta) x Q queries of
     15-20 aa x 5 modes x off/on, for N in --sizes-b and Q in --q-sizes
  C  the same at N = --n-c (1M by default) for --q-c queries
  D  k rule: the shipped rule (forward's: one k from the shortest entry, indel short/rest
     split) vs per-length k, at N = --n-de, Q = --q-de. Per-length runs the shipped
     reversed pass once per entry length, so each pass picks that length's own k with
     the shipped code itself; its query index is rebuilt per length (milliseconds).
  E  minimum entry length (MIN_TARGET_LEN 2 / 5 / 8) on the tier-D DB and on human
  O  oracle: a query subsample against the smallest tier-B DB, compared with the brute
     force definition used by pepmatch/tests/test_short_targets_property.py

Each config runs in its own process (so peak RSS is that config's), with one untimed
warm-up and --repeats timed runs; the median is reported. Proteome index builds happen in
the warm-up, are timed separately and never count as search time. Stage times come from
wrapping the engine calls the matcher makes; forward and reversed calls are told apart by
the index they read (the reversed pass reads its temporary queries index).

Checks: flag-on forward rows == flag-off rows for the same config; every planted entry
whose edits fit the mode's budget is found; tier O has zero false and zero missed rows
for entries long enough for the seed guarantee (>= 2 * (edits + 1) residues).

  time-short-targets.py --workdir W --outdir R [--smoke] [--expect-so-sha256 SHA]
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
import resource
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
  'indel1': {'max_indels': 1},
  'indel2': {'max_indels': 2},
}
EDITS = {'exact': ('mm', 0), 'mm1': ('mm', 1), 'mm2': ('mm', 2), 'indel1': ('indel', 1), 'indel2': ('indel', 2)}
QUERY_SETS = ['mhc-ligands', 'neoepitopes', 'milk-peptides', 'coronavirus']
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
NOISE_LEN = 8   # reversed hits on entries shorter than this are reported as noise
FIELDS = ['tier', 'dataset', 'n_entries', 'n_queries', 'mode', 'flag', 'variant', 'status',
          'total_s', 'forward_s', 'reversed_s', 'pick_s', 'qindex_s', 'rev_search_s',
          'translate_merge_s', 'index_build_s', 'repeats_s', 'peak_mb', 'hit_rows',
          'forward_rows', 'added_rows', 'noise_rows', 'planted_in_budget', 'planted_found',
          'note']


def sha256(path, buf=1 << 20):
  h = hashlib.sha256()
  with open(path, 'rb') as f:
    while chunk := f.read(buf):
      h.update(chunk)
  return h.hexdigest()


def read_fasta(path):
  records, name, seq = [], None, []
  with open(path) as f:
    for line in f:
      line = line.strip()
      if line.startswith('>'):
        if name is not None:
          records.append((name, ''.join(seq)))
        name, seq = line[1:].split()[0], []
      elif line:
        seq.append(line)
  if name is not None:
    records.append((name, ''.join(seq)))
  return records


# ── data ─────────────────────────────────────────────────────────────────────────────

def make_epitope_db(proteome, n, seed, path):
  """N windows of human proteins, 8-15 aa for ~80% and 16-25 aa for the rest, in UniProt
  header form so the metadata join runs as it would on a real index. Returns the entries
  with where each came from, for planting."""
  rng = random.Random(seed)
  pool = [s for _, s in proteome if len(s) >= 40]
  entries = []
  with open(path, 'w') as f:
    for i in range(n):
      length = rng.randint(8, 15) if rng.random() < 0.8 else rng.randint(16, 25)
      prot = rng.choice(pool)
      start = rng.randrange(len(prot) - length + 1)
      name = f'EPI{i:07d}'
      f.write(f'>sp|{name}|{name} synthetic epitope OS=Homo sapiens OX=9606 GN={name} PE=1 SV=1\n'
              f'{prot[start:start + length]}\n')
      entries.append((name, prot, start, length))
  return entries


def _edit(rng, region, kind, n):
  """`region` with n interior edits: substitutions ('sub'), or insertions or deletions
  ('indel', one kind per region, since 2-indel search is homogeneous)."""
  if kind == 'sub':
    for _ in range(n):
      i = rng.randrange(len(region))
      region = region[:i] + rng.choice([a for a in AMINO_ACIDS if a != region[i]]) + region[i + 1:]
  elif rng.random() < 0.5:
    for _ in range(n):
      i = rng.randrange(1, len(region) - 1)
      region = region[:i] + region[i + 1:]
  else:
    for _ in range(n):
      i = rng.randrange(1, len(region))
      region = region[:i] + rng.choice(AMINO_ACIDS) + region[i:]
  return region


def make_queries(entries, proteome, q, seed, path):
  """Q queries of 15-20 aa: ~30% planted around a DB entry of at most 17 aa (exact, 1-2
  substitutions or 1-2 indels inside the entry's stretch, flanked by its own protein),
  the rest random human windows. The same queries serve every mode; the plant list says
  which modes must find each: those whose edits fit the budget and whose entry is still
  shorter than the query."""
  rng = random.Random(seed)
  pool = [s for _, s in proteome if len(s) >= 40]
  plantable = [e for e in entries if e[3] <= 17]
  kinds = [('exact', 0)] * 30 + [('sub', 1)] * 20 + [('sub', 2)] * 15 + [('indel', 1)] * 20 + [('indel', 2)] * 15
  plants = []
  with open(path, 'w') as f:
    for i in range(q):
      qid = f'Q{i:06d}'
      if plantable and rng.random() < 0.3:
        name, prot, start, length = rng.choice(plantable)
        kind, n = rng.choice(kinds)
        qlen = rng.randint(max(15, length + 1), 20)
        # Any window of qlen that holds the whole entry: lo in [start+length-qlen, start].
        lo = rng.randint(max(0, start + length - qlen), min(start, len(prot) - qlen))
        window, a = prot[lo:lo + qlen], start - lo
        region = window[a:a + length]
        if kind != 'exact':
          region = _edit(rng, region, kind, n)
        seq = window[:a] + region + window[a + length:]
        # Deletions can leave the query no longer than the entry; that pair is the forward
        # search's case, not containment, so the length pair is kept to judge each plant.
        plants.append((qid, name, kind, n, length, len(seq)))
      else:
        prot = rng.choice(pool)
        qlen = rng.randint(15, 20)
        s = rng.randrange(len(prot) - qlen + 1)
        seq = prot[s:s + qlen]
      f.write(f'>{qid}\n{seq}\n')
  with open(str(path) + '.plants.json', 'w') as f:
    json.dump(plants, f)


def in_budget(mode, kind, n):
  mtype, edits = EDITS[mode]
  if kind == 'exact':
    return True
  return (kind == 'sub' and mtype == 'mm' and n <= edits) or (kind == 'indel' and mtype == 'indel' and n <= edits)


# ── one config, in its own process ───────────────────────────────────────────────────

def run_child(cfg):
  import polars as pl
  import pepmatch.matcher as M
  from pepmatch import Matcher

  stage = {}

  def add(key, secs):
    stage[key] = stage.get(key, 0.0) + secs

  def wrap(name, classify):
    real = getattr(M, name)

    def timed(*args, **kwargs):
      t = time.perf_counter()
      try:
        return real(*args, **kwargs)
      finally:
        add(classify(args), time.perf_counter() - t)
    setattr(M, name, timed)
    return real

  def is_queries_index(path):
    return Path(str(path)).name.startswith('queries')
  wrap('rs_preprocess', lambda a: 'qindex' if is_queries_index(a[2]) else 'index_build')
  wrap('rs_match', lambda a: 'rev_search' if is_queries_index(a[0]) else 'forward_engine')
  wrap('rs_indel_match', lambda a: 'rev_search' if is_queries_index(a[0]) else 'forward_engine')
  real_short = wrap('rs_short_proteins', lambda a: 'pick')
  wrapped_short = M.rs_short_proteins

  counts = {}
  real_search, real_merge = Matcher._short_target_search, Matcher._merge_short_targets

  def per_length(self, pepidx_path):
    # The shipped pass, once per entry length: its own k rule then picks that length's k.
    longest = max(len(s) for _, s in self.query)
    _, seqs = real_short(pepidx_path, M.MIN_TARGET_LEN, longest)
    frames = []
    for length in sorted(set(map(len, seqs))):
      M.rs_short_proteins = lambda p, lo, hi, L=length: real_short(p, L, L + 1)
      try:
        df = real_search(self, pepidx_path)
      finally:
        M.rs_short_proteins = wrapped_short
      if df is not None:
        frames.append(df)
    return pl.concat(frames) if frames else None

  def timed_search(self, pepidx_path):
    t = time.perf_counter()
    try:
      return per_length(self, pepidx_path) if cfg['variant'] == 'per-length' else real_search(self, pepidx_path)
    finally:
      add('reversed_search_total', time.perf_counter() - t)

  def timed_merge(self, forward, short):
    t = time.perf_counter()
    counts['forward_rows'] = forward.filter(pl.col('Matched Sequence').is_not_null()).height
    try:
      return real_merge(self, forward, short)
    finally:
      add('merge', time.perf_counter() - t)
  Matcher._short_target_search, Matcher._merge_short_targets = timed_search, timed_merge

  if cfg.get('floor'):
    M.MIN_TARGET_LEN = cfg['floor']

  def once():
    stage.clear()
    counts.clear()
    t = time.perf_counter()
    with contextlib.redirect_stdout(io.StringIO()):
      df = Matcher(query=cfg['queries'], proteome_file=cfg['proteome'], preprocessed_files_path=cfg['index_dir'],
                   output_format='dataframe', short_targets=cfg['flag'] == 'on', **MODES[cfg['mode']]).match()
    return df, time.perf_counter() - t, dict(stage), dict(counts)

  df, _, warm, warm_counts = once()   # builds any missing proteome index; not timed as search
  totals, stages = [], []
  for _ in range(cfg['repeats']):
    _, secs, st, _ = once()
    totals.append(secs)
    stages.append(st)

  def med(key):
    return statistics.median(s.get(key, 0.0) for s in stages)
  total = statistics.median(totals)
  rev_total = med('reversed_search_total') + med('merge')
  hits = df.filter(pl.col('Matched Sequence').is_not_null())
  row = dict(status='ok', total_s=total, forward_s=total - rev_total, reversed_s=rev_total, pick_s=med('pick'),
             qindex_s=med('qindex'), rev_search_s=med('rev_search'),
             translate_merge_s=rev_total - med('pick') - med('qindex') - med('rev_search'),
             index_build_s=warm.get('index_build', 0.0), repeats_s=';'.join(f'{t:.4f}' for t in totals),
             peak_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, hit_rows=hits.height)
  if cfg['flag'] == 'on':
    forward = warm_counts.get('forward_rows', hits.height)
    added = hits.slice(forward)
    # A search that hit nothing leaves the column typed Null, which has no string methods.
    entry_len = pl.col('Matched Sequence').cast(pl.Utf8).str.len_chars()
    row.update(forward_rows=forward, added_rows=added.height, noise_rows=added.filter(entry_len < NOISE_LEN).height)
  if cfg.get('plants'):
    with open(cfg['plants']) as f:
      plants = [p for p in json.load(f) if in_budget(cfg['mode'], p[2], p[3]) and p[4] < p[5]]
    found = {(qid, pid.split('.')[0]) for qid, pid in hits.select('Query ID', 'Protein ID').rows()}
    row.update(planted_in_budget=len(plants), planted_found=sum((p[0], p[1]) in found for p in plants))
  print('RESULT ' + json.dumps(row))


# ── oracle subsample ─────────────────────────────────────────────────────────────────

def _oracle_one(args):
  sys.path.insert(0, str(REPO / 'pepmatch' / 'tests'))
  from test_short_targets_property import expected_rows
  query, entries, mode = args
  return expected_rows([query], entries, mode)


def run_oracle(queries_path, db_path, mode, n, workdir, threads):
  """Brute-force definition vs the engine, n queries against the whole DB."""
  from multiprocessing import Pool
  from pepmatch import Matcher
  queries = [s for _, s in read_fasta(queries_path)][:n]
  entries = {header.split('|')[1]: seq for header, seq in read_fasta(db_path)}
  with contextlib.redirect_stdout(io.StringIO()):
    df = Matcher(query=queries, proteome_file=str(db_path), preprocessed_files_path=str(workdir / 'index'),
                 output_format='dataframe', short_targets=True, **MODES[mode]).match()
  found = {(r['Query Sequence'], r['Protein ID'], r['Matched Sequence'], r['Index start'], r['Index end'],
            r['Query start'], r['Query end']) for r in df.iter_rows(named=True) if r['Matched Sequence'] is not None}
  with Pool(threads) as pool:
    parts = pool.map(_oracle_one, [(q, entries, mode) for q in queries])
  expected = {}
  for part in parts:
    for row, guaranteed in part.items():
      expected[row] = expected.get(row, False) or guaranteed
  missed = {row for row, guaranteed in expected.items() if guaranteed} - found
  return dict(mode=mode, queries=len(queries), engine_rows=len(found), oracle_rows=len(expected),
              false_rows=len(found - set(expected)), missed_rows=len(missed))


# ── driver ───────────────────────────────────────────────────────────────────────────

def provenance(args):
  import pepmatch
  import polars
  from pepmatch import _rs
  so = Path(_rs.__file__).resolve()
  so_sha = sha256(so)
  if args.expect_so_sha256 and so_sha != args.expect_so_sha256:
    sys.exit(f'FATAL: engine sha256 {so_sha[:12]} != the build job\'s {args.expect_so_sha256[:12]}')
  if not hasattr(_rs, 'rs_short_proteins'):
    sys.exit('FATAL: this engine has no rs_short_proteins: not the short-targets branch')
  head = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=REPO, capture_output=True, text=True).stdout.strip()
  cpu = ''
  with contextlib.suppress(Exception):
    cpu = next(line.split(':', 1)[1].strip() for line in open('/proc/cpuinfo') if line.startswith('model name'))
  return {'host': platform.node(), 'cpu': cpu, 'cpus_available': len(os.sched_getaffinity(0)),
          'python': platform.python_version(), 'polars': polars.__version__,
          'pepmatch': str(Path(pepmatch.__file__).parent), 'git_head': head, 'engine_sha256': so_sha,
          'started': time.strftime('%Y-%m-%dT%H:%M:%S%z')}


def main():
  ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
  ap.add_argument('--workdir', type=Path, required=True)
  ap.add_argument('--outdir', type=Path, required=True)
  ap.add_argument('--proteome', type=Path, default=BENCH / 'proteomes' / 'human.fasta')
  ap.add_argument('--queries-dir', type=Path, default=BENCH / 'queries')
  ap.add_argument('--tiers', default='A,B,D,E,C,O')
  ap.add_argument('--sizes-b', default='10000,100000')
  ap.add_argument('--q-sizes', default='100,1000,10000')
  ap.add_argument('--n-c', type=int, default=1_000_000)
  ap.add_argument('--q-c', default='1000')
  ap.add_argument('--modes-c', default=','.join(MODES))
  ap.add_argument('--n-de', type=int, default=100_000)
  ap.add_argument('--q-de', type=int, default=1000)
  ap.add_argument('--oracle-n', type=int, default=200)
  ap.add_argument('--repeats', type=int, default=3)
  ap.add_argument('--cap-min', type=float, default=45)
  ap.add_argument('--seed', type=int, default=20261106)
  ap.add_argument('--threads', type=int, default=len(os.sched_getaffinity(0)))
  ap.add_argument('--smoke', action='store_true', help='tiny sizes, 1 repeat, a ~3,500-protein human subset')
  ap.add_argument('--expect-so-sha256', default='')
  ap.add_argument('--child', help=argparse.SUPPRESS)
  args = ap.parse_args()
  if args.child:
    return run_child(json.loads(args.child))

  if args.smoke:
    args.sizes_b, args.q_sizes, args.n_c, args.q_c = '1000', '50,200', 5000, '100'
    args.n_de, args.q_de, args.oracle_n, args.repeats, args.cap_min = 2000, 100, 20, 1, 10
  args.workdir.mkdir(parents=True, exist_ok=True)
  args.outdir.mkdir(parents=True, exist_ok=True)
  info = provenance(args)
  index_dir = args.workdir / 'index'
  index_dir.mkdir(exist_ok=True)
  tiers = set(args.tiers.split(','))

  proteome = read_fasta(args.proteome)
  data = {'proteome': (str(args.proteome), sha256(args.proteome))}
  human_path = args.proteome
  query_set_paths = {qs: args.queries_dir / f'{qs}-test.fasta' for qs in QUERY_SETS}
  if args.smoke:   # small human index that keeps every short protein, so the reversed pass has work
    human_path = args.workdir / 'human-smoke.fasta'
    keep = [r for r in proteome if len(r[1]) < 26] + [r for r in proteome if len(r[1]) >= 26][:3000]
    with open(human_path, 'w') as f:
      for name, seq in keep:
        f.write(f'>{name} OS=Homo sapiens OX=9606 PE=1 SV=1\n{seq}\n')
    for qs, path in list(query_set_paths.items()):
      small = args.workdir / f'{qs}-smoke.fasta'
      with open(small, 'w') as f:
        for name, seq in read_fasta(path)[:50]:
          f.write(f'>{name}\n{seq}\n')
      query_set_paths[qs] = small

  dbs, qsets = {}, {}

  def db(n):
    if n not in dbs:
      path = args.workdir / f'epitopes-{n}.fasta'
      dbs[n] = (path, make_epitope_db(proteome, n, args.seed + n, path))
      data[path.name] = (str(path), sha256(path))
    return dbs[n]

  def queries(n, q):
    if (n, q) not in qsets:
      qpath = args.workdir / f'queries-{n}-{q}.fasta'
      make_queries(db(n)[1], proteome, q, args.seed + 7 * n + q, qpath)
      qsets[(n, q)] = qpath
      data[qpath.name] = (str(qpath), sha256(qpath))
    return qsets[(n, q)]

  configs = []

  def add(tier, dataset, n, q, mode, flag, proteome_path, query_path, plants=None, variant='shipped', floor=None):
    configs.append(dict(tier=tier, dataset=dataset, n_entries=n, n_queries=q, mode=mode, flag=flag, variant=variant,
                        floor=floor, proteome=str(proteome_path), queries=str(query_path), plants=plants,
                        index_dir=str(index_dir), repeats=args.repeats))

  if 'A' in tiers:
    for qs in QUERY_SETS:
      for mode in MODES:
        for flag in ('off', 'on'):
          add('A', f'human/{qs}', 'human', qs, mode, flag, human_path, query_set_paths[qs])
  for tier, sizes, qsizes, modes in (('B', args.sizes_b, args.q_sizes, ','.join(MODES)),
                                     ('C', str(args.n_c), args.q_c, args.modes_c)):
    if tier in tiers:
      for n in map(int, sizes.split(',')):
        for q in map(int, qsizes.split(',')):
          qpath = queries(n, q)
          for mode in modes.split(','):
            for flag in ('off', 'on'):
              add(tier, f'epitopes-{n}', n, q, mode, flag, db(n)[0], qpath, plants=str(qpath) + '.plants.json')
  if 'D' in tiers or 'E' in tiers:
    qpath = queries(args.n_de, args.q_de)
    plants = str(qpath) + '.plants.json'
  if 'D' in tiers:
    for mode in MODES:
      for variant in ('shipped', 'per-length'):
        add('D', f'epitopes-{args.n_de}', args.n_de, args.q_de, mode, 'on', db(args.n_de)[0], qpath, plants, variant)
  if 'E' in tiers:
    for floor in (2, 5, 8):
      for mode in MODES:
        add('E', f'epitopes-{args.n_de}', args.n_de, args.q_de, mode, 'on', db(args.n_de)[0], qpath, plants,
            f'floor={floor}', floor)
      for qs in QUERY_SETS:
        for mode in MODES:
          add('E', f'human/{qs}', 'human', qs, mode, 'on', human_path, query_set_paths[qs], None, f'floor={floor}', floor)

  with open(args.outdir / 'run-info.txt', 'w') as f:
    for k, v in info.items():
      f.write(f'{k}: {v}\n')
    f.write(f'threads: {args.threads}\nrepeats: {args.repeats}\ncap_min: {args.cap_min}\nseed: {args.seed}\n'
            f'smoke: {args.smoke}\nconfigs: {len(configs)}\n')
    for name, (path, digest) in data.items():
      f.write(f'data {name}: {digest}  {path}\n')

  results_path = args.outdir / 'short-target-runs.tsv'
  timed_out = set()
  with open(results_path, 'w', newline='') as out:
    w = csv.DictWriter(out, fieldnames=FIELDS, delimiter='\t', extrasaction='ignore')
    w.writeheader()
    for i, cfg in enumerate(configs, 1):
      base = {k: cfg[k] for k in ('tier', 'dataset', 'n_entries', 'n_queries', 'mode', 'flag', 'variant')}
      kind = (cfg['tier'], cfg['mode'], cfg['flag'], cfg['variant'])
      if cfg['n_entries'] != 'human' and any(k == kind and cfg['n_entries'] >= n and cfg['n_queries'] >= q
                                             for k, n, q in timed_out):
        w.writerow({**base, 'status': 'skipped', 'note': 'a smaller config of this kind hit the cap'})
        continue
      t0 = time.perf_counter()
      try:
        proc = subprocess.run([sys.executable, __file__, '--workdir', str(args.workdir), '--outdir', str(args.outdir),
                               '--child', json.dumps(cfg)], capture_output=True, text=True, timeout=args.cap_min * 60)
        line = next((ln for ln in proc.stdout.splitlines() if ln.startswith('RESULT ')), None)
        if proc.returncode or not line:
          row = {**base, 'status': 'error', 'note': (proc.stderr.strip().splitlines() or ['?'])[-1][:300]}
        else:
          row = {**base, **json.loads(line[len('RESULT '):])}
      except subprocess.TimeoutExpired:
        row = {**base, 'status': 'timeout', 'note': f'over {args.cap_min:g} min'}
        if cfg['n_entries'] != 'human':
          timed_out.add((kind, cfg['n_entries'], cfg['n_queries']))
      w.writerow(row)
      out.flush()
      print(f'[{i}/{len(configs)}] {cfg["tier"]} {cfg["dataset"]} q={cfg["n_queries"]} {cfg["mode"]} '
            f'{cfg["flag"]} {cfg["variant"]}: {row["status"]} ({time.perf_counter() - t0:.1f}s)', flush=True)

  if 'O' in tiers:
    n0, q0 = int(args.sizes_b.split(',')[0]), int(args.q_sizes.split(',')[-1])
    with open(args.outdir / 'short-target-oracle.tsv', 'w', newline='') as out:
      w = csv.DictWriter(out, fieldnames=['mode', 'queries', 'engine_rows', 'oracle_rows', 'false_rows',
                                          'missed_rows'], delimiter='\t')
      w.writeheader()
      for mode in MODES:
        n = args.oracle_n if EDITS[mode][0] == 'mm' else max(1, args.oracle_n // 2)
        r = run_oracle(queries(n0, q0), db(n0)[0], mode, n, args.workdir, args.threads)
        w.writerow(r)
        out.flush()
        print(f'oracle {mode}: {r}', flush=True)

  summarize(results_path, args.outdir)


def summarize(results_path, outdir):
  """Pair flag off/on rows into the table for Dan, checking forward rows are unchanged."""
  with open(results_path) as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
  off = {(r['tier'], r['dataset'], r['n_queries'], r['mode']): r for r in rows
         if r['flag'] == 'off' and r['status'] == 'ok'}
  lines = ['tier\tdataset\tqueries\tmode\toff_s\ton_s\tadded_s\tadded_pct\treversed_s\trows_added\tnoise_rows'
           '\tpeak_mb_on\tplanted_found\tforward_rows_unchanged']
  for r in rows:
    if r['flag'] != 'on' or r['variant'] != 'shipped' or r['status'] != 'ok' or r['tier'] in ('D', 'E'):
      continue
    o = off.get((r['tier'], r['dataset'], r['n_queries'], r['mode']))
    on_s = float(r['total_s'])
    off_s = float(o['total_s']) if o else float('nan')
    planted = f"{r['planted_found']}/{r['planted_in_budget']}" if r['planted_in_budget'] else ''
    lines.append(f"{r['tier']}\t{r['dataset']}\t{r['n_queries']}\t{r['mode']}\t{off_s:.3f}\t{on_s:.3f}"
                 f"\t{on_s - off_s:.3f}\t{100 * (on_s - off_s) / off_s:.1f}\t{float(r['reversed_s']):.3f}"
                 f"\t{r['added_rows']}\t{r['noise_rows']}\t{float(r['peak_mb']):.0f}\t{planted}"
                 f"\t{(o['hit_rows'] == r['forward_rows']) if o else ''}")
  lines.append('')
  lines.append('tier\tdataset\tqueries\tmode\tvariant\ton_s\treversed_s\trows_added\tnoise_rows\tplanted_found')
  for r in rows:
    if r['tier'] in ('D', 'E') and r['status'] == 'ok':
      planted = f"{r['planted_found']}/{r['planted_in_budget']}" if r['planted_in_budget'] else ''
      lines.append(f"{r['tier']}\t{r['dataset']}\t{r['n_queries']}\t{r['mode']}\t{r['variant']}"
                   f"\t{float(r['total_s']):.3f}\t{float(r['reversed_s']):.3f}\t{r['added_rows']}"
                   f"\t{r['noise_rows']}\t{planted}")
  lines.append('')
  lines.append(f'# configs: {len(rows)}, ok: {sum(r["status"] == "ok" for r in rows)}, errors: '
               f'{sum(r["status"] == "error" for r in rows)}, timeouts: {sum(r["status"] == "timeout" for r in rows)}, '
               f'skipped: {sum(r["status"] == "skipped" for r in rows)}')
  for r in rows:
    if r['status'] in ('error', 'timeout'):
      lines.append(f"# {r['status']}: {r['tier']} {r['dataset']} q={r['n_queries']} {r['mode']} {r['flag']} "
                   f"{r['variant']}: {r['note']}")
  (outdir / 'short-target-summary.tsv').write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
  main()
