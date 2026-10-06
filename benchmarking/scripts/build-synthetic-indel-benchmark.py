#!/usr/bin/env python3
"""Fabricate an indel benchmark dataset from the human proteome.

WHY THIS EXISTS
  The scaling figure (PEPMatch vs Brute Force over 100 -> 1M queries) needs a query set
  far larger than any real indel dataset we have, with exact recall ground truth. COSMIC
  is 3,510 queries, all 9-mers, and yields ZERO 2-indel ground truth; CEDAR is 50 unique
  queries and 1-indel only. So we sample real peptide windows out of the proteome and
  inject controlled interior indels: the source window IS the answer, which gives exact
  ground truth, while the residues stay real so aligner composition/E-value statistics
  (and therefore the BLAST/DIAMOND/MMseqs2 comparison) remain representative. Randomly
  generated peptides would give neither.

WHAT IT GUARANTEES
  * Every emitted query has >=1 real, oracle-confirmed match in the proteome.
  * Edits are INTERIOR only. Query-terminal indels are blocked BY DESIGN in both the
    engine (match.rs is_terminal_deletion / the dfs base case) and the committed oracles,
    so planting one would fabricate an expectation no correct tool can ever satisfy.
  * Ground truth is PROTEOME-WIDE: every placement in all proteins, not just the planted
    one -- found with the same sound pigeonhole prefilter the Brute Force baseline uses,
    and spot-audited against an exhaustive no-prefilter scan (--audit-sample).
  * Subsets are strictly nested: the N=100 set is the first 100 queries of the N=1M set,
    so a point on the scaling curve differs from its neighbour only in COUNT.
  * Byte-identical output for a given (seed, proteome, size) -- see the manifest.

GROUND TRUTH IS BUILT FROM THE COMMITTED ORACLE, NEVER A REIMPLEMENTATION. It reuses
benchmarking/methods/brute_force.py wholesale, so the expected rows are by construction
in the exact shape (accession parse, 1-based Index start, exact-matches-dropped) that the
baseline reports and that benchmarking.recall() joins on.

Usage
  build-synthetic-indel-benchmark.py --indels 1 --n-max 1000000 \
      --sizes 100,1000,5000,10000,100000,1000000
  build-synthetic-indel-benchmark.py --indels 1 --n-max 1000 --cost-estimate
  build-synthetic-indel-benchmark.py --indels 1 --print-config
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
import random
import subprocess
import sys
import time
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BENCH = REPO / 'benchmarking'
PROTEOME_DEFAULT = BENCH / 'proteomes' / 'human.fasta'

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
STANDARD = set(AMINO_ACIDS)

# Measured from benchmarking/queries/cedar-indel-test.fasta (69 records). Note the shape:
# ~45% of CEDAR sits at length 25 because those are long candidate peptides, not because
# 25-mers are biologically typical. Kept as the default so the synthetic set mirrors the
# real dataset it stands in for, but --length-dist uniform is available when a flat
# spread over the supported range is wanted instead.
CEDAR_HISTOGRAM = {
  8: 3, 9: 7, 10: 7, 11: 1, 12: 3, 13: 1, 14: 2, 15: 2,
  16: 1, 17: 2, 18: 1, 20: 1, 22: 2, 23: 2, 25: 31, 26: 3,
}

# Per-query ground-truth cost, core-seconds, measured on human.fasta with the same
# prefilter+oracle path this script uses. Cost is driven by the shortest pigeonhole piece
# (len // (indels+1)): a 2-mer seed hits ~70% of the proteome, so short queries with 2
# indels are catastrophically expensive. Used only by --cost-estimate.
COST_TABLE = {
  1: {8: 1.96, 9: 0.80, 10: 0.11, 11: 0.03, 12: 0.41, 13: 0.02, 14: 0.03,
      15: 0.01, 16: 0.02, 17: 0.055, 18: 0.01, 19: 0.01, 20: 0.01, 21: 0.01,
      22: 0.01, 23: 0.02, 24: 0.02, 25: 0.02, 26: 0.02},
  2: {8: 64.7, 9: 15.4, 10: 13.1, 11: 9.6, 12: 1.43, 13: 0.88, 14: 1.70,
      15: 0.22, 16: 0.10, 17: 0.04, 18: 0.011, 19: 0.008, 20: 0.006, 21: 0.005,
      22: 0.004, 23: 0.003, 24: 0.002, 25: 0.001, 26: 0.001},
}
# Observed ratio between the single-core microbenchmark above and real 32-core cluster
# runs (COSMIC: 3,510 9-mers in 310 s wall = 2.79 core-s/query vs 0.80 predicted).
# Scheduling overhead plus the heavy tail of multi-placement queries.
HEAVY_TAIL_FACTOR = 3.3


def load_brute_force():
  """Load benchmarking/methods/brute_force.py by path.

  methods/ has no __init__.py, and this module is the single source of the proteome
  parse, the pigeonhole prefilter, and the committed-oracle call. Re-implementing any of
  it here would risk a ground truth the baseline cannot reproduce -- most insidiously a
  different accession string, which turns every join into a silent 0% recall.
  """
  name = 'bf_for_generator'
  spec = importlib.util.spec_from_file_location(name, BENCH / 'methods' / 'brute_force.py')
  mod = importlib.util.module_from_spec(spec)
  # Register BEFORE exec so the module is importable by name. multiprocessing pickles a
  # worker function by qualified name, so without this Pool(_search_one) dies with
  # "Can't pickle ... import of module 'bf_for_generator' failed".
  sys.modules[name] = mod
  spec.loader.exec_module(mod)
  return mod


def sha256_file(path, _buf=1 << 20):
  h = hashlib.sha256()
  with open(path, 'rb') as f:
    while chunk := f.read(_buf):
      h.update(chunk)
  return h.hexdigest()


def git_blob_hash(path):
  try:
    out = subprocess.run(['git', 'hash-object', str(path)], capture_output=True,
                         text=True, cwd=str(REPO), timeout=15)
    return out.stdout.strip() or None
  except Exception:
    return None


def build_length_sampler(dist, min_len, max_len):
  """Return (draw(rng) -> int, weights dict) for the requested length distribution."""
  if dist == 'uniform':
    weights = {L: 1 for L in range(min_len, max_len + 1)}
  elif dist == 'cedar':
    weights = {L: c for L, c in CEDAR_HISTOGRAM.items() if min_len <= L <= max_len}
  else:
    with open(dist) as f:
      weights = {int(k): float(v) for k, v in json.load(f).items()
                 if min_len <= int(k) <= max_len}
  if not weights:
    raise SystemExit(f'length distribution {dist!r} has no mass inside [{min_len},{max_len}]')

  lengths = sorted(weights)
  total = float(sum(weights[L] for L in lengths))
  cumulative, acc = [], 0.0
  for L in lengths:
    acc += weights[L] / total
    cumulative.append((acc, L))

  def draw(rng):
    u = rng.random()
    for threshold, L in cumulative:
      if u <= threshold:
        return L
    return cumulative[-1][1]

  return draw, weights


def delete_interior(rng, seq, n):
  """Remove n residues from seq, never the first or last position.

  Replaces tests/_with_deletions, which draws from randrange(len(seq)) and so can delete
  a terminal residue. A terminal deletion is barred by the oracle
  (_terminal_deletion_blocked: d == 0 or d == query_len - 1) and by match.rs, so a query
  planted that way would have no findable match and would silently become a decoy with
  no ground-truth row.
  """
  for _ in range(n):
    if len(seq) < 3:
      raise ValueError('sequence too short for an interior deletion')
    d = rng.randrange(1, len(seq) - 1)
    seq = seq[:d] + seq[d + 1:]
  return seq


def insert_interior(rng, seq, n):
  """Insert n residues into seq at interior positions (mirrors tests/_with_insertions,
  whose randrange(1, len(seq)) is already interior-safe)."""
  for _ in range(n):
    i = rng.randrange(1, len(seq))
    seq = seq[:i] + rng.choice(AMINO_ACIDS) + seq[i:]
  return seq


# ---------------------------------------------------------------------------
# Phase 3 worker: exhaustive, NO-prefilter ground truth for one query.
# Reads the proteome from the forked brute_force globals (copy-on-write), same as
# _search_one. Deliberately does not use _pieces/_proteins_containing -- the whole point
# is to check the prefilter did not silently drop a real match.
# ---------------------------------------------------------------------------
_BF = None


def _audit_one(query):
  rows = set()
  for j in range(len(_BF._SEQS)):
    for start, matched in _BF._oracle(query, _BF._SEQS[j], _BF._N):
      if len(matched) != len(query):
        rows.add((query, matched, _BF._ACCS[j], start + 1))
  # Return the query alongside its rows. The caller collects these under imap_unordered,
  # which yields in COMPLETION order -- so results MUST be keyed by this query, never
  # paired positionally against the input list.
  return query, rows


def main():
  p = argparse.ArgumentParser(
    description='Fabricate a synthetic indel benchmark dataset from the human proteome.',
    formatter_class=argparse.RawDescriptionHelpFormatter,
  )
  p.add_argument('--indels', type=int, choices=[1, 2], required=True)
  p.add_argument('--n-max', type=int, default=100000)
  p.add_argument('--sizes', default=None,
                 help='Comma-separated subset sizes; default is every power-of-ten '
                      'ladder point up to --n-max, plus 5000 (the aligner table size).')
  p.add_argument('--min-len', type=int, default=None,
                 help='Default 8 for 1 indel. For 2 indels short queries are ruinously '
                      'expensive (an 8-mer splits into 2,2,4 and a 2-mer seed hits ~70%% '
                      'of the proteome), so consider 12.')
  p.add_argument('--max-len', type=int, default=25)
  p.add_argument('--length-dist', default='cedar',
                 help="'cedar' (default, mirrors the real dataset; ~45%% at length 25), "
                      "'uniform', or a path to a JSON {length: weight} map.")
  p.add_argument('--seed', type=int, default=20260823)
  p.add_argument('--threads', type=int, default=os.cpu_count() or 1)
  p.add_argument('--proteome', default=str(PROTEOME_DEFAULT))
  p.add_argument('--out-prefix', default=None, help='Default synth-<n>indel')
  p.add_argument('--name-prefix', default=None, help='Default synth_<n>indel')
  p.add_argument('--audit-sample', type=int, default=200,
                 help='Queries to re-derive with an exhaustive no-prefilter scan. This '
                      'is the check that the ground truth is real and not just whatever '
                      'the prefilter happened to return. 0 disables (not recommended).')
  p.add_argument('--rejected-cap', type=int, default=10000)
  p.add_argument('--cost-estimate', action='store_true',
                 help='Project ground-truth cost from the length distribution, then exit.')
  p.add_argument('--print-config', action='store_true',
                 help='Print the benchmarking_parameters.json blocks, then exit.')
  p.add_argument('--dry-run', action='store_true',
                 help='Sample and plant, but skip ground truth and writing.')
  args = p.parse_args()

  n = args.indels
  min_len = args.min_len if args.min_len is not None else (8 if n == 1 else 12)
  max_len = args.max_len
  prefix = args.out_prefix or f'synth-{n}indel'
  name_prefix = args.name_prefix or f'synth_{n}indel'

  # A 2-indel query must stay >= 6 after edits (the engine hard-rejects below that), and
  # the prefilter needs len >= n+1 non-empty pieces.
  floor = max(n + 1, 3, 6 if n == 2 else 0)
  if min_len < floor:
    raise SystemExit(f'--min-len {min_len} is below the floor {floor} for {n} indel(s)')
  if max_len < min_len:
    raise SystemExit('--max-len is below --min-len')

  if args.sizes:
    sizes = sorted({int(s) for s in args.sizes.split(',') if s.strip()})
  else:
    sizes = sorted({s for s in [100, 1000, 5000, 10000, 100000, 1000000]
                    if s <= args.n_max} | {args.n_max})
  if sizes[-1] > args.n_max:
    raise SystemExit(f'--sizes asks for {sizes[-1]} but --n-max is {args.n_max}')
  if sizes[-1] < args.n_max:
    # Ground truth is computed for all n_max queries but only the requested subsets are
    # emitted, so this silently throws away compute -- at 1M that is hours. Warn loudly
    # rather than fail: it is occasionally what you want (e.g. deliberately sampling).
    waste = 100.0 * (args.n_max - sizes[-1]) / args.n_max
    print(f'WARNING: largest --sizes entry is {sizes[-1]:,} but --n-max is '
          f'{args.n_max:,}. Ground truth will be computed for all {args.n_max:,} queries '
          f'and {waste:.1f}% of it discarded unwritten. Set --n-max {sizes[-1]} to avoid '
          f'that, or add {args.n_max} to --sizes.', flush=True)

  draw_length, weights = build_length_sampler(args.length_dist, min_len, max_len)

  def label_for(size):
    if size >= 1_000_000 and size % 1_000_000 == 0:
      return f'{size // 1_000_000}m'
    if size >= 1000 and size % 1000 == 0:
      return f'{size // 1000}k'
    return str(size)

  if args.print_config:
    blocks = {}
    for size in sizes:
      blocks[f'{name_prefix}_{label_for(size)}'] = {
        'lengths': sorted(weights),
        'mismatches': 1,
        'indels': n,
        'query': f'queries/{prefix}-{label_for(size)}.fasta',
        'proteome': f'proteomes/{Path(args.proteome).name}',
        'expected': f'expected/{prefix}-{label_for(size)}-expected.tsv',
      }
    print(json.dumps(blocks, indent=2))
    return

  if args.cost_estimate:
    table = COST_TABLE[n]
    total_w = sum(weights.values())
    per_query = sum(weights[L] / total_w * table.get(L, table[max(table)])
                    for L in weights) * HEAVY_TAIL_FACTOR
    print(f'Projected ground-truth cost for {n}-indel, lengths {min_len}-{max_len} '
          f'({args.length_dist}):')
    print(f'  ~{per_query:.2f} core-s/query (incl. {HEAVY_TAIL_FACTOR}x heavy-tail factor)')
    print(f'  {"queries":>10}  {"core-hours":>11}  {"wall @%d cores" % args.threads:>16}')
    for size in sizes:
      ch = per_query * size / 3600
      print(f'  {size:>10,}  {ch:>11.2f}  {ch / args.threads:>13.2f} h')
    print('\nCost is paid again by each benchmark pass (timing, and each memory phase).')
    return

  t_start = time.time()
  print(f'{"=" * 72}\n  SYNTHETIC {n}-INDEL BENCHMARK\n{"=" * 72}')
  print(f'  proteome    {args.proteome}')
  print(f'  n_max       {args.n_max:,}   sizes {sizes}')
  print(f'  lengths     {min_len}-{max_len} ({args.length_dist})')
  print(f'  seed        {args.seed}   threads {args.threads}')
  print()

  # ---- Phase 0: proteome index (reuses the baseline's parse verbatim) -------------
  print('[0/5] building proteome index (reusing methods/brute_force.py)...')
  global _BF
  _BF = load_brute_force()
  indexer = _BF.Benchmarker(
    benchmark='synthetic', query='/dev/null', proteome=args.proteome,
    lengths=[], max_mismatches=0, method_parameters={'threads': args.threads}, indels=n,
  )
  indexer._build_index()
  n_proteins = len(_BF._SEQS)
  n_residues = sum(_BF._LENS)
  print(f'      {n_proteins:,} proteins, {n_residues:,} residues')

  # ---- Phase 1: sample real windows and plant interior indels ---------------------
  print(f'[1/5] sampling {args.n_max:,} queries...')
  rng = random.Random(args.seed)
  tally = Counter()
  accepted = []           # (query, matched_window, accession, start0, match_type)
  seen = set()
  rejected_rows = []
  concat_len = len(_BF._CONCAT)

  while len(accepted) < args.n_max:
    L = draw_length(rng)
    # Strict alternation, not a coin flip: guarantees EVERY nested prefix is exactly
    # 50/50, including the N=100 subset where sampling noise would otherwise show.
    match_type = 'insertion_match' if len(accepted) % 2 == 0 else 'deletion_match'
    # deletion_match  => matched is SHORTER than query  => plant by inserting into a (L-n) window
    # insertion_match => matched is LONGER  than query  => plant by deleting from a (L+n) window
    # (matcher.py: kind = 'd' if len(matched) < len(query) else 'i' -- the label is a
    # property of the length delta, so building it the other way inverts every label.)
    w = L + n if match_type == 'insertion_match' else L - n

    u = rng.randrange(concat_len)
    j = _BF.bisect.bisect_right(_BF._STARTS, u) - 1
    off = u - _BF._STARTS[j]

    # Require a residue of flank on each side: a window flush against a protein terminus
    # can make the planted edit indistinguishable from a protein-boundary case.
    if off < 1 or off + w > _BF._LENS[j] - 1:
      tally['window_crosses_protein_edge'] += 1
      continue
    window = _BF._SEQS[j][off:off + w]
    if not set(window) <= STANDARD:
      tally['nonstandard_residue'] += 1        # X (8,536) and U (36) in human.fasta
      continue

    try:
      if match_type == 'insertion_match':
        query = delete_interior(rng, window, n)
      else:
        query = insert_interior(rng, window, n)
    except ValueError:
      tally['too_short_to_edit'] += 1
      continue

    if len(query) != L:
      tally['length_mismatch_bug'] += 1
      continue
    if query in seen:
      tally['duplicate_query'] += 1
      continue

    # THE GATE. Ask the committed oracle whether the placement we just planted is
    # actually recoverable in its source protein. This is what makes "no terminal edit,
    # no accidental exact match, no index arithmetic slip" a mechanically enforced
    # property instead of a claim in a comment.
    if (off, window) not in _BF._oracle(query, _BF._SEQS[j], n):
      tally['planted_not_recoverable'] += 1
      if len(rejected_rows) < args.rejected_cap:
        rejected_rows.append({
          'Query Sequence': query, 'Window': window, 'Protein ID': _BF._ACCS[j],
          'Index start': off + 1, 'Match Type': match_type,
          'Reason': 'planted placement not returned by the oracle',
        })
      continue

    seen.add(query)
    accepted.append((query, window, _BF._ACCS[j], off, match_type))
    tally['ACCEPTED'] += 1

    if len(accepted) % 25000 == 0:
      print(f'      {len(accepted):,}/{args.n_max:,}  ({time.time() - t_start:.0f}s)')

  draws = sum(tally.values())
  print(f'      accepted {len(accepted):,} from {draws:,} draws')

  if args.dry_run:
    print('\n--dry-run: stopping before ground truth.')
    for reason, count in tally.most_common():
      print(f'  {count:>9,}  ({100 * count / draws:5.2f}%)  {reason}')
    return

  # ---- Phase 2: proteome-wide ground truth ----------------------------------------
  print(f'[2/5] ground truth over all {n_proteins:,} proteins '
        f'({args.threads} threads)...')
  queries_sorted = sorted(seen, key=len)   # shortest (costliest) first: see below
  truth = {}
  t0 = time.time()
  done = 0
  if args.threads > 1:
    with Pool(args.threads) as pool:
      # imap_unordered(chunksize=1) rather than pool.map: cost per query spans ~3 orders
      # of magnitude with length, and map's default chunking hands one worker a block of
      # short queries while the rest idle.
      for rows in pool.imap_unordered(_BF._search_one, queries_sorted, chunksize=1):
        if rows:
          truth.setdefault(rows[0][0], []).extend(rows)
        done += 1
        if done % 1000 == 0:
          rate = done / max(time.time() - t0, 1e-9)
          print(f'      {done:,}/{len(queries_sorted):,}  '
                f'{rate:.1f} q/s  eta {(len(queries_sorted) - done) / max(rate, 1e-9) / 60:.1f} min')
  else:
    for q in queries_sorted:
      rows = _BF._search_one(q)
      if rows:
        truth.setdefault(q, []).extend(rows)
      done += 1
  print(f'      {sum(len(v) for v in truth.values()):,} raw rows in {time.time() - t0:.0f}s')

  # Every planted placement must be in the prefiltered truth. Free, and it covers 100% of
  # queries -- a prefilter that lost planted matches would show up here immediately.
  missing = [(q, w, a, off) for q, w, a, off, _ in accepted
             if (q, w, a, off + 1) not in set(truth.get(q, []))]
  if missing:
    raise SystemExit(f'FATAL: {len(missing)} planted placements absent from ground truth. '
                     f'First: {missing[0]}')
  empty = [q for q, *_ in accepted if not truth.get(q)]
  if empty:
    raise SystemExit(f'FATAL: {len(empty)} queries have no expected rows. First: {empty[0]}')
  print('      all planted placements present; no query left without a match')

  # The audit runs in Phase 4, AFTER the subsets are on disk -- see the note there. It is
  # deliberately NOT a barrier before writing: computing the ground truth costs hours, and
  # a failure in a self-check must never be the reason that work is discarded.
  audit_report = {'sampled': 0, 'ok': True, 'note': 'skipped'}

  # ---- Phase 3: emit nested subsets ------------------------------------------------
  print('[3/5] writing nested subsets...')
  meta = {q: (w, a, off, mt) for q, w, a, off, mt in accepted}
  order = [q for q, *_ in accepted]          # acceptance order == query id order
  queries_dir, expected_dir = BENCH / 'queries', BENCH / 'expected'
  queries_dir.mkdir(exist_ok=True)
  expected_dir.mkdir(exist_ok=True)

  columns = ['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start',
             'Index end', 'Match Type', 'Planted', 'Query Length']
  emitted, prev_set = {}, None
  for size in sizes:
    lab = label_for(size)
    subset = order[:size]
    subset_set = set(subset)
    if prev_set is not None and not prev_set <= subset_set:
      raise SystemExit(f'FATAL: subset {lab} is not a superset of the previous size')
    prev_set = subset_set

    fasta_path = queries_dir / f'{prefix}-{lab}.fasta'
    with open(fasta_path, 'w') as f:
      for qid, q in enumerate(subset, start=1):
        f.write(f'>{qid}\n{q}\n')

    rows = []
    for q in subset:
      window, acc, off, match_type = meta[q]
      planted = (q, window, acc, off + 1)
      for row in dict.fromkeys(truth[q]):            # dedup, preserve discovery order
        rows.append({
          'Query Sequence': row[0],
          'Matched Sequence': row[1],
          'Protein ID': row[2],
          'Index start': row[3],
          'Index end': row[3] - 1 + len(row[1]),
          'Match Type': ('deletion_match' if len(row[1]) < len(row[0])
                         else 'insertion_match'),
          'Planted': 'yes' if row == planted else 'no',
          'Query Length': len(row[0]),
        })
    # Matched Sequence is part of the key because the other three do NOT form a total
    # order: one query can match at the SAME start in the SAME protein both as a
    # deletion and as an insertion (e.g. RHKGEMENALRYS at A0A0A0MRF6:466 -> RHKGEMENALRS
    # and RHKGEMENALRSYS). Both rows are legitimate, so without this the tie fell to
    # set-iteration order and the file was not byte-reproducible -- which would quietly
    # falsify the sha256 the manifest records.
    rows.sort(key=lambda r: (r['Query Sequence'], r['Protein ID'], r['Index start'],
                             r['Matched Sequence']))

    expected_path = expected_dir / f'{prefix}-{lab}-expected.tsv'
    with open(expected_path, 'w', newline='') as f:
      writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t')
      writer.writeheader()
      writer.writerows(rows)

    balance = Counter(meta[q][3] for q in subset)
    emitted[lab] = {
      'size': size,
      'queries_fasta': str(fasta_path.relative_to(BENCH)),
      'expected_tsv': str(expected_path.relative_to(BENCH)),
      'expected_rows': len(rows),
      'insertion_match': balance['insertion_match'],
      'deletion_match': balance['deletion_match'],
    }
    print(f'      {lab:>5}  {size:>9,} queries  {len(rows):>10,} expected rows  '
          f'({balance["insertion_match"]}/{balance["deletion_match"]} ins/del)')

  if rejected_rows:
    rej = expected_dir / f'{prefix}-rejected.tsv'
    with open(rej, 'w', newline='') as f:
      writer = csv.DictWriter(
        f, delimiter='\t',
        fieldnames=['Query Sequence', 'Window', 'Protein ID', 'Index start',
                    'Match Type', 'Reason'])
      writer.writeheader()
      writer.writerows(rejected_rows)

  # ---- Phase 4: exhaustive audit (breaks the circularity) -------------------------
  # Every subset is already on disk, so this gate can never cost the ground-truth compute.
  # It records its verdict in the manifest (audit.ok) and sets the process exit code; the
  # benchmark preflight refuses to run unless audit.ok is true. So a real prefilter fault
  # blocks the benchmarks WITHOUT throwing away hours of correct data.
  if args.audit_sample > 0:
    k = min(args.audit_sample, len(accepted))
    # Stratify across lengths so the audit covers the cheap long queries AND the short
    # ones where the prefilter does the most work.
    by_len = {}
    for q, *_ in accepted:
      by_len.setdefault(len(q), []).append(q)
    audit_rng = random.Random(args.seed + 1)
    pick, lengths_cycle = [], sorted(by_len)
    while len(pick) < k:
      progressed = False
      for L in lengths_cycle:
        if by_len[L] and len(pick) < k:
          pick.append(by_len[L].pop(audit_rng.randrange(len(by_len[L]))))
          progressed = True
      if not progressed:
        break
    print(f'[4/5] exhaustive no-prefilter audit on {len(pick)} queries...')
    t0 = time.time()
    # KEY results by the query each _audit_one returns -- NEVER pair positionally against
    # `pick`. imap_unordered yields in completion order, so the i-th result is not the
    # i-th query; positional pairing falsely failed 195/200 at 1M scale while the ground
    # truth was perfect. dict() over (query, rows) pairs is order-independent.
    if args.threads > 1:
      with Pool(args.threads) as pool:
        audit = dict(pool.imap_unordered(_audit_one, pick, chunksize=1))
    else:
      audit = dict(_audit_one(q) for q in pick)
    mismatches = [q for q in pick if audit.get(q, set()) != set(truth.get(q, []))]
    secs = round(time.time() - t0, 1)
    if mismatches:
      audit_report = {'sampled': len(pick), 'ok': False, 'seconds': secs,
                      'note': f'{len(mismatches)} of {len(pick)} queries disagreed with the '
                              f'exhaustive scan; first {mismatches[0]}'}
      print(f'      !! {len(mismatches)}/{len(pick)} DISAGREED ({secs}s) -- see FATAL below',
            flush=True)
    else:
      audit_report = {'sampled': len(pick), 'ok': True, 'seconds': secs,
                      'note': f'exhaustive scan agreed on all {len(pick)} queries'}
      print(f'      agreed on all {len(pick)} queries ({secs}s)')
  else:
    print('[4/5] audit SKIPPED (--audit-sample 0) -- ground truth is unverified')

  # ---- Phase 5: manifest ------------------------------------------------------------
  print('[5/5] manifest...')
  manifest = {
    'generated_by': str(Path(__file__).relative_to(REPO)),
    'generator_git_blob': git_blob_hash(Path(__file__)),
    'oracle_git_blobs': {
      name: git_blob_hash(REPO / 'pepmatch' / 'tests' / f'{name}.py')
      for name in ('indel_brute_force', 'indel2_brute_force')
    },
    'args': vars(args),
    'seed': args.seed,
    'indels': n,
    'length_range': [min_len, max_len],
    'length_distribution': args.length_dist,
    'length_weights': {str(k): v for k, v in sorted(weights.items())},
    'proteome': {
      'path': str(Path(args.proteome).name),
      'sha256': sha256_file(args.proteome),
      'proteins': n_proteins,
      'residues': n_residues,
    },
    'reconciliation': dict(tally.most_common()),
    'total_draws': draws,
    'observed_length_histogram': dict(sorted(Counter(len(q) for q, *_ in accepted).items())),
    'audit': audit_report,
    'files': {},
  }
  for lab, info in emitted.items():
    info = dict(info)
    info['queries_sha256'] = sha256_file(BENCH / info['queries_fasta'])
    info['expected_sha256'] = sha256_file(BENCH / info['expected_tsv'])
    manifest['files'][lab] = info
  manifest_path = expected_dir / f'{prefix}-manifest.json'
  with open(manifest_path, 'w') as f:
    json.dump(manifest, f, indent=2)

  # ---- Reconciliation: every draw lands in exactly one bucket ----------------------
  print(f'\nRECONCILIATION -- every one of the {draws:,} draws is accounted for:\n')
  for reason, count in tally.most_common():
    print(f'  {count:>9,}  ({100 * count / draws:5.2f}%)  {reason}')
  print(f'  {"-" * 9}')
  print(f'  {draws:>9,}           TOTAL draws')
  print(f'\nmanifest: {manifest_path.relative_to(BENCH)}')
  print(f'done in {(time.time() - t_start) / 60:.1f} min')

  # Gate LAST, so the files and the manifest (recording audit.ok=false) are already
  # persisted. A non-zero exit stops the self-certifying sbatch before the benchmarks run,
  # and the recorded audit status makes the benchmark preflight refuse the dataset too --
  # but the hours of ground truth survive on disk for inspection either way.
  if not audit_report['ok']:
    print(f'\nFATAL: audit failed -- {audit_report["note"]}. Files were written and the '
          f'manifest records audit.ok=false; investigate before trusting these numbers.',
          file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  main()
