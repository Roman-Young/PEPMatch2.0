"""rs_short_target_match against the composite it replaced: read the entries out with
rs_short_proteins, search them through rs_match / rs_indel_match against the queries index,
and drop the misses. Same entries, same per-peptide search, so the hits must be the same rows;
only where the entries and the misses live differs. Exact search emits an entry's hits in hash
order, so rows are compared as sorted lists."""
import random
import pytest
from pepmatch._rs import (rs_preprocess, rs_short_proteins, rs_short_lengths, rs_short_target_match,
                          rs_match, rs_indel_match)

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
REPEATS = ['AAAAAA', 'EKEKEK', 'QQQQ', 'LSLSLS', 'GGGGG']
MODES = {'exact': (False, 0), 'mm1': (False, 1), 'mm2': (False, 2), 'indel1': (True, 1), 'indel2': (True, 2)}


def _seq(rng, n):
  return ''.join(rng.choice(AMINO_ACIDS) for _ in range(n))


def _mutate(rng, s, indels, n):
  s = list(s)
  for _ in range(n):
    i = rng.randrange(len(s))
    if not indels:
      s[i] = rng.choice(AMINO_ACIDS)
    elif rng.random() < 0.5 and len(s) > 2:
      del s[i]
    else:
      s.insert(i, rng.choice(AMINO_ACIDS))
  return ''.join(s)


def _corpus(rng, indels, edits):
  """Queries of 6-18 aa (a quarter from repeat tracts) and entries of 1-20 aa: random ones,
  windows cut from queries, the same windows with up to `edits` edits, and repeats of both."""
  queries = []
  for _ in range(40):
    n = rng.randint(6, 18)
    if rng.random() < 0.25:
      tract = ''.join(rng.sample(REPEATS, len(REPEATS)))
      start = rng.randrange(len(tract) - n + 1)
      queries.append(tract[start:start + n])
    else:
      queries.append(_seq(rng, n))
  entries = [_seq(rng, rng.randint(1, 20)) for _ in range(60)]
  for _ in range(120):
    q = rng.choice(queries)
    n = rng.randint(2, len(q))
    start = rng.randrange(len(q) - n + 1)
    window = q[start:start + n]
    entries.append(window if rng.random() < 0.4 else _mutate(rng, window, indels, rng.randint(1, max(1, edits))))
  entries += rng.sample(entries, 20)   # duplicate entries: one row set per protein number
  return queries, [e for e in entries if e]


def _write(path, seqs):
  with open(path, 'w') as f:
    for i, s in enumerate(seqs):
      f.write(f'>sp|E{i:05d}|E{i:05d} entry OS=Test OX=1 GN=E{i} PE=1 SV=1\n{s}\n')


def _reference(proteome_idx, queries_idx, lo, hi, k, edits, indels):
  nums, seqs = rs_short_proteins(proteome_idx, lo, hi)
  entries = [(str(n), s) for n, s in zip(nums, seqs)]
  if indels:
    cols = rs_indel_match(queries_idx, entries, edits)
  else:
    cols = rs_match(queries_idx, entries, k, edits)
  ids, seq, win, qnum, n, _, qs, qe = cols
  return sorted((int(i), s, w, q, e, a, b) for i, s, w, q, e, a, b in zip(ids, seq, win, qnum, n, qs, qe) if w is not None)


@pytest.mark.parametrize('mode', MODES)
@pytest.mark.parametrize('seed', range(4))
def test_fused_search_returns_the_composite_hits(tmp_path, mode, seed):
  indels, edits = MODES[mode]
  rng = random.Random(f'{mode}-{seed}')
  queries, entries = _corpus(rng, indels, edits)
  longest = max(map(len, queries))
  proteome, queries_fa = tmp_path / 'proteome.fasta', tmp_path / 'queries.fasta'
  _write(proteome, entries)
  _write(queries_fa, queries)
  proteome_idx = str(tmp_path / 'proteome_2mers.pepidx')
  rs_preprocess(str(proteome), 2, proteome_idx)

  split = 3 * (edits + 1)
  # The search's own ranges (whole span; indel short / rest split) plus ranges that cut through
  # it, an empty one, and floors above the entries' minimum.
  ranges = [(2, longest), (2, min(split, longest)), (split, longest), (5, 9), (1, 3), (10, 10), (8, longest)]
  total = 0
  for k in (2, 3, 4):
    queries_idx = str(tmp_path / f'queries_{k}mers.pepidx')
    rs_preprocess(str(queries_fa), k, queries_idx)
    for lo, hi in ranges:
      fused = rs_short_target_match(proteome_idx, queries_idx, lo, hi, k, edits, indels)
      got = sorted(zip(*fused))
      assert got == _reference(proteome_idx, queries_idx, lo, hi, k, edits, indels), (k, lo, hi)
      total += len(got)
  assert total > 0   # the corpus plants windows of the queries, so a mode that finds nothing is broken


def test_fused_search_reports_rows_in_entry_order(tmp_path):
  # Python sorts the hits anyway, but the contract is protein order, as rs_short_proteins gives.
  queries = ['MKVLATGGWQMK', 'RRHFDPQRST']
  entries = ['WQMK', 'XXXX', 'RRHF', 'LATG', 'PQRS']
  _write(tmp_path / 'p.fasta', entries)
  _write(tmp_path / 'q.fasta', queries)
  rs_preprocess(str(tmp_path / 'p.fasta'), 2, str(tmp_path / 'p.pepidx'))
  rs_preprocess(str(tmp_path / 'q.fasta'), 2, str(tmp_path / 'q.pepidx'))
  nums, seqs, wins, qnums, n, qs, qe = rs_short_target_match(
    str(tmp_path / 'p.pepidx'), str(tmp_path / 'q.pepidx'), 2, 12, 2, 0, False)
  assert nums == [1, 3, 4, 5]
  assert seqs == wins == ['WQMK', 'RRHF', 'LATG', 'PQRS']
  assert qnums == [1, 2, 1, 2]
  assert list(zip(qs, qe)) == [(9, 12), (1, 4), (4, 7), (6, 9)]
  assert n == [0, 0, 0, 0]


def test_lengths_and_fused_search_select_the_same_entries(tmp_path):
  rng = random.Random(5)
  entries = [_seq(rng, rng.randint(1, 20)) for _ in range(200)]
  _write(tmp_path / 'p.fasta', entries)
  rs_preprocess(str(tmp_path / 'p.fasta'), 2, str(tmp_path / 'p.pepidx'))
  for lo, hi in [(2, 21), (5, 9), (9, 9), (0, 100)]:
    expected = {}
    for s in entries:
      if lo <= len(s) < hi:
        expected[len(s)] = expected.get(len(s), 0) + 1
    assert rs_short_lengths(str(tmp_path / 'p.pepidx'), lo, hi) == sorted(expected.items())
