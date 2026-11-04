"""Randomized checks of the short-target search against an independent definition.

With short_targets on, a hit is either the forward search's (the whole query on a window of
an entry, unchanged from today) or a containment (the whole entry on a window of a longer
query). Both are enumerated directly here -- Hamming distance over windows, or the
indel2_brute_force oracle with the roles swapped -- never by running the engine.

The check_* functions return a list of problems so scripts/validate_short_targets.py can run
them over thousands of seeds; the tests below run a fixed handful."""
import random
import re
import pytest
from pepmatch import Matcher
from pepmatch.matcher import MIN_TARGET_LEN
from indel2_brute_force import brute_force_search

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
MODES = {
  'exact': ('mm', 0),
  'mm1': ('mm', 1),
  'mm2': ('mm', 2),
  'indel1': ('indel', 1),
  'indel2': ('indel', 2),
}
REPEAT_BLOCKS = ['AAAAAA', 'EKEKEK', 'QQQQ', 'LSLSLS', 'GGGGG']


def _random_sequence(rng, length):
  return ''.join(rng.choice(AMINO_ACIDS) for _ in range(length))


def _query(rng, longest):
  # Repeat tracts are where seeds are least unique and indel placement is ambiguous.
  length = rng.randint(8, longest)
  if rng.random() < 0.25:
    seq = ''.join(rng.sample(REPEAT_BLOCKS, len(REPEAT_BLOCKS)))
    start = rng.randrange(len(seq) - length + 1)
    return seq[start:start + length]
  return _random_sequence(rng, length)


def _edited(rng, seq, kind, edits):
  """`seq` with `edits` changes of the search's kind: substitutions for a mismatch search,
  for an indel search interior insertions or interior deletions, never both."""
  if kind == 'mm':
    for _ in range(edits):
      i = rng.randrange(len(seq))
      seq = seq[:i] + rng.choice(AMINO_ACIDS) + seq[i + 1:]
  elif rng.random() < 0.5:
    for _ in range(edits):
      if len(seq) > 2:
        i = rng.randrange(1, len(seq) - 1)
        seq = seq[:i] + seq[i + 1:]
  else:
    for _ in range(edits):
      i = rng.randrange(1, len(seq))
      seq = seq[:i] + rng.choice(AMINO_ACIDS) + seq[i:]
  return seq


def make_case(rng, mode):
  """Queries, and entries that hit them both ways: long ones holding an edited query, short
  ones cut from an edited query window, plus unrelated ones and a query copied whole."""
  kind, edits = MODES[mode]
  # Forward exact search takes k from the shortest query and the index caps k at 16. Rows
  # are keyed on the query sequence, so the queries are kept distinct.
  queries = []
  while len(queries) < 5:
    q = _query(rng, 16 if mode == 'exact' else 22)
    if q not in queries:
      queries.append(q)
  entries = {}
  for i in range(3):
    core = _edited(rng, rng.choice(queries), kind, rng.randint(0, edits))
    entries[f'L{i}'] = _random_sequence(rng, rng.randint(2, 8)) + core + _random_sequence(rng, rng.randint(2, 8))
  for i in range(10):
    q = rng.choice(queries)
    length = rng.randint(2, len(q) - 1)
    start = rng.randrange(len(q) - length + 1)
    entries[f'S{i}'] = _edited(rng, q[start:start + length], kind, rng.randint(0, edits))
  entries['U0'] = _random_sequence(rng, rng.randint(4, 12))
  entries['E0'] = rng.choice(queries)
  return queries, entries


def write_fasta(path, entries):
  path.write_text(''.join(f'>{name}\n{seq}\n' for name, seq in entries.items()))
  return path


def _windows(short, long_, kind, edits):
  """(0-based start, window) for every window of `long_` that `short` aligns to in full."""
  if kind == 'mm':
    n = len(short)
    return {
      (s, long_[s:s + n]) for s in range(len(long_) - n + 1)
      if sum(a != b for a, b in zip(short, long_[s:s + n])) <= edits
    }
  return brute_force_search(short, long_, edits)


def expected_rows(queries, entries, mode):
  """{row: guaranteed} with row = (query, protein ID, matched, index start, index end,
  query start, query end). guaranteed marks rows the engine must find: the searched peptide
  -- the query going forward, the entry in reverse -- holds edits + 1 disjoint 2-mer seeds."""
  kind, edits = MODES[mode]
  floor = 2 * (edits + 1)
  rows = {}
  for q in queries:
    for name, t in entries.items():
      pid = f'{name}.1'
      for s, w in _windows(q, t, kind, edits):
        rows[(q, pid, w, s + 1, s + len(w), 1, len(q))] = len(q) >= floor
      if MIN_TARGET_LEN <= len(t) < len(q):
        for s, w in _windows(t, q, kind, edits):
          row = (q, pid, t, 1, len(t), s + 1, s + len(w))
          rows[row] = rows.get(row, False) or len(t) >= floor
  return rows


def _label_rebuilds(window, offset, label, matched):
  """Whether an Indel Positions label, read along the full query, turns the query window it
  covers into `matched`. Entries are applied together ('d: C[3], d: E[5]'), except when each
  entry alone accounts for every edit: then they are alternative placements of one chunk in
  a periodic repeat ('d: LA[8], d: AL[9]'), and each must work on its own."""
  edits = re.findall(r'([di]): ([A-Z]+)\[(\d+)(?:,(\d+))?\]', label)
  n = abs(len(matched) - len(window))
  if len(edits) > 1 and all(len(residues) == n for _, residues, _, _ in edits):
    return all(matched in _apply(window, offset, [e]) for e in edits)
  return matched in _apply(window, offset, edits)


def _apply(window, offset, edits):
  """Every sequence the edits produce on the window, trying each range end."""
  options = [[int(lo), int(hi or lo)] for _, _, lo, hi in edits]
  results = set()
  for choice in range(1 << len(edits)):
    seq = window
    shift = 0
    for e, (kind, residues, _, _) in enumerate(edits):
      pos = options[e][(choice >> e) & 1] - offset - 1 + shift   # 0-based in the edited window
      if kind == 'd':
        if seq[pos:pos + len(residues)] != residues:
          break
        seq = seq[:pos] + seq[pos + len(residues):]
        shift -= len(residues)
      else:
        seq = seq[:pos] + residues + seq[pos:]
        shift += len(residues)
    else:
      results.add(seq)
  return results


def check_case(tmp_path, seed, mode):
  """Problems found for one random case; empty when the engine agrees with the definition."""
  rng = random.Random(seed * 31 + list(MODES).index(mode))
  queries, entries = make_case(rng, mode)
  kind, edits = MODES[mode]
  proteome = write_fasta(tmp_path / f'entries_{seed}_{mode}.fasta', entries)
  df = Matcher(
    query=queries, proteome_file=str(proteome), preprocessed_files_path=str(tmp_path),
    output_format='dataframe', short_targets=True,
    **({'max_mismatches': edits} if kind == 'mm' else {'max_indels': edits}),
  ).match()

  problems = []
  edit_col, pos_col = df.columns[8], df.columns[9]
  found = []
  for row in df.iter_rows(named=True):
    if row['Matched Sequence'] is None:
      continue
    q, m, qs, qe = row['Query Sequence'], row['Matched Sequence'], row['Query start'], row['Query end']
    found.append((q, row['Protein ID'], m, row['Index start'], row['Index end'], qs, qe))
    # Positions must read along the full query: rebuilding the matched sequence from the
    # query window and the reported edits catches a label from the wrong side or offset.
    window = q[qs - 1:qe]
    if kind == 'mm':
      want = '[' + ', '.join(str(qs + i) for i, (a, b) in enumerate(zip(window, m)) if a != b) + ']'
      if row[pos_col] != want:
        problems.append(f'{mode} seed {seed}: positions {row[pos_col]} != {want} for {q} / {m} at {qs}-{qe}')
    elif not _label_rebuilds(window, qs - 1, row[pos_col], m):
      problems.append(f'{mode} seed {seed}: label {row[pos_col]} does not turn {window} (query {qs}-{qe}) into {m}')
    if row[edit_col] != (sum(a != b for a, b in zip(window, m)) if kind == 'mm' else abs(len(m) - len(window))):
      problems.append(f'{mode} seed {seed}: edit count {row[edit_col]} wrong for {q} / {m} at {qs}-{qe}')

  if len(found) != len(set(found)):
    problems.append(f'{mode} seed {seed}: duplicate rows')
  expected = expected_rows(queries, entries, mode)
  for row in set(found) - set(expected):
    problems.append(f'{mode} seed {seed}: reported but not a valid alignment: {row}')
  for row in {r for r, guaranteed in expected.items() if guaranteed} - set(found):
    problems.append(f'{mode} seed {seed}: missed: {row}')
  return problems


def check_symmetry(tmp_path, seed, mode):
  """Searching query A against short entries B gives the same containment alignments as
  searching each B against A: Bjoern's target of identical results whichever side a
  sequence is on. Both directions pick k by the same rule from the same B's, so this holds
  even for B's too short for the seed guarantee."""
  rng = random.Random(10_000 + seed * 31 + list(MODES).index(mode))
  kind, edits = MODES[mode]
  a = _query(rng, 16 if mode == 'exact' else 22)
  bs = {}
  for i in range(8):
    length = rng.randint(MIN_TARGET_LEN, len(a) - 1)
    start = rng.randrange(len(a) - length + 1)
    b = _edited(rng, a[start:start + length], kind, rng.randint(0, edits))
    if MIN_TARGET_LEN <= len(b) < len(a):
      bs[f'B{i}'] = b
  if not bs:
    return []
  budget = {'max_mismatches': edits} if kind == 'mm' else {'max_indels': edits}

  forward_db = write_fasta(tmp_path / f'bs_{seed}_{mode}.fasta', bs)
  one = Matcher(query=[a], proteome_file=str(forward_db), preprocessed_files_path=str(tmp_path),
                output_format='dataframe', short_targets=True, **budget).match()
  edit_col = one.columns[8]
  a_side = {
    (r['Protein ID'][:-2], r['Query start'], r['Query end'], r[edit_col])
    for r in one.iter_rows(named=True)
    if r['Matched Sequence'] is not None and r['Matched Sequence'] == bs[r['Protein ID'][:-2]]
    and r['Index start'] == 1 and r['Index end'] == len(r['Matched Sequence'])
  }

  reverse_db = write_fasta(tmp_path / f'a_{seed}_{mode}.fasta', {'A': a})
  names = list(bs)
  two = Matcher(query=[bs[n] for n in names], proteome_file=str(reverse_db), preprocessed_files_path=str(tmp_path),
                output_format='dataframe', **budget).match()
  b_side = {
    (names[int(r['Query ID']) - 1], r['Index start'], r['Index end'], r[edit_col])
    for r in two.iter_rows(named=True) if r['Matched Sequence'] is not None
  }
  if a_side != b_side:
    return [f'{mode} seed {seed}: A={a} B={bs}: only searching A {sorted(a_side - b_side)}, '
            f'only searching the Bs {sorted(b_side - a_side)}']
  return []


@pytest.mark.parametrize('mode', list(MODES))
@pytest.mark.parametrize('seed', range(15))
def test_short_targets_match_the_definition(tmp_path, seed, mode):
  assert check_case(tmp_path, seed, mode) == []


@pytest.mark.parametrize('mode', list(MODES))
@pytest.mark.parametrize('seed', range(10))
def test_containment_is_the_same_from_either_side(tmp_path, seed, mode):
  assert check_symmetry(tmp_path, seed, mode) == []
