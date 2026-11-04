import pytest
import polars as pl
import polars.testing as plt

from pathlib import Path

from pepmatch import Matcher, Preprocessor
from pepmatch._rs import rs_metadata, rs_short_proteins

# Lengths 18, 5, 2, 0, 10. The fourth record has no sequence line, so it is stored as an
# empty protein; it must still take protein number 4 so the fifth keeps number 5.
SHORT_PROTEOME = (
  '>sp|P11111|LONG_HUMAN Long protein OS=Homo sapiens OX=9606 GN=LNG PE=1 SV=1\n'
  'MKVLATGGWQMKRRDEYH\n'
  '>sp|P22222|FIVE_HUMAN Five-mer OS=Homo sapiens OX=9606 GN=FIV PE=1 SV=1\n'
  'RRHFD\n'
  '>sp|P33333|TWO_HUMAN Two-mer OS=Homo sapiens OX=9606 GN=TWO PE=1 SV=1\n'
  'AC\n'
  '>sp|P44444|EMPTY_HUMAN Empty OS=Homo sapiens OX=9606 GN=EMP PE=1 SV=1\n'
  '>sp|P55555|TEN_HUMAN Ten-mer OS=Homo sapiens OX=9606 GN=TEN PE=1 SV=1\n'
  'PQRSTWYVDE\n'
)


@pytest.fixture
def short_proteome(tmp_path) -> Path:
  fasta = tmp_path / 'short.fasta'
  fasta.write_text(SHORT_PROTEOME)
  for k in (2, 8):
    Preprocessor(str(fasta), preprocessed_files_path=str(tmp_path)).preprocess(k)
  return fasta


def index_at(proteome, k):
  return str(proteome.with_name(f'short_{k}mers.pepidx'))


def test_short_proteins_numbering_matches_metadata(short_proteome):
  nums, seqs = rs_short_proteins(index_at(short_proteome, 8), 0, 100)
  assert nums == [1, 2, 3, 4, 5]
  assert seqs == ['MKVLATGGWQMKRRDEYH', 'RRHFD', 'AC', '', 'PQRSTWYVDE']
  ids = dict(zip(*rs_metadata(index_at(short_proteome, 8), nums)[:2]))
  assert [ids[n] for n in nums] == ['P11111', 'P22222', 'P33333', 'P44444', 'P55555']


def test_short_proteins_bounds_are_min_inclusive_max_exclusive(short_proteome):
  index = index_at(short_proteome, 8)
  assert rs_short_proteins(index, 5, 10) == ([2], ['RRHFD'])
  assert rs_short_proteins(index, 5, 11) == ([2, 5], ['RRHFD', 'PQRSTWYVDE'])
  assert rs_short_proteins(index, 1, 5) == ([3], ['AC'])


def test_short_proteins_include_proteins_shorter_than_k(short_proteome):
  # The 5-mer and 2-mer have no 8-mers, so they are absent from the k-mer table, but the
  # reversed search needs them precisely because they are short.
  assert rs_short_proteins(index_at(short_proteome, 8), 2, 8) == ([2, 3], ['RRHFD', 'AC'])


def test_short_proteins_do_not_depend_on_k(short_proteome):
  assert rs_short_proteins(index_at(short_proteome, 2), 0, 100) == rs_short_proteins(index_at(short_proteome, 8), 0, 100)


# ── The short-target search: proteome entries shorter than the query that lie inside it. ──
# Expected rows below are worked out by hand from the query's side, never by running the
# engine: Matched Sequence is the whole entry, Index start/end span all of it, Query
# start/end say where it sits in the query, and edit positions count along the query.

QUERY = 'MKTAYIAKQRQISFVK'   # 16 residues: M1 K2 T3 A4 Y5 I6 A7 K8 Q9 R10 Q11 I12 S13 F14 V15 K16


def write_entries(path, entries):
  path.write_text(''.join(
    f'>sp|{name}|{name} {name} protein OS=Test OX=1 GN={name} PE=1 SV=1\n{seq}\n'
    for name, seq in entries.items()
  ))
  return path


def search(tmp_path, entries, query=None, short_targets=True, **kwargs):
  """Search `query` (default [QUERY]) against a FASTA of `entries` (name -> sequence)."""
  fasta = write_entries(tmp_path / 'db.fasta', entries)
  return Matcher(
    query=query or [QUERY], proteome_file=str(fasta), preprocessed_files_path=str(tmp_path),
    output_format='dataframe', short_targets=short_targets, **kwargs,
  ).match()


def hits(df):
  """(Protein ID, Matched, edits, positions, Index start, Index end, Query start, Query end)."""
  edit_col, pos_col = df.columns[8], df.columns[9]
  h = df.filter(pl.col('Matched Sequence').is_not_null())
  return sorted(zip(
    h['Protein ID'], h['Matched Sequence'], h[edit_col], h[pos_col],
    h['Index start'], h['Index end'], h['Query start'], h['Query end'],
  ))


@pytest.mark.parametrize('target, q_start, q_end', [
  (QUERY[0:10], 1, 10),    # prefix
  (QUERY[2:12], 3, 12),    # interior: holds only one of the query's 8-mer tiles
  (QUERY[6:16], 7, 16),    # suffix
  (QUERY[4:12], 5, 12),    # 8-mer
  (QUERY[3:8], 4, 8),      # 5-mer
  (QUERY[6:8], 7, 8),      # 2-mer, the exact-search floor
])
def test_exact_entry_inside_query(tmp_path, target, q_start, q_end):
  df = search(tmp_path, {'T': target})
  assert hits(df) == [('T.1', target, 0, '[]', 1, len(target), q_start, q_end)]


@pytest.mark.parametrize('target, mismatches, positions, q_start, q_end', [
  ('KTAYIAWQRQ', 1, '[8]', 2, 11),       # K8 -> W, entry starts at query position 2
  ('YIWKQRQI', 1, '[7]', 5, 12),         # A7 -> W, entry starts at query position 5
  ('WTAYIAKQRQ', 1, '[2]', 2, 11),       # mismatch on the entry's own first residue
  ('TAWIAKQWQI', 2, '[5, 10]', 3, 12),   # Y5 -> W and R10 -> W
])
def test_mismatched_entry_positions_count_along_the_query(tmp_path, target, mismatches, positions, q_start, q_end):
  df = search(tmp_path, {'T': target}, max_mismatches=mismatches)
  assert hits(df) == [('T.1', target, mismatches, positions, 1, len(target), q_start, q_end)]


@pytest.mark.parametrize('target, indels, label, q_start, q_end', [
  # The entry carries a W the query lacks: an insertion before query position 9 (Q9). Read
  # from the entry's side it would be 'd: W[8]'.
  ('KTAYIAKWQRQISFV', 1, 'i: W[9]', 2, 15),
  ('KTAYIKQRQIS', 1, 'd: A[7]', 2, 13),             # A7 missing from the entry
  ('TAYIWAKQRQWISF', 2, 'i: W[7], i: W[12]', 3, 14),
  ('KTAYAKQRISFV', 2, 'd: I[6], d: Q[11]', 2, 15),
])
def test_indel_entry_labels_count_along_the_query(tmp_path, target, indels, label, q_start, q_end):
  df = search(tmp_path, {'T': target}, max_indels=indels)
  assert hits(df) == [('T.1', target, indels, label, 1, len(target), q_start, q_end)]


def test_entry_twice_in_query_reports_both_places(tmp_path):
  df = search(tmp_path, {'T': 'AKQRQISF'}, query=['AKQRQISFAKQRQISF'])
  assert hits(df) == [
    ('T.1', 'AKQRQISF', 0, '[]', 1, 8, 1, 8),
    ('T.1', 'AKQRQISF', 0, '[]', 1, 8, 9, 16),
  ]


def test_forward_row_wins_and_other_placements_are_kept(tmp_path):
  # A9 inside A10 with one indel has six placements (indel2_brute_force agrees): the whole
  # query with one A missing from the entry, two exact 9-residue windows, and three
  # 8-residue windows the entry covers with one extra A. The whole-query one is the forward
  # search's alignment and is kept once, in its forward row; the other five are new. Every
  # valid placement is reported, as everywhere else in PEPMatch.
  df = search(tmp_path, {'T': 'A' * 9}, query=['A' * 10], max_indels=1)
  assert hits(df) == [
    ('T.1', 'A' * 9, 0, '[]', 1, 9, 1, 9),
    ('T.1', 'A' * 9, 0, '[]', 1, 9, 2, 10),
    ('T.1', 'A' * 9, 1, 'd: A[2,9]', 1, 9, 1, 10),
    ('T.1', 'A' * 9, 1, 'i: A[2,8]', 1, 9, 1, 8),
    ('T.1', 'A' * 9, 1, 'i: A[3,9]', 1, 9, 2, 9),
    ('T.1', 'A' * 9, 1, 'i: A[4,10]', 1, 9, 3, 10),
  ]
  assert df.filter((pl.col('Query start') == 1) & (pl.col('Query end') == 10)).height == 1


def test_entry_shorter_only_by_its_indel_stays_a_forward_row(tmp_path):
  target = QUERY[:7] + QUERY[8:]   # K8 missing: 15 residues, found by the forward search too
  df = search(tmp_path, {'T': target}, max_indels=1)
  assert hits(df) == [('T.1', target, 1, 'd: K[8]', 1, 15, 1, 16)]


def test_entry_as_long_as_query_is_forward_only(tmp_path):
  assert hits(search(tmp_path, {'T': QUERY})) == [('T.1', QUERY, 0, '[]', 1, 16, 1, 16)]


def test_query_inside_long_protein_is_unchanged(tmp_path):
  df = search(tmp_path, {'L': 'GGGG' + QUERY + 'GGGG'})
  assert hits(df) == [('L.1', QUERY, 0, '[]', 5, 20, 1, 16)]


@pytest.mark.parametrize('target, kwargs', [
  ('WWWW' + QUERY[:6], {}),                 # overhang: neither sequence contains the other
  (QUERY[2:12] + 'W', {'max_indels': 1}),   # extra residue on the entry's end: an overhang of one
])
def test_overlap_that_is_not_containment_is_not_reported(tmp_path, target, kwargs):
  df = search(tmp_path, {'T': target}, **kwargs)
  assert df.height == 1 and df['Matched Sequence'][0] is None


def test_query_with_no_hits_keeps_one_miss_row(tmp_path):
  df = search(tmp_path, {'T': 'WWWWWWWW'})
  assert df.height == 1
  assert df.select('Query start', 'Query end').row(0) == (None, None)


def test_miss_row_replaced_by_short_target_hit(tmp_path):
  df = search(tmp_path, {'T': QUERY[2:12]})
  assert df.height == 1 and df['Matched Sequence'][0] == QUERY[2:12]


def test_duplicate_entries_each_get_a_row(tmp_path):
  df = search(tmp_path, {'T1': QUERY[2:12], 'T2': QUERY[2:12]})
  assert [h[0] for h in hits(df)] == ['T1.1', 'T2.1']


def test_duplicate_query_ids_keep_their_own_sequences(tmp_path):
  queries = tmp_path / 'queries.fasta'
  queries.write_text(f'>q\n{QUERY}\n>q\nGGGGGGGGWWWWWWWW\n')
  df = search(tmp_path, {'T': QUERY[2:12]}, query=str(queries))
  rows = df.select('Query ID', 'Query Sequence', 'Matched Sequence').rows()
  assert ('q', QUERY, QUERY[2:12]) in rows
  assert ('q', 'GGGGGGGGWWWWWWWW', None) in rows
  assert len(rows) == 2


def test_empty_and_one_residue_entries_keep_numbering(tmp_path):
  df = search(tmp_path, {'E0': '', 'E1': 'K', 'T': QUERY[2:12]})
  h = df.filter(pl.col('Matched Sequence').is_not_null())
  assert h.select('Protein ID', 'Gene', 'Protein Name').rows() == [('T.1', 'T', 'T protein')]


def test_sequence_version_off(tmp_path):
  df = search(tmp_path, {'T': QUERY[2:12]}, sequence_version=False)
  assert hits(df)[0][0] == 'T'


def test_discontinuous_queries_get_null_query_span(tmp_path):
  df = search(tmp_path, {'T': QUERY[2:12], 'L': 'GGGG' + QUERY + 'GGGG'}, query=[QUERY, 'K6, T7'])
  disc = df.filter(pl.col('Query Sequence') == 'K6, T7')
  assert disc.height >= 1
  assert disc['Query start'].is_null().all() and disc['Query end'].is_null().all()
  assert ('T.1', QUERY[2:12]) in df.select('Protein ID', 'Matched Sequence').rows()


@pytest.mark.parametrize('fmt', ['csv', 'tsv', 'json', 'xlsx'])
def test_file_output_carries_query_span(tmp_path, fmt):
  fasta = write_entries(tmp_path / 'db.fasta', {'T': QUERY[2:12]})
  out = tmp_path / 'results'
  Matcher(query=[QUERY], proteome_file=str(fasta), preprocessed_files_path=str(tmp_path),
          output_format=fmt, output_name=str(out), short_targets=True).match()
  path = Path(f'{out}.{fmt}')
  assert path.exists()
  if fmt in ('csv', 'tsv'):
    df = pl.read_csv(path, separator='\t' if fmt == 'tsv' else ',')
    assert df.select('Query start', 'Query end').row(0) == (3, 12)


def test_no_short_entries_only_adds_the_span(tmp_path):
  entries = {'L1': 'GGGG' + QUERY + 'GGGG', 'L2': 'WWWWWWWWWWWWWWWWWWWWWWWW'}
  off = search(tmp_path, entries, short_targets=False)
  on = search(tmp_path, entries)
  plt.assert_frame_equal(on.drop('Query start', 'Query end'), off)


@pytest.mark.parametrize('kwargs', [{'max_mismatches': 2}, {'max_indels': 2}])
def test_entries_below_the_seed_guarantee_are_searched_with_a_warning(tmp_path, capsys, kwargs):
  # With two edits, 2-mer seeds guarantee recall only from 2 * (2 + 1) = 6 residues. As for
  # a short forward query, a 5-residue entry is still searched and the user warned: an exact
  # placement is always found, only some placements with edits may be missed.
  # (With two indels S6 also aligns to query 4-11 with QR missing from it, a second valid
  # placement, so only which entries are found is checked here.)
  df = search(tmp_path, {'F5': QUERY[3:8], 'S6': QUERY[3:9]}, **kwargs)
  assert {h[0] for h in hits(df)} == {'F5.1', 'S6.1'}
  out = capsys.readouterr().out
  assert 'complete recall is not guaranteed for proteome entries shorter than 6 residues (lengths: [5])' in out


def test_entries_longer_than_index_k_limit(tmp_path):
  # The k-mer table holds at most 16 residues per k-mer, so a 25-residue entry is searched
  # at k=16, not at its own length. The forward k is explicit because auto-k on these
  # 40-mers would exceed 16 too.
  query = 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRV'
  target = query[5:30]
  df = search(tmp_path, {'T': target}, query=[query], k=5)
  assert hits(df) == [('T.1', target, 0, '[]', 1, 25, 6, 30)]


@pytest.mark.parametrize('kwargs', [{'best_match': True}, {'counts_only': True}])
def test_short_targets_rejects_unsupported_modes(tmp_path, kwargs):
  with pytest.raises(ValueError, match='short_targets'):
    search(tmp_path, {'T': QUERY[2:12]}, **kwargs)


@pytest.mark.parametrize('query', [' ' + QUERY, QUERY + '\n', QUERY[:8] + ' ' + QUERY[8:]])
def test_short_targets_rejects_whitespace_in_a_query(tmp_path, query):
  # The queries are indexed through a FASTA file whose reader trims each line, so a leading
  # space would move every Query start/end by one.
  with pytest.raises(ValueError, match='whitespace'):
    search(tmp_path, {'T': QUERY[2:12]}, query=[query])
  assert search(tmp_path, {'T': QUERY[2:12]}, query=[query], short_targets=False).height >= 1


INVARIANT_ENTRIES = {
  'L1': 'GGGG' + QUERY + 'GGGG',
  'L2': 'SHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ',
  'T1': QUERY[2:12],
  'T2': 'TAYIWAKQRQ',
  'T3': 'LEERLGLIEV',
  'T4': 'AKQRQ',
}
INVARIANT_QUERIES = [QUERY, 'RQLEERLGLIEVQAPI', 'GGGGGGGGGGGG', 'TAYIAKQRQISF']


@pytest.mark.parametrize('kwargs', [{}, {'max_mismatches': 1}, {'max_indels': 1}, {'max_indels': 2}])
def test_forward_rows_are_unchanged(tmp_path, kwargs):
  off = search(tmp_path, INVARIANT_ENTRIES, query=INVARIANT_QUERIES, short_targets=False, **kwargs)
  on = search(tmp_path, INVARIANT_ENTRIES, query=INVARIANT_QUERIES, **kwargs)
  hit = pl.col('Matched Sequence').is_not_null()

  # Forward hits come first and untouched. Compared as sorted rows: forward exact search
  # returns one query's several hits in no fixed order, so two runs can list them differently.
  off_hits = off.filter(hit)
  forward = on.filter(hit).head(off_hits.height).drop('Query start', 'Query end')
  plt.assert_frame_equal(forward.sort(forward.columns), off_hits.sort(off_hits.columns))

  # Every added row is a whole entry, placed inside a longer query.
  added = on.filter(hit).slice(off_hits.height)
  assert added.height > 0
  assert (added['Index start'] == 1).all()
  assert (added['Index end'] == added['Matched Sequence'].str.len_chars()).all()
  assert (added['Matched Sequence'].str.len_chars() < added['Query Sequence'].str.len_chars()).all()

  # A query keeps its miss row unless the short-target search found it something.
  for qid, qseq in off.filter(~hit).select('Query ID', 'Query Sequence').rows():
    mine = on.filter((pl.col('Query ID') == qid) & (pl.col('Query Sequence') == qseq))
    assert mine.height == 1 or mine.filter(hit).height == mine.height

  # The added rows are sorted, so they come out identically on every run.
  again = search(tmp_path, INVARIANT_ENTRIES, query=INVARIANT_QUERIES, **kwargs).filter(hit)
  plt.assert_frame_equal(again.slice(off_hits.height), added)
