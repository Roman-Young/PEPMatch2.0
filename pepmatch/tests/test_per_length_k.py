"""Every query is searched at its own k, the largest its length allows for the edit budget,
whatever else is in the batch. These tests pin the rule, that a mixed batch returns exactly
what its queries return alone and in the order given, and that best match's budgets keep
the recall guarantee."""
import random

import polars as pl
import polars.testing as plt
import pytest
from pepmatch import Matcher
from pepmatch.benchmarker import Benchmarker
from pepmatch.matcher import best_match_budgets, k_for

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'


def test_k_is_the_largest_that_keeps_edits_plus_one_disjoint_seeds():
  assert [k_for(9, e) for e in range(5)] == [9, 4, 3, 2, 2]
  assert [k_for(15, e) for e in range(5)] == [15, 7, 5, 3, 3]
  assert k_for(40, 0) == k_for(40, 1) == 16     # capped at the longest k-mer an index stores
  assert k_for(3, 2) == 2                        # floored at the smallest index
  for length in range(2, 60):
    for edits in range(6):
      k = k_for(length, edits)
      assert 2 <= k <= 16
      if length >= 2 * (edits + 1):
        assert length // k >= edits + 1          # the pigeonhole guarantee


def test_best_match_budgets_cover_every_budget_under_the_guarantee():
  for length in range(2, 40):
    budgets = best_match_budgets(length)
    assert budgets == sorted(set(budgets)) and budgets[-1] == length
    guaranteed = [m for m in budgets if length >= 2 * (m + 1)]
    previous = -1
    for m in guaranteed:
      # one k serves every budget skipped since the last round, so none is searched worse
      assert {k_for(length, skipped) for skipped in range(previous + 1, m + 1)} == {k_for(length, m)}
      previous = m
    # past the guarantee nothing is skipped
    assert budgets[len(guaranteed):] == list(range(previous + 1, length + 1))


def corpus(seed):
  rng = random.Random(seed)
  proteins = [''.join(rng.choice(AMINO_ACIDS) for _ in range(rng.randint(40, 90))) for _ in range(12)]
  queries = []
  for _ in range(30):
    length = rng.randint(6, 24)
    source = rng.choice(proteins)
    start = rng.randrange(len(source) - length + 1)
    window = list(source[start:start + length])
    for position in rng.sample(range(length), rng.randint(0, 2)):
      window[position] = rng.choice(AMINO_ACIDS)
    queries.append(''.join(window))
  queries += [''.join(rng.choice(AMINO_ACIDS) for _ in range(rng.randint(6, 24))) for _ in range(4)]
  return proteins, queries


def write_proteome(tmp_path, proteins):
  path = tmp_path / 'proteome.fasta'
  path.write_text(''.join(f'>p{i + 1}\n{seq}\n' for i, seq in enumerate(proteins)))
  return path


def search(fasta, tmp_path, queries, **mode):
  return Matcher(query=queries, proteome_file=str(fasta), preprocessed_files_path=str(tmp_path), **mode).match()


MODES = [dict(), dict(max_mismatches=1), dict(max_mismatches=2), dict(max_indels=1), dict(max_indels=2)]


@pytest.mark.parametrize('mode', MODES)
def test_a_mixed_batch_returns_what_each_query_returns_alone(tmp_path, mode):
  # Each query's rows, and their order, depend only on that query; the batch is those
  # rows in the order the queries were given.
  proteins, queries = corpus(1)
  fasta = write_proteome(tmp_path, proteins)
  batch = search(fasta, tmp_path, queries, **mode)
  alone = []
  for i, query in enumerate(queries):
    frame = search(fasta, tmp_path, [query], **mode)
    alone.append(frame.with_columns(pl.lit(str(i + 1)).alias('Query ID')))
  expected = pl.concat(alone, how='vertical_relaxed')
  text = ['Matched Sequence', 'Indel Positions'] if 'max_indels' in mode else []
  plt.assert_frame_equal(batch, expected.with_columns(pl.col(text).cast(pl.Utf8)) if text else expected)
  assert batch['Matched Sequence'].is_not_null().sum() > 10


@pytest.mark.parametrize('mode', [dict(), dict(max_mismatches=2)])
def test_counts_of_a_mixed_batch_are_each_query_counted_alone(tmp_path, mode):
  proteins, queries = corpus(2)
  fasta = write_proteome(tmp_path, proteins)
  batch = search(fasta, tmp_path, queries, counts_only=True, **mode)
  alone = [search(fasta, tmp_path, [q], counts_only=True, **mode).with_columns(pl.lit(str(i + 1)).alias('Query ID'))
           for i, q in enumerate(queries)]
  plt.assert_frame_equal(batch, pl.concat(alone))
  assert batch.height >= 8


def test_each_length_is_searched_at_its_own_k(tmp_path, capsys):
  proteins, _ = corpus(3)
  fasta = write_proteome(tmp_path, proteins)
  queries = [proteins[0][:8], proteins[1][:15], proteins[2][:8], proteins[3][:24]]
  capsys.readouterr()
  df = search(fasta, tmp_path, queries, max_mismatches=1)
  out = capsys.readouterr().out
  assert 'Searching 2 peptides against proteome (k=4,' in out     # the two 8-mers
  assert 'Searching 1 peptides against proteome (k=7,' in out     # the 15-mer
  assert 'Searching 1 peptides against proteome (k=12,' in out    # the 24-mer
  assert sorted(p.name for p in tmp_path.glob('*.pepidx')) == [
    'proteome_12mers.pepidx', 'proteome_4mers.pepidx', 'proteome_7mers.pepidx']
  assert df.filter(pl.col('Mismatches') == 0)['Query Sequence'].unique(maintain_order=True).to_list() == queries


def test_a_short_query_does_not_change_how_a_long_one_is_searched(tmp_path, capsys):
  proteins, _ = corpus(4)
  fasta = write_proteome(tmp_path, proteins)
  long_query = proteins[0][5:25]
  search(fasta, tmp_path, [long_query], max_mismatches=1)
  capsys.readouterr()
  search(fasta, tmp_path, [long_query, 'ACDE'], max_mismatches=1)
  out = capsys.readouterr().out
  assert '(k=10,' in out and '(k=2,' in out


def test_only_queries_too_short_for_two_mer_seeds_lose_the_guarantee(tmp_path, capsys):
  proteins, _ = corpus(5)
  fasta = write_proteome(tmp_path, proteins)
  search(fasta, tmp_path, [proteins[0][:6], proteins[0][:20]], max_mismatches=2)
  assert 'not guaranteed' not in capsys.readouterr().out      # 6 = 2 * (2 + 1)
  search(fasta, tmp_path, [proteins[0][:5], proteins[0][:20]], max_mismatches=2)
  assert 'shorter than 6 residues (lengths: [5])' in capsys.readouterr().out


def test_best_match_of_mixed_lengths_finds_each_peptides_fewest_mismatches(tmp_path):
  proteins, queries = corpus(6)
  fasta = write_proteome(tmp_path, proteins)
  df = search(fasta, tmp_path, queries, best_match=True)
  assert df['Query Sequence'].to_list() == queries     # one row each, in the order given
  for query, found in zip(queries, df['Mismatches']):
    fewest = min(
      sum(a != b for a, b in zip(query, protein[i:i + len(query)]))
      for protein in proteins for i in range(len(protein) - len(query) + 1)
    )
    if len(query) >= 2 * (fewest + 1):     # within the recall guarantee the best is always found
      assert found == fewest, query
    else:
      assert found is None or found >= fewest, query


@pytest.mark.parametrize('max_mismatches', [0, 2, -1])   # -1 is the benchmark parameters' best match
def test_the_benchmark_hook_builds_every_index_its_search_reads(tmp_path, monkeypatch, capsys, max_mismatches):
  # Index builds are timed apart from the search, so none may be left for the search to do,
  # including for lengths the benchmark parameters do not list.
  proteins, queries = corpus(7)
  fasta = write_proteome(tmp_path, proteins)
  query_file = tmp_path / 'queries.txt'
  query_file.write_text('\n'.join(queries) + '\n')
  monkeypatch.chdir(tmp_path)
  bench = Benchmarker('test', query_file, fasta, lengths=[9], max_mismatches=max_mismatches, method_parameters={})
  bench.preprocess_proteome()
  capsys.readouterr()
  assert len(bench.search()) >= len(queries)
  assert 'No preprocessed file' not in capsys.readouterr().out
