"""Building indexes ahead of time takes the search's own options, not a k: Matcher.preprocess()
and pepmatch-preprocess build exactly the indexes that search then reads."""
import sys

import pytest
from pepmatch import Matcher
from pepmatch.matcher import index_ks
from pepmatch.shell import run_matcher, run_preprocessor

FASTA = (
  '>sp|P11111|AAAA_HUMAN Alpha protein OS=Homo sapiens OX=9606 GN=AAA PE=1 SV=1\n'
  'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV\n'
  '>sp|P22222|BBBB_HUMAN Beta protein OS=Homo sapiens OX=9606 GN=BBB PE=1 SV=1\n'
  'MSDNGPQNQRNAPRITFGGPSDSTGSNQNGERSGARSKQRRPQGLPNNTASWFTALTQHGKEDLKFPRGQ\n'
)
QUERIES = ['MKTAYIAK', 'QISFVKSHF', 'SDNGPQNQRNAPRIT', 'LIEVQAPILSRVGDGTQDNLSGAEK']


def setup(tmp_path):
  fasta = tmp_path / 'p.fasta'
  fasta.write_text(FASTA)
  queries = tmp_path / 'q.txt'
  queries.write_text('\n'.join(QUERIES) + '\n')
  return fasta, queries


def built(tmp_path):
  return sorted(int(p.name.split('_')[1].removesuffix('mers.pepidx')) for p in tmp_path.glob('*.pepidx'))


def cli(monkeypatch, *argv):
  monkeypatch.setattr(sys, 'argv', list(argv))


def test_index_ks_is_one_k_per_distinct_rule_result():
  assert index_ks([9]) == [9]
  assert index_ks([8, 9, 15, 25], 2) == [2, 3, 5, 8]
  assert index_ks(range(8, 16), 2) == [2, 3, 4, 5]
  assert index_ks([40, 50]) == [16]                    # capped at the longest index
  assert index_ks([9], best_match=True) == [2, 3, 4, 9]
  assert index_ks([]) == []


@pytest.mark.parametrize('options', [
  {}, {'max_mismatches': 1}, {'max_mismatches': 3}, {'max_indels': 1}, {'max_indels': 2},
  {'best_match': True}, {'max_mismatches': 2, 'counts_only': True},
])
def test_preprocess_builds_exactly_what_the_search_reads(tmp_path, capsys, options):
  fasta, queries = setup(tmp_path)
  matcher = Matcher(str(queries), str(fasta), preprocessed_files_path=str(tmp_path), **options)
  paths = matcher.preprocess()
  assert built(tmp_path) == matcher.index_ks() and len(paths) == len(matcher.index_ks())
  capsys.readouterr()
  matcher.match()
  assert 'No preprocessed file' not in capsys.readouterr().out   # nothing left for the search to build
  assert built(tmp_path) == matcher.index_ks()


def test_preprocess_matches_a_search_that_builds_on_first_use(tmp_path):
  fasta, queries = setup(tmp_path)
  (tmp_path / 'a').mkdir()
  (tmp_path / 'b').mkdir()
  Matcher(str(queries), str(fasta), 2, preprocessed_files_path=str(tmp_path / 'a')).preprocess()
  Matcher(str(queries), str(fasta), 2, preprocessed_files_path=str(tmp_path / 'b')).match()
  assert built(tmp_path / 'a') == built(tmp_path / 'b')
  for index in (tmp_path / 'a').glob('*.pepidx'):
    assert index.read_bytes() == (tmp_path / 'b' / index.name).read_bytes()


def test_preprocess_keeps_indexes_that_are_already_current(tmp_path, capsys):
  fasta, queries = setup(tmp_path)
  matcher = Matcher(str(queries), str(fasta), 1, preprocessed_files_path=str(tmp_path))
  matcher.preprocess()
  before = {p.name: p.stat().st_mtime_ns for p in tmp_path.glob('*.pepidx')}
  capsys.readouterr()
  matcher.preprocess()
  out = capsys.readouterr().out
  assert 'No preprocessed file' not in out and 'indexes ready' in out
  assert before == {p.name: p.stat().st_mtime_ns for p in tmp_path.glob('*.pepidx')}


def test_discontinuous_epitopes_need_the_2mer_index(tmp_path):
  fasta, _ = setup(tmp_path)
  matcher = Matcher(['M1, K2, A4', 'QISFVKSHF'], str(fasta), preprocessed_files_path=str(tmp_path))
  assert matcher.index_ks() == [2, 9]
  matcher.preprocess()
  assert built(tmp_path) == [2, 9]


def test_cli_preprocess_takes_the_search_command_line(tmp_path, monkeypatch, capsys):
  fasta, queries = setup(tmp_path)
  search = ['-q', str(queries), '-p', str(fasta), '-m', '2', '-P', str(tmp_path)]
  cli(monkeypatch, 'pepmatch-preprocess', *search)
  run_preprocessor()
  assert built(tmp_path) == [2, 3, 5, 8]
  capsys.readouterr()
  cli(monkeypatch, 'pepmatch-match', *search, '-o', str(tmp_path / 'out'))
  run_matcher()
  assert 'No preprocessed file' not in capsys.readouterr().out
  assert built(tmp_path) == [2, 3, 5, 8]


@pytest.mark.parametrize('lengths, flags, expected', [
  ('9', [], [9]),
  ('8-15', ['-m', '2'], [2, 3, 4, 5]),
  ('8-11,15', ['-i', '1'], [4, 5, 7]),
  ('9', ['-b'], [2, 3, 4, 9]),
])
def test_cli_preprocess_by_lengths_when_there_is_no_query_yet(tmp_path, monkeypatch, lengths, flags, expected):
  fasta, _ = setup(tmp_path)
  cli(monkeypatch, 'pepmatch-preprocess', '-p', str(fasta), '-l', lengths, *flags, '-P', str(tmp_path), '-n', 'named')
  run_preprocessor()
  assert built(tmp_path) == expected
  assert all(p.name.startswith('named_') for p in tmp_path.glob('*.pepidx'))


def test_cli_preprocess_still_reads_the_old_proteome_flag(tmp_path, monkeypatch):
  fasta, _ = setup(tmp_path)
  cli(monkeypatch, 'pepmatch-preprocess', '--proteome', str(fasta), '-l', '9', '-P', str(tmp_path))
  run_preprocessor()
  assert built(tmp_path) == [9]


@pytest.mark.parametrize('argv', [
  ['-k', '5'],                                   # there is no k to give any more
  [],                                            # neither a query nor lengths
  ['-q', 'q.txt', '-l', '9'],                    # both
  ['-l', '0'], ['-l', '9-'], ['-l', 'nine'], ['-l', '15-8'],
])
def test_cli_preprocess_rejects_a_bad_command_line(tmp_path, monkeypatch, argv):
  fasta, _ = setup(tmp_path)
  cli(monkeypatch, 'pepmatch-preprocess', '-p', str(fasta), '-P', str(tmp_path), *argv)
  with pytest.raises(SystemExit):
    run_preprocessor()
  assert built(tmp_path) == []


def test_cli_preprocess_by_lengths_checks_the_search_options(tmp_path, monkeypatch):
  fasta, _ = setup(tmp_path)
  cli(monkeypatch, 'pepmatch-preprocess', '-p', str(fasta), '-l', '9', '-m', '1', '-i', '1', '-P', str(tmp_path))
  with pytest.raises(ValueError, match='mutually exclusive'):
    run_preprocessor()
