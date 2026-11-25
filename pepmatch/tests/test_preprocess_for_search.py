"""Preprocessing takes no k. A search reads only the indexes in SEARCH_KS (k = 2 to 5), so
building those once prepares a proteome for every search: any lengths, any mode."""
import sys

import pytest
from pepmatch import Matcher, Preprocessor
from pepmatch.matcher import index_ks, k_for
from pepmatch.preprocessor import SEARCH_KS
from pepmatch.shell import run_matcher, run_preprocessor

FASTA = (
  '>sp|P11111|AAAA_HUMAN Alpha protein OS=Homo sapiens OX=9606 GN=AAA PE=1 SV=1\n'
  'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV\n'
  '>sp|P22222|BBBB_HUMAN Beta protein OS=Homo sapiens OX=9606 GN=BBB PE=1 SV=1\n'
  'MSDNGPQNQRNAPRITFGGPSDSTGSNQNGERSGARSKQRRPQGLPNNTASWFTALTQHGKEDLKFPRGQ\n'
)
QUERIES = ['MKT', 'MKTAYIAK', 'QISFVKSHF', 'SDNGPQNQRNAPRIT', 'LIEVQAPILSRVGDGTQDNLSGAEK',
           'RQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ', 'M1, K2, A4']


def setup(tmp_path):
  fasta = tmp_path / 'p.fasta'
  fasta.write_text(FASTA)
  queries = tmp_path / 'q.txt'
  queries.write_text('\n'.join(QUERIES) + '\n')
  return fasta, queries


def built(tmp_path):
  return sorted(int(p.name.split('_')[1].removesuffix('mers.pepidx')) for p in tmp_path.glob('*.pepidx'))


def test_no_search_reads_an_index_outside_the_search_ks():
  assert SEARCH_KS == (2, 3, 4, 5)
  for length in range(1, 80):
    for edits in range(8):
      assert k_for(length, edits) in SEARCH_KS
    assert set(index_ks([length], best_match=True)) <= set(SEARCH_KS)
  assert index_ks([9]) == [5] and index_ks([4]) == [4]
  assert index_ks([8, 9, 15, 25], 2) == [2, 3, 5]
  assert index_ks([40, 50]) == [5]
  assert index_ks([]) == []


@pytest.mark.parametrize('options', [
  {}, {'max_mismatches': 1}, {'max_mismatches': 2}, {'max_mismatches': 3},
  {'best_match': True}, {'max_mismatches': 2, 'counts_only': True},
])
def test_preprocess_leaves_nothing_for_any_search_to_build(tmp_path, capsys, options):
  fasta, queries = setup(tmp_path)
  paths = Preprocessor(str(fasta), preprocessed_files_path=str(tmp_path)).preprocess()
  assert built(tmp_path) == [2, 3, 4, 5] and len(paths) == 4
  capsys.readouterr()
  Matcher(str(queries), str(fasta), preprocessed_files_path=str(tmp_path), **options).match()
  assert 'No preprocessed file' not in capsys.readouterr().out
  assert built(tmp_path) == [2, 3, 4, 5]


@pytest.mark.parametrize('indels', [1, 2])
def test_preprocess_leaves_nothing_for_an_indel_search_to_build(tmp_path, capsys, indels):
  fasta, _ = setup(tmp_path)
  Preprocessor(str(fasta), preprocessed_files_path=str(tmp_path)).preprocess()
  capsys.readouterr()
  Matcher(QUERIES[:-1], str(fasta), max_indels=indels, preprocessed_files_path=str(tmp_path)).match()
  assert 'No preprocessed file' not in capsys.readouterr().out
  assert built(tmp_path) == [2, 3, 4, 5]


def test_preprocess_builds_what_a_search_builds_on_first_use(tmp_path):
  fasta, queries = setup(tmp_path)
  (tmp_path / 'a').mkdir()
  (tmp_path / 'b').mkdir()
  Preprocessor(str(fasta), preprocessed_files_path=str(tmp_path / 'a')).preprocess()
  Matcher(str(queries), str(fasta), 2, preprocessed_files_path=str(tmp_path / 'b')).match()
  assert set(built(tmp_path / 'b')) <= set(built(tmp_path / 'a'))
  for index in (tmp_path / 'b').glob('*.pepidx'):
    assert index.read_bytes() == (tmp_path / 'a' / index.name).read_bytes()


def test_preprocess_keeps_indexes_that_are_already_current(tmp_path, capsys):
  fasta, _ = setup(tmp_path)
  preprocessor = Preprocessor(str(fasta), preprocessed_files_path=str(tmp_path))
  preprocessor.preprocess()
  before = {p.name: p.stat().st_mtime_ns for p in tmp_path.glob('*.pepidx')}
  capsys.readouterr()
  preprocessor.preprocess()
  out = capsys.readouterr().out
  assert 'No preprocessed file' not in out and 'indexes ready' in out
  assert before == {p.name: p.stat().st_mtime_ns for p in tmp_path.glob('*.pepidx')}


def test_preprocess_rebuilds_after_the_fasta_changes(tmp_path, capsys):
  fasta, _ = setup(tmp_path)
  preprocessor = Preprocessor(str(fasta), preprocessed_files_path=str(tmp_path))
  preprocessor.preprocess()
  fasta.write_text(FASTA + '>sp|P33333|CCCC_HUMAN Gamma OS=Homo sapiens OX=9606 GN=CCC PE=1 SV=1\nMKTAYIAKQWWWWW\n')
  capsys.readouterr()
  preprocessor.preprocess()
  assert capsys.readouterr().out.count('changed since') == 4
  df = Matcher(['KQWWWWW'], str(fasta), preprocessed_files_path=str(tmp_path)).match()
  assert df['Protein ID'].to_list() == ['P33333.1']


def test_preprocess_with_a_k_still_builds_that_one_index(tmp_path):
  fasta, _ = setup(tmp_path)
  Preprocessor(str(fasta), preprocessed_files_path=str(tmp_path)).preprocess(k=9)
  assert built(tmp_path) == [9]


def test_cli_preprocess_takes_only_the_proteome(tmp_path, monkeypatch, capsys):
  fasta, queries = setup(tmp_path)
  monkeypatch.setattr(sys, 'argv', ['pepmatch-preprocess', '-p', str(fasta), '-P', str(tmp_path), '-n', 'named'])
  run_preprocessor()
  assert sorted(p.name for p in tmp_path.glob('*.pepidx')) == [f'named_{k}mers.pepidx' for k in (2, 3, 4, 5)]
  capsys.readouterr()
  monkeypatch.setattr(sys, 'argv', ['pepmatch-match', '-q', str(queries), '-p', str(fasta), '-m', '2',
                                    '-P', str(tmp_path), '-n', 'named', '-o', str(tmp_path / 'out')])
  run_matcher()
  assert 'No preprocessed file' not in capsys.readouterr().out
  assert len(list(tmp_path.glob('*.pepidx'))) == 4


@pytest.mark.parametrize('command, argv', [
  (run_preprocessor, ['-k', '5']),
  (run_matcher, ['-q', 'q.txt', '-k', '5']),
])
def test_cli_has_no_k_option(tmp_path, monkeypatch, command, argv):
  fasta, _ = setup(tmp_path)
  monkeypatch.chdir(tmp_path)
  monkeypatch.setattr(sys, 'argv', ['pepmatch', '-p', str(fasta), '-P', str(tmp_path), *argv])
  with pytest.raises(SystemExit):
    command()
  assert built(tmp_path) == []
