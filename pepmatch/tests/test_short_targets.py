import pytest

from pathlib import Path

from pepmatch import Preprocessor
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
