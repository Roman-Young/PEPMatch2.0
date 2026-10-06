from .preprocessor import Preprocessor
from .matcher import Matcher


def _query_lengths(query_path):
  lengths = []
  with open(query_path) as f:
    for line in f:
      line = line.strip()
      if line and not line.startswith('>'):
        lengths.append(len(line))
  return lengths


class Benchmarker:
  def __init__(
    self, benchmark: str, query: str, proteome: str, lengths: list, max_mismatches: int,
    method_parameters: dict, indels: int = 0
  ):
    self.query = str(query)
    self.proteome = str(proteome)
    self.lengths = lengths
    self.max_mismatches = max_mismatches
    self.indels = indels

    # Indel search derives its own pigeonhole k inside Matcher and is mutually
    # exclusive with mismatch search, so it must NOT run the mismatch k ladder below.
    if indels > 0:
      self.best_match = False
      self.mm = 0
      self.k = 0
      return

    if max_mismatches == -1:
      self.best_match = True
      self.mm = 0
      self.k = 0
    elif max_mismatches == 0:
      self.best_match = False
      self.mm = 0
      self.k = min(lengths)
    else:
      self.best_match = False
      self.mm = max_mismatches
      if benchmark == 'coronavirus':
        self.k = 2
      elif benchmark == 'neoepitopes':
        self.k = 3
      else:
        self.k = min(lengths) // (max_mismatches + 1)
        if self.k < 2:
          self.k = 2

  def __str__(self):
    return 'PEPMatch'

  def _indel_k_values(self):
    """The k(s) indel_search will actually use, so preprocessing builds those indexes
    and index-build time is charged to proteome preprocessing rather than to search.
    Mirrors Matcher.indel_search's query partition (matcher.py): queries too short for
    k>=3 run on a k=2 table, the rest at max(2, min_len // (indels + 1))."""
    n = self.indels
    lengths = _query_lengths(self.query)
    ks = set()
    short = [L for L in lengths if L // (n + 1) < 3]
    rest = [L for L in lengths if L // (n + 1) >= 3]
    if short:
      ks.add(2)
    if rest:
      ks.add(max(2, min(rest) // (n + 1)))
    return ks

  def preprocess_proteome(self):
    if self.indels > 0:
      for k in self._indel_k_values():
        Preprocessor(self.proteome).preprocess(k=k)
    elif self.best_match:
      for k in [15, 7, 3, 2]:
        Preprocessor(self.proteome).preprocess(k=k)
    else:
      Preprocessor(self.proteome).preprocess(k=self.k)

  def preprocess_query(self):
    raise TypeError('PEPMatch does not preprocess queries.')

  def search(self):
    if self.indels > 0:
      df = Matcher(
        query=self.query,
        proteome_file=self.proteome,
        max_indels=self.indels,
        output_format='dataframe',
        sequence_version=False
      ).match()
    else:
      df = Matcher(
        query=self.query,
        proteome_file=self.proteome,
        max_mismatches=self.mm,
        k=self.k,
        best_match=self.best_match,
        output_format='dataframe',
        sequence_version=False
      ).match()

    return df.select(
      ['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start']
    ).to_pandas()
