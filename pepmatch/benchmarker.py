from .preprocessor import Preprocessor
from .matcher import Matcher

class Benchmarker:
  def __init__(self, benchmark, query, proteome, lengths, max_mismatches, method_parameters):
    self.proteome = str(proteome)
    self.best_match = max_mismatches == -1   # the benchmark parameters' code for best match
    self.mm = 0 if self.best_match else max_mismatches
    self.matcher = Matcher(
      query=str(query),
      proteome_file=self.proteome,
      max_mismatches=self.mm,
      best_match=self.best_match,
      output_format='dataframe',
      sequence_version=False
    )

  def preprocess_proteome(self):
    """Builds every index the search reads, from the lengths actually in the query file
    (`lengths` in the benchmark parameters may leave some out). Always rebuilds, so the
    benchmark times the build."""
    for k in self.matcher.index_ks():
      Preprocessor(self.proteome).preprocess(k=k)

  def preprocess_query(self):
    raise TypeError('PEPMatch does not preprocess queries.')

  def search(self):
    return self.matcher.match().select(
      ['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start']
    ).to_pandas()
