from .preprocessor import Preprocessor
from .matcher import Matcher, best_match_budgets, k_for

class Benchmarker:
  def __init__(self, benchmark, query, proteome, lengths, max_mismatches, method_parameters):
    self.query = str(query)
    self.proteome = str(proteome)
    self.lengths = lengths
    self.best_match = max_mismatches == -1   # the benchmark parameters' code for best match
    self.mm = 0 if self.best_match else max_mismatches

  def index_ks(self):
    """Every k the search will read an index for: each query length at its own k."""
    if self.best_match:
      return sorted({k_for(length, m) for length in self.lengths for m in best_match_budgets(length)})
    return sorted({k_for(length, self.mm) for length in self.lengths})

  def preprocess_proteome(self):
    for k in self.index_ks():
      Preprocessor(self.proteome).preprocess(k=k)

  def preprocess_query(self):
    raise TypeError('PEPMatch does not preprocess queries.')

  def search(self):
    df = Matcher(
      query=self.query,
      proteome_file=self.proteome,
      max_mismatches=self.mm,
      best_match=self.best_match,
      output_format='dataframe',
      sequence_version=False
    ).match()

    return df.select(
      ['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start']
    ).to_pandas()
