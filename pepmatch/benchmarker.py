from .preprocessor import Preprocessor
from .matcher import Matcher, best_match_budgets, k_for

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

  def index_ks(self):
    """Every k the search reads an index for, from the lengths actually in the query file
    (`lengths` in the benchmark parameters may leave some out)."""
    lengths = {len(seq) for _, seq in self.matcher.query}
    if self.best_match:
      return sorted({k_for(length, m) for length in lengths for m in best_match_budgets(length)})
    return sorted({k_for(length, self.mm) for length in lengths})

  def preprocess_proteome(self):
    for k in self.index_ks():
      Preprocessor(self.proteome).preprocess(k=k)

  def preprocess_query(self):
    raise TypeError('PEPMatch does not preprocess queries.')

  def search(self):
    return self.matcher.match().select(
      ['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start']
    ).to_pandas()
