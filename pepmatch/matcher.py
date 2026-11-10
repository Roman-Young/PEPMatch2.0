import os
import tempfile
import polars as pl
from pathlib import Path
from Bio import SeqIO
from ._rs import rs_preprocess, rs_match, rs_discontinuous, rs_metadata, rs_match_counts, rs_indel_match, rs_short_lengths, rs_short_target_match

VALID_OUTPUT_FORMATS = ['dataframe', 'polars', 'pandas', 'csv', 'tsv', 'xlsx', 'json']
DATAFRAME_OUTPUT_FORMATS = frozenset({'dataframe', 'polars', 'pandas'})
FASTA_EXTENSIONS = {
  '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa', '.frn'
}

# Shortest proteome entry the short-target search looks for inside the queries: a 1-residue
# entry cannot hold even one 2-mer seed. A placeholder until the floor is settled with the PI.
MIN_TARGET_LEN = 2
# preprocess.rs packs each k-mer into a u128, one byte per residue.
MAX_INDEX_K = 16

def output_matches(df: pl.DataFrame, output_format: str, output_name: str) -> None:
  """Write the results to disk in `output_format`, appending the extension if absent."""
  path = output_name.__str__()
  if not path.lower().endswith(f".{output_format}"):
    path += f".{output_format}"
  if output_format == 'csv':
    df.write_csv(path)
  elif output_format == 'tsv':
    df.write_csv(path, separator='\t')
  elif output_format == 'xlsx':
    df.write_excel(path)
  elif output_format == 'json':
    df.write_json(path)


def _return_dataframe(df: pl.DataFrame, output_format: str):
  """Return the finished frame as Polars ('dataframe' is an alias for it) or as pandas.

  The pandas conversion needs pyarrow, so pandas-only installs fall back to building
  the frame column by column."""
  if output_format in ('dataframe', 'polars'):
    return df
  try:
    import pandas as pd
  except ImportError as e:
    raise ImportError(
      "output_format='pandas' requires pandas. "
      "Install with: pip install 'pepmatch[pandas]'"
    ) from e
  try:
    return df.to_pandas()
  except ModuleNotFoundError:
    return pd.DataFrame({col: df[col].to_list() for col in df.columns})


def _indel_placements(query, matched):
  """Every valid placement of the edits between `query` and `matched`, each a tuple of
  (1-based query position, residue) per edit; a repeat makes placements equivalent.
  Query-terminal deletions and boundary insertions are barred, and an inserted
  residue's position counts the query residues preceding it.

  Removing n positions from the longer string splits it into n + 1 stretches, and
  stretch t must equal the shorter string shifted by t residues. Placements are read
  off precomputed agreement runs for each shift, left to right so they come out in
  lexicographic order (_indel_edits reads placements[0]), instead of rebuilding the
  string for every combination of positions. With one edit that reduces to two numbers:
  position c is valid when the strings agree before it (common prefix) and after it
  (common suffix), so the valid positions are one run, read off those two lengths."""
  L, M = len(query), len(matched)
  n = abs(M - L)
  if n == 0:
    return []
  longer, shorter = (query, matched) if M < L else (matched, query)
  N = len(longer)
  if n == 1:
    prefix = 0
    while prefix < N - 1 and longer[prefix] == shorter[prefix]:
      prefix += 1
    suffix = 0
    while suffix < N - 1 and longer[N - 1 - suffix] == shorter[N - 2 - suffix]:
      suffix += 1
    combos = [(c,) for c in range(max(1, N - 1 - suffix), min(prefix, N - 2) + 1)]
  else:
    # runs[s][x]: consecutive positions from x where longer[y + s] == shorter[y].
    runs = []
    for s in range(n + 1):
      run = [0] * (N + 1)
      for x in range(len(shorter) - 1, -1, -1):
        if x + s < N and longer[x + s] == shorter[x]:
          run[x] = run[x + 1] + 1
      runs.append(run)

    combos = []

    def place(t, start, combo):
      if t == n:
        if runs[n][start - n] >= N - start:
          combos.append(combo)
        return
      for c in range(max(start, 1), N - 1):   # never an end residue: the terminal rule
        if c - start > runs[t][start - t]:
          break
        place(t + 1, c + 1, combo + (c,))

    place(0, 0, ())
  if M < L:
    return [tuple((i + 1, query[i]) for i in combo) for combo in combos]
  return [tuple((k - u + 1, matched[k]) for u, k in enumerate(combo)) for combo in combos]


def _indel_edits(query, matched):
  """Entries for the Indel Positions column as (kind, residues, low, high).

  One edit gives one position range. Two edits at the same site collapse to a chunk,
  ranged only across starts that remove identical residues (in a periodic repeat the
  missing pair changes along the range); independent edits each get their own range."""
  placements = _indel_placements(query, matched)
  if not placements:
    return []
  kind = 'd' if len(matched) < len(query) else 'i'

  if len(placements[0]) == 1:
    positions = sorted(p for ((p, _),) in placements)
    residue = placements[0][0][1]
    return [(kind, residue, positions[0], positions[-1])]

  chunks = sorted({
    (p1, r1 + r2) for (p1, r1), (p2, r2) in placements
    if (p2 == p1 + 1 if kind == 'd' else p2 == p1)
  })
  if chunks:
    entries = []
    for start, residues in chunks:
      if entries and entries[-1][1] == residues and entries[-1][3] + 1 == start:
        k, r, low, _ = entries[-1]
        entries[-1] = (k, r, low, start)
      else:
        entries.append((kind, residues, start, start))
    return entries

  entries = []
  for slot in (0, 1):
    positions = sorted({pl[slot][0] for pl in placements})
    entries.append((kind, placements[0][slot][1], positions[0], positions[-1]))
  return entries


def format_indel_positions(query, matched, offset=0):
  """Render the Indel Positions column, e.g. 'd: A[6]', 'd: AA[2,4]', 'd: C[3], d: E[5]',
  or '[]' for an exact match. A 1-based range collapses when the placement is unambiguous.
  `offset` shifts every position, for when `query` is a window starting `offset` residues
  into the full query."""
  edits = _indel_edits(query, matched)
  if not edits:
    return '[]'
  return ', '.join(
    f'{kind}: {residues}[{low + offset}]' if low == high
    else f'{kind}: {residues}[{low + offset},{high + offset}]'
    for kind, residues, low, high in edits
  )


class Matcher:
  """Searches query peptides against a preprocessed proteome index
  and returns matches as a DataFrame (Polars by default, or pandas) or output file."""

  def __init__(
    self,
    query,
    proteome_file,
    max_mismatches=0,
    max_indels=0,
    k=0,
    preprocessed_files_path='.',
    best_match=False,
    output_format='dataframe',
    output_name='',
    sequence_version=True,
    counts_only=False,
    short_targets=False,
  ):

    if k == 0:
      self.k = 0
      self.k_specified = False
    else:
      self.k = k
      self.k_specified = True

    if k != 0 and k < 2:
      raise ValueError('k must be >= 2.')

    if max_indels > 2:
      raise ValueError('max_indels > 2 is not yet supported. Only max_indels<=2 has been validated.')

    if max_indels > 0 and max_mismatches > 0:
      raise ValueError('max_indels and max_mismatches are mutually exclusive.')

    if max_indels > 0 and best_match:
      raise ValueError('max_indels and best_match are not yet supported together.')

    if max_indels > 0 and counts_only:
      raise ValueError('max_indels and counts_only are not yet supported together.')

    self.max_indels = max_indels

    if output_format not in VALID_OUTPUT_FORMATS:
      raise ValueError(
        f'Invalid output format. Choose from: {VALID_OUTPUT_FORMATS}'
      )

    if counts_only and best_match:
      raise ValueError('counts_only is not supported with best_match.')

    if short_targets and best_match:
      raise ValueError('short_targets and best_match are not yet supported together.')

    if short_targets and counts_only:
      raise ValueError('short_targets and counts_only are not yet supported together.')

    self.proteome_file = str(proteome_file)
    self.proteome_name = Path(proteome_file).stem
    self.max_mismatches = max_mismatches
    self.preprocessed_files_path = preprocessed_files_path
    self.best_match = best_match
    self.output_format = output_format
    self.sequence_version = sequence_version
    self.counts_only = counts_only
    self.short_targets = short_targets
    self.output_name = output_name or 'PEPMatch_results'

    self.query = self._parse_query(query)
    self.discontinuous_epitopes = self._find_discontinuous_epitopes()
    if self.max_indels > 0 and self.discontinuous_epitopes:
      raise ValueError('max_indels and discontinuous epitopes are not supported together.')
    self.query = self._clean_query()
    if not self.query and not self.discontinuous_epitopes:
      raise ValueError('Query is empty.')

    # The short-target search indexes the queries through a FASTA file, whose reader trims
    # each line, so a space or line break in a list query would shift every Query start/end.
    if short_targets:
      for qid, seq in self.query:
        if any(c.isspace() for c in seq):
          raise ValueError(f'Query {qid} contains whitespace, which short_targets cannot place; remove it.')

    if self.k_specified and self.query:
      shortest = min(len(seq) for _, seq in self.query)
      if self.k > shortest:
        raise ValueError(
          f'k ({self.k}) cannot exceed the shortest query peptide length ({shortest}); '
          f'shorter peptides would be silently skipped.'
        )

  def _parse_query(self, query):
    """Normalize the query into (id, sequence) pairs from a list, .txt file, or FASTA."""
    if isinstance(query, list):
      return [(str(i + 1), seq.upper()) for i, seq in enumerate(query)]

    path = Path(query)
    if not path.is_file():
      raise FileNotFoundError(f'Query file not found: {query}')

    if path.suffix.lower() in FASTA_EXTENSIONS:
      return [
        (record.id, str(record.seq).upper())
        for record in SeqIO.parse(str(path), 'fasta')
      ]
    elif path.suffix.lower() == '.txt':
      with open(path) as f:
        lines = [line.strip() for line in f if line.strip()]
      return [(str(i + 1), seq.upper()) for i, seq in enumerate(lines)]

    raise ValueError(
      f'Unsupported query format: {path.suffix}. '
      f'Use a Python list, .txt file, or FASTA file for the query.'
    )

  def _find_discontinuous_epitopes(self):
    """Pick out the queries written as residue/position lists, e.g. 'N1, A3, L5'."""
    discontinuous_epitopes = {}
    for query_id, peptide in self.query:
      try:
        epitope = [(x[0], int(x[1:])) for x in peptide.split(', ')]
        discontinuous_epitopes[query_id] = epitope
      except (ValueError, IndexError):
        continue
    return discontinuous_epitopes

  def _clean_query(self):
    """Drop the discontinuous epitopes from the linear query set."""
    discontinuous_ids = set(self.discontinuous_epitopes.keys())
    return [
      (qid, seq) for qid, seq in self.query
      if qid not in discontinuous_ids
    ]

  def match(self):
    """Run the search, returning a DataFrame or writing the results to `output_name`."""
    linear_df = pl.DataFrame()
    discontinuous_df = pl.DataFrame()

    if self.counts_only:
      k = self.k if self.k_specified else self._auto_k(self.max_mismatches)
      self._recall_warning(k, self.max_mismatches)
      df = self._counts_to_dataframe(self._search_counts(k, self.max_mismatches))
      if self.output_format in DATAFRAME_OUTPUT_FORMATS:
        return _return_dataframe(df, self.output_format)
      return output_matches(df, self.output_format, self.output_name)

    if self.query:
      if self.max_indels > 0:
        linear_df = self.indel_search()
        pepidx_path = self._indel_pepidx_path
      elif self.best_match and self.k_specified:
        results, pepidx_path = self._search(self.k, self.max_mismatches)
        linear_df = self._best_match_filter(self._to_dataframe(results, pepidx_path))
      elif self.best_match:
        linear_df = self.best_match_search()
      else:
        k = self.k if self.k_specified else self._auto_k(self.max_mismatches)
        self._recall_warning(k, self.max_mismatches)
        results, pepidx_path = self._search(k, self.max_mismatches)
        linear_df = self._to_dataframe(results, pepidx_path)

      if self.short_targets:   # best_match is rejected in __init__, so pepidx_path is set
        linear_df = self._merge_short_targets(linear_df, self._short_target_search(pepidx_path))

    if self.discontinuous_epitopes:
      pepidx_path = self._pepidx_path(2)
      if not os.path.isfile(pepidx_path):
        rs_preprocess(self.proteome_file, 2, pepidx_path)
      epitopes = [
        (qid, residues) for qid, residues in self.discontinuous_epitopes.items()
      ]
      results = rs_discontinuous(pepidx_path, epitopes, self.max_mismatches)
      discontinuous_df = self._to_dataframe(results, pepidx_path)
      if self.short_targets:   # a position list has no query span; null keeps the concat aligned
        discontinuous_df = discontinuous_df.with_columns(
          pl.lit(None, dtype=pl.Int64).alias('Query start'),
          pl.lit(None, dtype=pl.Int64).alias('Query end'),
        ).select(self._span_columns(self._final_columns(False)))

    dfs = [d for d in [linear_df, discontinuous_df] if d.height > 0]
    df = pl.concat(dfs, how="vertical_relaxed") if dfs else linear_df

    if self.output_format in DATAFRAME_OUTPUT_FORMATS:
      return _return_dataframe(df, self.output_format)
    output_matches(df, self.output_format, self.output_name)

  def indel_search(self):
    """Indel search, partitioned by the k each group of queries needs.

    max_indels=2 is homogeneous: two insertions or two deletions, never one of each.
    Queries too short for k=3 run against a 2-mer table on their own rather than
    dragging the whole batch down to k=2, where 2-mer seeds explode. Every group still
    runs at k <= its own pigeonhole optimum, so the match set equals a single-k run.
    Each group's hits are labeled from the index that found them: protein numbers are
    that index's own FASTA ordinals, so labeling them from another group's index would
    name the wrong protein whenever the two were built from different versions of the
    proteome. vertical_relaxed lets an all-miss group's null column stack under another
    group's String column, and the text columns are cast to String at the end so a batch
    that hits nothing has the same schema as one that does."""
    n = self.max_indels
    short = [(q, s) for q, s in self.query if len(s) // (n + 1) < 3]
    rest = [(q, s) for q, s in self.query if len(s) // (n + 1) >= 3]

    groups = {}
    if short:
      groups.setdefault(self._clamp_k(2), []).extend(short)
    if rest:
      rest_optimal = max(2, min(len(s) for _, s in rest) // (n + 1))
      groups.setdefault(self._clamp_k(rest_optimal), []).extend(rest)

    frames = []
    for k, queries in sorted(groups.items()):
      self._recall_warning(k, n, queries)
      pepidx_path = self._pepidx_path(k)
      if not os.path.isfile(pepidx_path):
        print(f"No preprocessed file found for k={k}, building index now "
              f"(this may take a moment)...")
        rs_preprocess(self.proteome_file, k, pepidx_path)
      print(f"Searching {len(queries)} peptides against {self.proteome_name} "
            f"(k={k}, max_indels={n})...")
      results = rs_indel_match(pepidx_path, queries, n)
      frames.append(self._to_dataframe(results, pepidx_path, is_indels=True))
    # Any group's index serves the short-target search: target sequences and protein numbers
    # are the same at every k.
    self._indel_pepidx_path = pepidx_path
    df = pl.concat(frames, how='vertical_relaxed')
    return df.with_columns(pl.col('Matched Sequence', 'Indel Positions').cast(pl.Utf8))

  def _clamp_k(self, optimal):
    """Resolve k for one query group: an explicit k is honored up to `optimal`, above
    which too few disjoint seeds remain to guarantee recall, so it is clamped down."""
    if self.k_specified and self.k > optimal:
      print(f"Requested k={self.k} exceeds k={optimal}, the largest k that "
            f"guarantees complete recall for these queries (pigeonhole); "
            f"using k={optimal}.")
      return optimal
    return self.k if self.k_specified else optimal

  def _recall_warning(self, k, edits, queries=None, what='query peptides', lengths=None):
    """Warn where k cannot guarantee complete recall. The seed guarantee needs
    len >= k * (edits + 1), so a shorter peptide may have a real match that every seed
    misses. Shared by the mismatch and indel paths, the latter passing its own group, and by
    the short-target search, whose searched peptides are proteome entries it passes as
    `lengths` only."""
    if lengths is None:
      lengths = [len(s) for _, s in (self.query if queries is None else queries)]
    unguaranteed = sorted({length for length in lengths if length < k * (edits + 1)})
    if unguaranteed:
      print(f"k={k}, edits={edits}: complete recall is not guaranteed for {what} "
            f"shorter than {k * (edits + 1)} residues (lengths: {unguaranteed}); "
            f"a match may exist but be missed.")

  def _auto_k(self, edits):
    """Largest k that still guarantees recall (pigeonhole): floor(L / (edits + 1)) taken
    on the shortest query peptide, so every peptide splits into edits+1 disjoint seeds.
    Floored at 2, the preprocessor minimum."""
    min_len = min(len(seq) for _, seq in self.query)
    return max(2, min_len // (edits + 1))

  def _pepidx_path(self, k):
    return os.path.join(
      self.preprocessed_files_path, f'{self.proteome_name}_{k}mers.pepidx'
    )

  def _search(self, k, max_mismatches, peptides=None):
    """Search `peptides` (default: the whole query) at k, building the index if absent.
    Returns the hit columns and the .pepidx they came from, so the metadata join reads
    the index this search used rather than rediscovering one by name."""
    pepidx_path = self._pepidx_path(k)
    if not os.path.isfile(pepidx_path):
      print(f"Preprocessing {self.proteome_name} with k={k}...")
      rs_preprocess(self.proteome_file, k, pepidx_path)
    query = peptides or self.query
    print(f"Searching {len(query)} peptides against {self.proteome_name} (k={k}, max_mismatches={max_mismatches})...")
    return rs_match(pepidx_path, query, k, max_mismatches), pepidx_path

  def _short_target_search(self, pepidx_path):
    """Hits where a proteome entry shorter than the query lies inside it, within the same
    edit budget. The forward search lays the whole query on a protein, so an entry shorter
    than the query never fits and is missed. Here the roles swap: the queries become a small
    index and each short entry is searched against it with the same engine. The hits are
    then reported from the query's side: Matched Sequence is the whole entry, Index start/end
    span all of it, Query start/end say where it landed in the query, and edit positions are
    numbered along the query. Returns None when nothing was found."""
    is_indels = self.max_indels > 0
    edits = self.max_indels if is_indels else self.max_mismatches
    longest = max(len(seq) for _, seq in self.query)
    # (length, count) of the entries to search. The entries themselves stay in the index:
    # rs_short_target_match reads them there and hands back only hits, so nothing that scales
    # with the number of entries crosses into Python.
    lengths = rs_short_lengths(pepidx_path, MIN_TARGET_LEN, longest)
    if not lengths:
      return None

    def lengths_in(lo, hi):
      return [length for length, _ in lengths if lo <= length < hi]

    # k is chosen the way the forward search chooses it, with the entries in the queries'
    # place: one k from the shortest entry, as _auto_k does, or for indels indel_search's
    # split, where entries too short for k=3 (length // (edits + 1) < 3) get their own 2-mer
    # index. Unlike there it is capped at MAX_INDEX_K, since entries run up to the longest
    # query's length. Each group is a length range [lo, hi); an empty one is skipped.
    def k_for(lo, hi):
      return min(MAX_INDEX_K, max(2, min(lengths_in(lo, hi)) // (edits + 1)))
    if is_indels:
      split = 3 * (edits + 1)
      short, rest = (MIN_TARGET_LEN, min(split, longest)), (max(MIN_TARGET_LEN, split), longest)
      groups = []
      if lengths_in(*short):
        groups.append((2, *short))
      if lengths_in(*rest):
        groups.append((k_for(*rest), *rest))
    else:
      groups = [(k_for(MIN_TARGET_LEN, longest), MIN_TARGET_LEN, longest)]

    hits = []
    with tempfile.TemporaryDirectory() as tmp:
      queries_fasta = os.path.join(tmp, 'queries.fasta')
      with open(queries_fasta, 'w') as f:
        for i, (_, seq) in enumerate(self.query):   # protein number i + 1 is query i
          f.write(f'>{i + 1}\n{seq}\n')

      print(f"Searching {sum(count for _, count in lengths)} entries of {self.proteome_name} shorter than the "
            f"longest query inside the queries ({'max_indels' if is_indels else 'max_mismatches'}={edits})...")
      for k, lo, hi in groups:
        self._recall_warning(k, edits, what='proteome entries', lengths=lengths_in(lo, hi))
        queries_index = os.path.join(tmp, f'queries_{k}mers.pepidx')
        rs_preprocess(queries_fasta, k, queries_index)
        cols = rs_short_target_match(pepidx_path, queries_index, lo, hi, k, edits, is_indels)
        for t_num, t_seq, win, q_num, n, qs, qe in zip(*cols):
          query_id, query_seq = self.query[q_num - 1]
          # An entry as long as the query is the forward search's to find.
          if len(t_seq) < len(query_seq):
            hits.append((q_num, t_num, qs, qe, query_id, query_seq, t_seq, win, n))

    if not hits:
      return None
    hits.sort(key=lambda h: h[:4])   # query order, then protein, then where in the query

    labels = {}
    positions = []
    for _, _, qs, _, _, _, t_seq, win, _ in hits:
      if is_indels:
        # Labelled on the query window the entry aligned to, so an ambiguous repeat range
        # stays inside that window, then shifted to count from the query's first residue.
        key = (win, t_seq, qs)
        if key not in labels:
          labels[key] = format_indel_positions(win, t_seq, qs - 1)
        positions.append(labels[key])
      else:
        positions.append('[' + ', '.join(str(qs + i) for i, (a, b) in enumerate(zip(win, t_seq)) if a != b) + ']')

    # Through _to_dataframe for the metadata join and Protein ID versioning, with the
    # positions already computed: is_indels=False takes them as given, and the two columns
    # are renamed back for an indel search.
    cols = (
      [h[4] for h in hits], [h[5] for h in hits], [h[6] for h in hits], [h[1] for h in hits],
      [h[8] for h in hits], positions, [1] * len(hits), [len(h[6]) for h in hits],
    )
    df = self._to_dataframe(cols, pepidx_path)
    if is_indels:
      df = df.rename({'Mismatches': 'Indels', 'Mutated Positions': 'Indel Positions'})
    return df.with_columns(
      pl.Series('Query start', [h[2] for h in hits], dtype=pl.Int64),
      pl.Series('Query end', [h[3] for h in hits], dtype=pl.Int64),
    ).select(self._span_columns(self._final_columns(is_indels)))

  def _merge_short_targets(self, forward, short):
    """Add the short-target hits to the forward results without changing any forward row.
    Forward hits span the whole query. A short-target hit the forward search already
    reported -- an entry shorter than the query only by the indels it absorbs -- is kept
    in its forward row. A query's miss row goes once it has a short-target hit."""
    is_indels = self.max_indels > 0
    columns = self._span_columns(self._final_columns(is_indels))
    hit = pl.col('Matched Sequence').is_not_null()
    forward = forward.with_columns(
      pl.when(hit).then(pl.lit(1, dtype=pl.Int64)).alias('Query start'),
      pl.when(hit).then(pl.col('Query Sequence').str.len_chars().cast(pl.Int64)).alias('Query end'),
    ).select(columns)
    if short is None:
      return forward

    # The query span is part of the key: an entry can sit at several places in one query,
    # and only the placement covering the whole query is the forward row's alignment.
    key = ['Query ID', 'Query Sequence', 'Protein ID', 'Matched Sequence',
           'Index start', 'Index end', 'Query start', 'Query end']
    # A forward search that hit nothing leaves its text columns typed Null, which cannot be
    # joined against strings.
    forward_keys = forward.filter(hit).select(key).cast(dict(short.select(key).schema))
    short = short.join(forward_keys, on=key, how='anti', maintain_order='left')

    gained = short.select('Query ID', 'Query Sequence').unique()
    forward = forward.with_row_index('_row')
    resolved = forward.filter(~hit).join(gained, on=['Query ID', 'Query Sequence'], how='semi').get_column('_row')
    forward = forward.filter(~pl.col('_row').is_in(resolved.implode())).drop('_row')
    return pl.concat([forward, short], how='vertical_relaxed')

  def _span_columns(self, columns):
    """`columns` with Query start / Query end placed after Index end."""
    i = columns.index('Index end') + 1
    return columns[:i] + ['Query start', 'Query end'] + columns[i:]

  def _search_counts(self, k, max_mismatches):
    """Count hits per (peptide, mismatch level) at k, building the index if absent."""
    pepidx_path = self._pepidx_path(k)
    if not os.path.isfile(pepidx_path):
      print(f"Preprocessing {self.proteome_name} with k={k}...")
      rs_preprocess(self.proteome_file, k, pepidx_path)
    print(f"Counting {len(self.query)} peptides against {self.proteome_name} (k={k}, max_mismatches={max_mismatches})...")
    return rs_match_counts(pepidx_path, self.query, k, max_mismatches)

  def _counts_to_dataframe(self, cols):
    """Aggregate counts: one row per (query peptide, mismatch level) with a hit.
    Memory is O(unique queries), independent of total hit count."""
    qid, qseq, mm, count = cols
    if not qid:
      return pl.DataFrame(schema={
        'Query ID': pl.Utf8, 'Query Sequence': pl.Utf8,
        'Mismatches': pl.Int64, 'Count': pl.UInt64,
      })
    return pl.DataFrame({
      'Query ID': qid,
      'Query Sequence': qseq,
      'Mismatches': pl.Series(mm, dtype=pl.Int64),
      'Count': pl.Series(count, dtype=pl.UInt64),
    })

  def best_match_search(self):
    """Closest match per peptide, walking k down and the mismatch budget up until
    nothing is left unmatched."""
    peptides_remaining = self.query.copy()
    acc = tuple([] for _ in range(8))   # 8 columnar accumulators
    pepidx_path = None

    def collect_matched(cols):
      matched_ids = set()
      matched_col = cols[2]
      for i in range(len(cols[0])):
        if matched_col[i] is not None:
          for j in range(8):
            acc[j].append(cols[j][i])
          matched_ids.add(cols[0][i])
      return matched_ids

    initial_k = min(len(seq) for _, seq in peptides_remaining)
    ks = [initial_k]
    while initial_k > 2:
      initial_k //= 2
      ks.append(max(initial_k, 2))
    ks = sorted(set(ks), reverse=True)

    for k in ks:
      if not peptides_remaining:
        break
      min_len = min(len(seq) for _, seq in peptides_remaining)
      max_mm = (min_len // k) - 1
      if max_mm < 0:
        max_mm = 0
      results, pepidx_path = self._search(k, max_mm, peptides_remaining)
      matched_ids = collect_matched(results)
      peptides_remaining = [
        (qid, seq) for qid, seq in peptides_remaining if qid not in matched_ids
      ]
      print(f"  -> k={k}, mismatches<={max_mm}: {len(peptides_remaining)} remaining")

    if peptides_remaining:
      max_mm = min(len(seq) for _, seq in peptides_remaining) // 2
      while peptides_remaining:
        shortest_len = min(len(seq) for _, seq in peptides_remaining)
        if max_mm >= shortest_len:
          break
        max_mm += 1
        results, pepidx_path = self._search(2, max_mm, peptides_remaining)
        matched_ids = collect_matched(results)
        peptides_remaining = [
          (qid, seq) for qid, seq in peptides_remaining if qid not in matched_ids
        ]
        print(f"  -> k=2, mismatches<={max_mm}: {len(peptides_remaining)} remaining")

    for qid, seq in peptides_remaining:   # leftover unmatched
      for j, val in enumerate((qid, seq, None, None, None, "[]", None, None)):
        acc[j].append(val)

    df = self._to_dataframe(acc, pepidx_path)
    return self._best_match_filter(df)

  def _best_match_filter(self, df):
    """Keep one row per query, breaking ties in order: fewest mismatches, gene priority,
    canonical (unhyphenated) Protein ID, SwissProt reviewed, protein existence level,
    non-fragment. Unmatched rows pass through untouched."""
    matched_df = df.filter(pl.col("Matched Sequence").is_not_null())
    unmatched_df = df.filter(pl.col("Matched Sequence").is_null())

    if matched_df.height > 0:
      matched_df = matched_df.with_columns(
        pl.col("Mismatches").min().over("Query ID").alias("_min_mm")
      ).filter(
        pl.col("Mismatches") == pl.col("_min_mm")
      ).drop("_min_mm")

      matched_df = matched_df.with_columns(
        pl.col("Gene Priority").max().over("Query ID").alias("_max_gp")
      ).filter(
        pl.col("Gene Priority") == pl.col("_max_gp")
      ).drop("_max_gp")

      matched_df = (
        matched_df.with_columns(
          (~pl.col("Protein ID").str.contains("-")).cast(pl.Int8).alias("_is_canonical")
        ).with_columns(
          pl.col("_is_canonical").max().over("Query ID").alias("_max_canonical")
        ).filter(
          pl.col("_is_canonical") == pl.col("_max_canonical")
        ).drop(["_is_canonical", "_max_canonical"])
      )

      matched_df = matched_df.with_columns(
        pl.col("SwissProt Reviewed").cast(pl.Int8).max().over("Query ID").alias("_max_reviewed")
      ).filter(
        pl.col("SwissProt Reviewed").cast(pl.Int8) == pl.col("_max_reviewed")
      ).drop("_max_reviewed")

      matched_df = matched_df.with_columns(
        pl.col("Protein Existence Level").min().over("Query ID").alias("_min_pe")
      ).filter(
        pl.col("Protein Existence Level") == pl.col("_min_pe")
      ).drop("_min_pe")

      matched_df = (
        matched_df.with_columns(
          (~pl.col("Protein Name").str.contains("Fragment")).alias("_not_fragment")
        ).with_columns(
          pl.col("_not_fragment").any().over("Query ID").alias("_has_non_fragment")
        ).filter(
          (pl.col("_not_fragment") & pl.col("_has_non_fragment")) |
          (~pl.col("_has_non_fragment"))
        ).drop(["_not_fragment", "_has_non_fragment"])
      )

      matched_df = matched_df.unique(subset=["Query ID"], keep="first")

    return pl.concat([matched_df, unmatched_df], how="vertical_relaxed")

  def _metadata_table(self, pepidx_path: str, protein_nums=None) -> pl.DataFrame:
    """Metadata for `protein_nums` (every protein when None), read from `pepidx_path`, the
    index the search just used. Protein numbers are that index's own FASTA ordinals, so any
    other file -- including one whose name merely shares this proteome's prefix -- would
    misattribute every hit. Searches pass only the proteins they hit: building the whole
    proteome's table (82k rows for human) cost more than the search itself on small queries.
    A hit on a protein this index lacks means the search's indexes disagree, which raises."""
    try:
      pnum, pid, name, species, taxon, gene, exist, seqver, geneprio, swiss = rs_metadata(pepidx_path, protein_nums)
    except ValueError as e:   # .pepidx files are reused by name, so a FASTA edited in place leaves stale ones
      raise ValueError(
        f"{e}. This proteome's .pepidx files in {self.preprocessed_files_path} were most likely built "
        f"from different versions of {self.proteome_file}; delete them to rebuild."
      ) from e
    str_cols = ['Protein ID','Protein Name','Species','Taxon ID','Gene',
                'Protein Existence Level','Sequence Version','Gene Priority','SwissProt Reviewed']
    m = pl.DataFrame({
      'protein_num': pl.Series(pnum, dtype=pl.UInt32),
      'Protein ID': pid, 'Protein Name': name, 'Species': species, 'Taxon ID': taxon,
      'Gene': gene, 'Protein Existence Level': exist, 'Sequence Version': seqver,
      'Gene Priority': geneprio, 'SwissProt Reviewed': swiss,
    }, schema_overrides={c: pl.Utf8 for c in str_cols})   # String even when nothing was hit
    m = m.with_columns(
      pl.when(pl.col(c) == '').then(None).otherwise(pl.col(c)).alias(c) for c in str_cols
    )
    return m.with_columns([
      pl.col('Protein Existence Level').cast(pl.Int64),
      pl.col('Sequence Version').cast(pl.Int64),
      pl.col('Gene Priority').cast(pl.Int64),
      pl.col('SwissProt Reviewed').cast(pl.Int64).cast(pl.Boolean),
    ])

  def _final_columns(self, is_indels):
    """Column order of the results. Indel search reports Indels / Indel Positions and
    every other mode Mismatches / Mutated Positions, never both -- an empty twin column
    would imply a search that never ran."""
    edit_col = 'Indels' if is_indels else 'Mismatches'
    pos_col = 'Indel Positions' if is_indels else 'Mutated Positions'
    return [
      'Query ID','Query Sequence','Matched Sequence','Protein ID','Protein Name','Species',
      'Taxon ID','Gene', edit_col, pos_col,'Index start','Index end',
      'Protein Existence Level','Gene Priority','SwissProt Reviewed',
    ]

  def _to_dataframe(self, cols, pepidx_path, is_indels=False):
    """Build the results frame from the columnar Rust output, joining protein metadata
    from `pepidx_path` once instead of cloning it into every hit row. Only the proteins
    the search hit are read, and the join keeps the hits' row order. rs_indel_match
    reuses the mismatch slot for its edit count, so only the column name differs by
    mode; indel positions are derived here from (query, matched) rather than in Rust,
    once per distinct pair since a matched peptide often recurs across proteins, and
    miss rows stay null."""
    qid, qseq, matched, pnum, mm, mutated, istart, iend = cols

    edit_col = 'Indels' if is_indels else 'Mismatches'
    pos_col = 'Indel Positions' if is_indels else 'Mutated Positions'
    final_columns = self._final_columns(is_indels)

    if not qid:
      schema = {c: pl.Utf8 for c in final_columns}
      for c in (edit_col, 'Index start', 'Index end', 'Protein Existence Level', 'Gene Priority'):
        schema[c] = pl.Int64
      schema['SwissProt Reviewed'] = pl.Boolean
      return pl.DataFrame(schema=schema)

    if is_indels:
      annotations = {}
      positions = []
      for q, m in zip(qseq, matched):
        if m is not None and (q, m) not in annotations:
          annotations[(q, m)] = format_indel_positions(q, m)
        positions.append(annotations[(q, m)] if m is not None else None)
    else:
      positions = mutated

    base = pl.DataFrame({
      'Query ID': qid,
      'Query Sequence': qseq,
      'Matched Sequence': matched,
      'protein_num': pl.Series(pnum, dtype=pl.UInt32),
      edit_col: pl.Series(mm, dtype=pl.Int64),
      pos_col: positions,
      'Index start': pl.Series(istart, dtype=pl.Int64),
      'Index end': pl.Series(iend, dtype=pl.Int64),
    })

    hit_proteins = base.get_column('protein_num').drop_nulls().unique().to_list()
    metadata = self._metadata_table(pepidx_path, hit_proteins)
    df = base.join(metadata, on='protein_num', how='left', maintain_order='left').drop('protein_num')

    if self.sequence_version:
      df = df.with_columns(
        pl.when(pl.col('Sequence Version').is_not_null())
          .then(pl.col('Protein ID') + '.' + pl.col('Sequence Version').cast(pl.Utf8))
          .otherwise(pl.col('Protein ID'))
          .alias('Protein ID')
      )

    return df.drop('Sequence Version').select(final_columns)
