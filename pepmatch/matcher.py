import os
import re
import polars as pl
from pathlib import Path
from Bio import SeqIO
from .preprocessor import MAX_INDEX_K, build_index, built_from, check_index
from ._rs import rs_match, rs_discontinuous, rs_metadata, rs_match_counts, rs_indel_match

VALID_OUTPUT_FORMATS = ['dataframe', 'polars', 'pandas', 'csv', 'tsv', 'xlsx', 'json']
DATAFRAME_OUTPUT_FORMATS = frozenset({'dataframe', 'polars', 'pandas'})
FASTA_EXTENSIONS = {
  '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa', '.frn'
}

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


def format_indel_positions(query, matched):
  """Render the Indel Positions column, e.g. 'd: A[6]', 'd: AA[2,4]', 'd: C[3], d: E[5]',
  or '[]' for an exact match. A 1-based range collapses when the placement is unambiguous."""
  edits = _indel_edits(query, matched)
  if not edits:
    return '[]'
  return ', '.join(
    f'{kind}: {residues}[{low}]' if low == high else f'{kind}: {residues}[{low},{high}]'
    for kind, residues, low, high in edits
  )


def k_for(length, edits):
  """The k a peptide of `length` is searched at with `edits` mismatches or indels: the
  largest that still splits it into edits + 1 disjoint k-mers, so at least one of them is
  untouched in any match (pigeonhole) and recall is complete. Floored at 2, the smallest
  index, where a peptide shorter than 2 * (edits + 1) loses that guarantee, and capped at
  MAX_INDEX_K, the longest k-mer an index can store (a smaller k only adds seeds)."""
  return min(MAX_INDEX_K, max(2, length // (edits + 1)))


def best_match_budgets(length):
  """The mismatch budgets best match tries for a peptide of `length`, in order. While the
  recall guarantee holds, each k the peptide can use is searched once, at the largest
  budget that k guarantees; past it, one more mismatch at a time."""
  budgets = []
  m = 0
  while length >= 2 * (m + 1):
    m = length // k_for(length, m) - 1
    budgets.append(m)
    m += 1
  return budgets + list(range(m, length + 1))


class Matcher:
  """Searches query peptides against a preprocessed proteome index
  and returns matches as a DataFrame (Polars by default, or pandas) or output file."""

  def __init__(
    self,
    query,
    proteome_file,
    max_mismatches=0,
    max_indels=0,
    preprocessed_files_path='.',
    best_match=False,
    output_format='dataframe',
    output_name='',
    sequence_version=True,
    counts_only=False,
    proteome_name='',
  ):

    if max_mismatches < 0 or max_indels < 0:
      raise ValueError('max_mismatches and max_indels must be >= 0.')

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

    self.proteome_file = str(proteome_file)
    # A custom name must match the one given to Preprocessor, which names the index files.
    self.proteome_name = proteome_name or Path(proteome_file).stem
    self.max_mismatches = max_mismatches
    self.preprocessed_files_path = preprocessed_files_path
    self.best_match = best_match
    self.output_format = output_format
    self.sequence_version = sequence_version
    self.counts_only = counts_only
    self.output_name = output_name or 'PEPMatch_results'

    self.query = self._parse_query(query)
    self.discontinuous_epitopes = self._find_discontinuous_epitopes()
    if self.max_indels > 0 and self.discontinuous_epitopes:
      raise ValueError('max_indels and discontinuous epitopes are not supported together.')
    self.query = self._clean_query()
    self.unsearchable = self._set_aside_unsearchable()
    if not self.query and not self.discontinuous_epitopes and not self.unsearchable:
      raise ValueError('Query is empty.')

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
    """Pick out the queries written as residue/position lists, e.g. 'N1, A3, L5' or
    'N1,A3,L5'. No peptide contains a comma, so either spacing is unambiguous. An epitope
    with a negative position is noted for _set_aside_unsearchable."""
    discontinuous_epitopes = {}
    self._negative_epitopes = set()
    for query_id, peptide in self.query:
      if len(peptide) < 2 or peptide[1] not in '-0123456789':   # an epitope starts with a residue and its position
        continue
      try:
        epitope = [(x[0], int(x[1:])) for x in re.split(r'\s*,\s*', peptide.strip())]
      except (ValueError, IndexError):
        continue
      if any(position < 0 for _, position in epitope):
        self._negative_epitopes.add((query_id, peptide))
        continue
      discontinuous_epitopes[query_id] = epitope
    return discontinuous_epitopes

  def _clean_query(self):
    """Drop the discontinuous epitopes from the linear query set."""
    discontinuous_ids = set(self.discontinuous_epitopes.keys())
    return [
      (qid, seq) for qid, seq in self.query
      if qid not in discontinuous_ids
    ]

  def _set_aside_unsearchable(self):
    """Remove the linear queries that cannot be searched -- empty, containing whitespace, or
    an epitope with a negative position -- and warn about them. They are reported as misses
    instead of stopping the run, and they take no part in choosing k, so one bad entry in a
    large file cannot drag every other query down to a small k."""
    searchable, unsearchable = [], []
    for qid, seq in self.query:
      if (qid, seq) in self._negative_epitopes:
        unsearchable.append((qid, seq, 'lists a negative position (positions are 1-based)'))
      elif seq.split() != [seq]:   # empty, or whitespace anywhere
        unsearchable.append((qid, seq, 'is empty or contains whitespace'))
      else:
        searchable.append((qid, seq))
    if unsearchable:
      print(f'Warning: {len(unsearchable)} query(s) cannot be searched and are reported as misses:')
      for qid, seq, reason in unsearchable[:10]:
        print(f'  {qid} ({seq!r}) {reason}')
      if len(unsearchable) > 10:
        print(f'  ... and {len(unsearchable) - 10} more')
    self.query = searchable
    return unsearchable

  def _miss_rows(self, queries, is_indels=False):
    """A results frame of miss rows for (id, sequence) `queries`, as the engine reports them."""
    empty = self._to_dataframe(tuple([] for _ in range(8)), None, is_indels)
    misses = pl.DataFrame(
      {'Query ID': [q for q, _ in queries], 'Query Sequence': [s for _, s in queries]},
      schema={'Query ID': pl.Utf8, 'Query Sequence': pl.Utf8},
    )
    if not is_indels:
      misses = misses.with_columns(pl.lit('[]').alias('Mutated Positions'))
    return pl.concat([empty, misses], how='diagonal_relaxed').select(empty.columns)

  def _add_unsearchable_misses(self, df):
    """`df` with a miss row for each query _set_aside_unsearchable removed."""
    if not self.unsearchable:
      return df
    misses = self._miss_rows([(q, s) for q, s, _ in self.unsearchable], self.max_indels > 0)
    return pl.concat([df, misses], how='vertical_relaxed') if df.width else misses

  def match(self):
    """Run the search, returning a DataFrame or writing the results to `output_name`."""
    linear_df = pl.DataFrame()
    discontinuous_df = pl.DataFrame()

    if self.counts_only:
      return self._finish(self._count())

    if self.query:
      linear_df = self.best_match_search() if self.best_match else self._search()
    linear_df = self._add_unsearchable_misses(linear_df)

    if self.discontinuous_epitopes:
      pepidx_path = self._index(2)
      epitopes = [
        (qid, residues) for qid, residues in self.discontinuous_epitopes.items()
      ]
      results = rs_discontinuous(pepidx_path, epitopes, self.max_mismatches)
      discontinuous_df = self._to_dataframe(results, pepidx_path)

    dfs = [d for d in [linear_df, discontinuous_df] if d.height > 0]
    df = pl.concat(dfs, how="vertical_relaxed") if dfs else linear_df
    return self._finish(df)

  def _finish(self, df):
    """Return the results as a DataFrame or write them to `output_name`."""
    if self.output_format in DATAFRAME_OUTPUT_FORMATS:
      return _return_dataframe(df, self.output_format)
    output_matches(df, self.output_format, self.output_name)

  def _groups(self, queries, edits):
    """`queries` split by the k each is searched at, smallest k first."""
    groups = {}
    for query in queries:
      groups.setdefault(k_for(len(query[1]), edits), []).append(query)
    return sorted(groups.items())

  def _search_groups(self, queries, edits, is_indels=False):
    """Search each peptide at its own k and return one results frame per k. Each frame is
    labeled from the index that found its hits: protein numbers are that index's own FASTA
    ordinals, so another group's index could name the wrong protein."""
    frames = []
    for k, group in self._groups(queries, edits):
      pepidx_path = self._index(k)
      print(f"Searching {len(group)} peptides against {self.proteome_name} "
            f"(k={k}, {'max_indels' if is_indels else 'max_mismatches'}={edits})...")
      if is_indels:
        results = rs_indel_match(pepidx_path, group, edits)
      else:
        results = rs_match(pepidx_path, group, k, edits)
      frames.append(self._to_dataframe(results, pepidx_path, is_indels))
    return frames

  def _in_query_order(self, frames):
    """The groups' frames as one, rows in the order the queries were given. Query IDs are
    still positions. vertical_relaxed lets an all-miss group's null columns stack under
    another group's String columns."""
    df = pl.concat(frames, how='vertical_relaxed')
    if len(frames) > 1:
      df = df.sort(pl.col('Query ID').cast(pl.Int64), maintain_order=True)
    return df

  def _with_query_ids(self, df):
    """`df` with the positional Query IDs replaced by the real ones."""
    ids = {str(i): qid for i, (qid, _) in enumerate(self.query)}
    return df.with_columns(pl.col('Query ID').replace_strict(ids, return_dtype=pl.Utf8))

  def _search(self):
    """Mismatch or indel search of every linear query, each at its own k.

    max_indels=2 is homogeneous: two insertions or two deletions, never one of each. The
    text columns are cast to String so an indel batch that hits nothing has the same schema
    as one that does."""
    is_indels = self.max_indels > 0
    edits = self.max_indels if is_indels else self.max_mismatches
    self._recall_warning(edits)
    df = self._with_query_ids(self._in_query_order(self._search_groups(self._positional_query(), edits, is_indels)))
    if is_indels:
      df = df.with_columns(pl.col('Matched Sequence', 'Indel Positions').cast(pl.Utf8))
    return df

  def _count(self):
    """Hits per (peptide, mismatch level), each peptide counted at its own k. One row per
    level with a hit; memory is O(unique queries), independent of the number of hits."""
    schema = {'Query ID': pl.Utf8, 'Query Sequence': pl.Utf8, 'Mismatches': pl.Int64, 'Count': pl.UInt64}
    if not self.query:   # counts are kept only for linear peptides
      return pl.DataFrame(schema=schema)
    edits = self.max_mismatches
    self._recall_warning(edits)
    frames = []
    for k, group in self._groups(self._positional_query(), edits):
      pepidx_path = self._index(k)
      print(f"Counting {len(group)} peptides against {self.proteome_name} (k={k}, max_mismatches={edits})...")
      qid, qseq, mm, count = rs_match_counts(pepidx_path, group, k, edits)
      frames.append(pl.DataFrame({'Query ID': qid, 'Query Sequence': qseq, 'Mismatches': mm, 'Count': count}, schema=schema))
    return self._with_query_ids(self._in_query_order(frames))

  def _recall_warning(self, edits):
    """Warn about peptides too short for the seed guarantee: below 2 * (edits + 1) residues
    there are not edits + 1 disjoint 2-mers, so a real match may have an edit in every seed."""
    unguaranteed = sorted({len(s) for _, s in self.query if len(s) < 2 * (edits + 1)})
    if unguaranteed:
      print(f"edits={edits}: complete recall is not guaranteed for query peptides shorter "
            f"than {2 * (edits + 1)} residues (lengths: {unguaranteed}); a match may exist "
            f"but be missed.")

  def _pepidx_path(self, k):
    return os.path.join(
      self.preprocessed_files_path, f'{self.proteome_name}_{k}mers.pepidx'
    )

  def _index(self, k):
    """Path of this proteome's k-mer index after checking it is complete, built if it is
    missing and rebuilt if the FASTA's content changed since (see built_from). An index
    whose FASTA is gone is used as it is."""
    path = self._pepidx_path(k)
    if os.path.isfile(path):
      check_index(path, k)
      if not os.path.isfile(self.proteome_file) or built_from(path, self.proteome_file):
        return path
      print(f"{self.proteome_file} changed since {path} was built; rebuilding it...")
    else:
      print(f"No preprocessed file found for k={k}, building index now (this may take a moment)...")
    build_index(self.proteome_file, k, path)
    check_index(path, k)
    return path

  def _positional_query(self):
    """The query with each peptide's position as its ID. Searches run on positions, since
    FASTA IDs need not be unique, and _with_query_ids puts the real IDs back."""
    return [(str(i), seq) for i, (_, seq) in enumerate(self.query)]

  def best_match_search(self):
    """Closest match per peptide: each is searched with a growing mismatch budget (see
    best_match_budgets), at the k that budget allows, until it has a hit. A round searches
    every peptide still unmatched, so peptides of one length and budget share a search."""
    remaining = self._positional_query()
    budgets = {qid: best_match_budgets(len(seq)) for qid, seq in remaining}
    frames = []
    round_number = 0
    while remaining:
      by_budget, exhausted = {}, []
      for qid, seq in remaining:
        if round_number < len(budgets[qid]):
          by_budget.setdefault(budgets[qid][round_number], []).append((qid, seq))
        else:
          exhausted.append((qid, seq))
      matched = set()
      for max_mm, peptides in sorted(by_budget.items()):
        for frame in self._search_groups(peptides, max_mm):
          hits = frame.filter(pl.col('Matched Sequence').is_not_null())
          frames.append(hits)
          matched.update(hits['Query ID'])
      if exhausted:   # no budget left: reported as misses
        frames.append(self._miss_rows(exhausted))
      remaining = [(qid, seq) for qid, seq in remaining if qid not in matched and round_number < len(budgets[qid])]
      print(f"  -> round {round_number + 1}: {len(remaining)} remaining")
      round_number += 1
    best = self._best_match_filter(pl.concat(frames, how='vertical_relaxed'))
    return self._with_query_ids(best.sort(pl.col('Query ID').cast(pl.Int64)))

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
