"""Brute-force baseline: the exhaustive, provably-complete reference for indel search.

It wraps the COMMITTED oracles (never a reimplementation):
  - pepmatch/tests/indel_brute_force.py   for 1 indel  (the oracle the expected files were built with)
  - pepmatch/tests/indel2_brute_force.py  for 2 indels (generic, homogeneous)
so its output reproduces the ground truth by construction.

A naive all-proteins scan is infeasible, so it uses a SOUND pigeonhole prefilter: split each
query into (indels + 1) disjoint pieces; a <=n-indel alignment perturbs at most n pieces, so at
least one piece survives verbatim in any matching window. Candidates are the proteins containing
one of those pieces -- never a false negative. tests/ is not an importable package, so the oracles
are loaded by file path.
"""
import bisect
import importlib.util
import os
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

_TESTS = Path(__file__).resolve().parents[2] / 'pepmatch' / 'tests'


def _load(name):
  spec = importlib.util.spec_from_file_location(name, _TESTS / f'{name}.py')
  mod = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(mod)
  return mod.brute_force_search


_ORACLE_1 = _load('indel_brute_force')    # brute_force_search(query, protein)
_ORACLE_2 = _load('indel2_brute_force')   # brute_force_search(query, protein, max_indels)


def _oracle(query, protein, n):
  if n == 1:
    return _ORACLE_1(query, protein)
  return _ORACLE_2(query, protein, n)


# Fork-shared proteome (Linux copy-on-write) so Pool workers need no pickling.
_ACCS, _SEQS, _STARTS, _LENS, _CONCAT, _N = [], [], [], [], '', 1


def _proteins_containing(sub):
  hits = set()
  if not sub:
    return hits
  i = _CONCAT.find(sub)
  while i != -1:
    j = bisect.bisect_right(_STARTS, i) - 1
    if i + len(sub) <= _STARTS[j] + _LENS[j]:   # fully inside protein j, not across the separator
      hits.add(j)
    i = _CONCAT.find(sub, i + 1)
  return hits


def _pieces(query, k):
  # The pigeonhole guarantee needs k NON-EMPTY pieces: n indels damage at most n of
  # them, so at least one survives verbatim. If len(query) < k the step collapses to 0
  # and every piece but the last is empty, leaving a single effective bucket -- the
  # guarantee degrades silently and true matches are dropped with no error, corrupting
  # the ground truth this baseline exists to provide. Unreachable for the shipped
  # 1-indel datasets (min query length 8); guarded so it fails loud if that changes.
  if len(query) < k:
    raise ValueError(
      f'Query {query!r} (length {len(query)}) is shorter than the {k} pigeonhole '
      f'pieces required for {k - 1} indel(s); the prefilter would be unsound.'
    )
  step = len(query) // k
  cuts = [i * step for i in range(k)] + [len(query)]
  return [query[cuts[i]:cuts[i + 1]] for i in range(k)]


def _search_one(query):
  # The proteome reaches Pool workers through fork's copy-on-write snapshot. Under a
  # 'spawn' start method the workers would re-import this module with empty globals and
  # every query would return no candidates -- 0% recall reported as a clean success.
  # Fail loud instead.
  if not _CONCAT:
    raise RuntimeError(
      'Proteome globals are empty in this worker: preprocess_proteome did not run here '
      '(non-fork multiprocessing start method?). Refusing to report a silent 0%.'
    )
  cands = set()
  for piece in _pieces(query, _N + 1):
    cands |= _proteins_containing(piece)
  rows = []
  for j in cands:
    for start, matched in _oracle(query, _SEQS[j], _N):
      if len(matched) != len(query):   # indel matches only (drop exact)
        rows.append((query, matched, _ACCS[j], start + 1))
  return rows


class Benchmarker:
  # search() reads the proteome structures that preprocess_proteome() builds in THIS
  # process's memory. Unlike PEPMatch's .pepidx or the aligners' on-disk databases,
  # there is no artifact to hand from one process to the next, so preprocessing and
  # search memory cannot be measured independently. Declared rather than inferred so
  # the harness reports "inseparable" instead of silently producing a wrong split.
  search_requires_preprocess = True

  def __init__(
    self, benchmark: str, query: str, proteome: str, lengths: list, max_mismatches: int,
    method_parameters: dict, indels: int = 0
  ):
    if indels <= 0:
      raise ValueError('Brute Force is an indel-only baseline.')
    self.query = str(query)
    self.proteome = str(proteome)
    self.indels = indels
    self.threads = int(method_parameters.get('threads', os.cpu_count() or 1))

  def __str__(self):
    return 'Brute Force'

  def preprocess_proteome(self):
    # Brute Force is the naive exhaustive reference: unlike the aligners (which build a
    # persisted database) it has no reusable index, so it reports N/A for preprocessing.
    # The in-memory prefilter is built lazily in search() instead -- see _build_index.
    raise TypeError('Brute Force builds its index at search time (no proteome preprocessing).')

  def _build_index(self):
    # Parse the proteome and build the pigeonhole prefilter in memory. Charged to search
    # time, and rebuilt every run (it is not a persisted artifact) -- which is also why
    # the memory table reports Brute Force's phases as inseparable.
    global _ACCS, _SEQS, _STARTS, _LENS, _CONCAT, _N
    _N = self.indels
    accs, seqs, header, chunks = [], [], None, []
    with open(self.proteome) as f:
      for line in f:
        line = line.rstrip('\n')
        if line.startswith('>'):
          if header is not None:
            accs.append(header)
            seqs.append(''.join(chunks))
          parts = line[1:].split('|')
          header = parts[1] if len(parts) >= 2 else line[1:].split()[0]
          chunks = []
        else:
          chunks.append(line)
      if header is not None:
        accs.append(header)
        seqs.append(''.join(chunks))
    starts, lens, cat, pos = [], [], [], 0
    for s in seqs:
      starts.append(pos)
      lens.append(len(s))
      cat.append(s)
      cat.append('\x00')     # separator: never an amino acid, blocks cross-protein matches
      pos += len(s) + 1
    _ACCS, _SEQS, _STARTS, _LENS, _CONCAT = accs, seqs, starts, lens, ''.join(cat)

  def preprocess_query(self):
    raise TypeError('Brute Force does not preprocess queries.')

  def search(self):
    self._build_index()   # build the prefilter here so its cost counts as search time
    queries = []
    with open(self.query) as f:
      for line in f:
        line = line.strip()
        if line and not line.startswith('>'):
          queries.append(line)
    queries = sorted(set(queries))

    if self.threads > 1:
      with Pool(self.threads) as pool:
        results = pool.map(_search_one, queries)
    else:
      results = [_search_one(q) for q in queries]

    rows = [row for sub in results for row in sub]
    return pd.DataFrame(rows, columns=['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start'])
