"""Naive brute force: the TRUE exhaustive baseline, with no prefilter of any kind.

Why this module exists (2026-08-31)
-----------------------------------
`methods/brute_force.py` is a *complete-recall reference*, not a naive baseline: it
narrows candidates with a pigeonhole prefilter (split the query into indels+1 disjoint
pieces, keep proteins containing one of them) before verifying with the oracle. That
prefilter is PEPMatch's own core trick, so "PEPMatch vs Brute Force" was really measuring
Rust + a persisted index against Python + a rebuilt filter -- the two methods shared their
candidate-generation strategy. A reader who sees "brute force" expects the dumb algorithm,
and that is not what they were being shown.

This module is the dumb algorithm. For every query it sweeps EVERY protein in the
proteome and asks the committed oracle whether the query aligns anywhere in it. No
seeding, no k-mers, no pigeonhole, no index. It shares no candidate-generation logic with
PEPMatch whatsoever.

Equivalence to "combinatorial expansion"
----------------------------------------
The baseline is often described as: expand each query into every string reachable by <=n
indels, then look each one up in the proteome. That is the same computation done in the
wasteful direction. For insertions the expansion is over the 20-letter amino-acid
alphabet, so a length-25 query with 2 insertions generates ~300 gap-position choices x 400
residue pairs ~= 120,000 variant strings -- and the overwhelming majority are strings that
appear nowhere in the proteome. Sweeping the proteome's OWN windows instead tests only
substrings that actually exist, and returns an identical answer set: both enumerate exactly
the windows W such that W reduces to the query under <=n homogeneous indels.

So this is the combinatorial-expansion baseline, evaluated in the direction that does not
waste 99.99% of its work on variants that cannot match. It is still O(proteome) per query,
which is the point: that is what brute force costs.

Correctness
-----------
Wraps the COMMITTED oracles (never a reimplementation), the same ones `brute_force.py`
uses and the same ones the synthetic ground truth was built with:
  - pepmatch/tests/indel_brute_force.py   (1 indel)
  - pepmatch/tests/indel2_brute_force.py  (2 indels, homogeneous)

This is byte-identical in OUTPUT to `brute_force.py` -- the prefilter it drops was sound
(zero false negatives), so removing it can only change runtime, never the row set. That
equivalence is not an assumption: `scripts/verify-naive-equivalence.py` proves it on real
datasets, and it is the same no-prefilter path the generator's Phase-4 audit already used
to certify the ground truth.
"""
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


# ONE algorithm for every indel count. The generic oracle (indel2_brute_force) takes
# max_indels and does a SINGLE positional sweep of the protein with an early-exit aligner,
# so a length-L query costs one pass per protein regardless of n. The 1-indel-specific
# oracle (indel_brute_force) instead rebuilds the query for each edit position and re-scans
# the whole protein once per position (~2L passes), which made the 1-indel baseline ~5x
# SLOWER than 2-indel purely as an implementation artifact -- an incoherent cost curve for
# a "brute force" (more edits should never be cheaper). Using the generic oracle for both
# removes that. Proven to return identical rows for n=1 vs the committed 1-indel oracle
# (scripts/verify-naive-equivalence.py --check-oracle-equivalence), and the sweep re-scores
# recall against the committed ground truth every run, so any divergence would be loud, not
# silent. The 1-indel oracle stays the ground-truth builder in the generator; only this
# timing baseline is unified.
_ORACLE = _load('indel2_brute_force')     # brute_force_search(query, protein, max_indels)


def _oracle(query, protein, n):
  return _ORACLE(query, protein, n)


# Fork-shared proteome (Linux copy-on-write) so Pool workers need no pickling.
# Deliberately NO _CONCAT / _STARTS here: those exist in brute_force.py only to support
# the substring prefilter this module refuses to use.
_ACCS, _SEQS, _N = [], [], 1


def _search_one(query):
  # The proteome reaches Pool workers through fork's copy-on-write snapshot. Under a
  # 'spawn' start method the workers would re-import this module with empty globals and
  # every query would return no candidates -- 0% recall reported as a clean success.
  # Fail loud instead. (Same guard as brute_force.py, same failure mode.)
  if not _SEQS:
    raise RuntimeError(
      'Proteome globals are empty in this worker: the proteome was not loaded here '
      '(non-fork multiprocessing start method?). Refusing to report a silent 0%.'
    )
  qlen = len(query)
  # The ONLY skip in this module, and it is not a prefilter: a protein shorter than the
  # shortest possible match window cannot contain one at any offset. It inspects a length,
  # never the residues, so it cannot discriminate between proteins that do and do not share
  # content with the query -- no candidate narrowing happens here.
  min_window = qlen - _N
  rows = []
  for j, seq in enumerate(_SEQS):
    if len(seq) < min_window:
      continue
    for start, matched in _oracle(query, seq, _N):
      if len(matched) != qlen:   # indel matches only (drop exact) -- matches brute_force.py
        rows.append((query, matched, _ACCS[j], start + 1))
  return rows


class Benchmarker:
  # search() reads the proteome that _load_proteome() builds in THIS process's memory.
  # There is no persisted artifact to hand from one process to the next, so preprocessing
  # and search memory cannot be measured independently. Declared rather than inferred so
  # the harness reports "inseparable" instead of silently producing a wrong split.
  search_requires_preprocess = True

  def __init__(
    self, benchmark: str, query: str, proteome: str, lengths: list, max_mismatches: int,
    method_parameters: dict, indels: int = 0
  ):
    if indels <= 0:
      raise ValueError('Naive Brute Force is an indel-only baseline.')
    self.query = str(query)
    self.proteome = str(proteome)
    self.indels = indels
    self.threads = int(method_parameters.get('threads', os.cpu_count() or 1))

  def __str__(self):
    return 'Brute Force (naive)'

  def preprocess_proteome(self):
    # Nothing is precomputed: there is no index, not even the in-memory substring table
    # brute_force.py builds. Parsing the FASTA is charged to search time, as it is there.
    raise TypeError('Naive Brute Force builds nothing (no proteome preprocessing).')

  def _load_proteome(self):
    # Same accession parse as methods/brute_force.py (`line[1:].split('|')[1]`), so the
    # Protein ID column is identical between the two baselines and their outputs can be
    # diffed row-for-row.
    global _ACCS, _SEQS, _N
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
    _ACCS, _SEQS = accs, seqs

  def preprocess_query(self):
    raise TypeError('Naive Brute Force does not preprocess queries.')

  def search(self):
    self._load_proteome()   # charged to search time (no persisted artifact)
    queries = []
    with open(self.query) as f:
      for line in f:
        line = line.strip()
        if line and not line.startswith('>'):
          queries.append(line)
    queries = sorted(set(queries))

    if self.threads > 1:
      # chunksize=1: per-query cost varies by orders of magnitude with query length (a
      # length-8 2-indel query is ~65 core-s, a length-25 is ~0.001), so the default
      # chunking would hand one worker a contiguous block of short queries and leave the
      # run tailing on a single core.
      with Pool(self.threads) as pool:
        results = pool.map(_search_one, queries, chunksize=1)
    else:
      results = [_search_one(q) for q in queries]

    rows = [row for sub in results for row in sub]
    return pd.DataFrame(rows, columns=['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start'])
