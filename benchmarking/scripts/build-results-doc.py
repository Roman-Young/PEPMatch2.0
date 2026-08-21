#!/usr/bin/env python3
"""Build the paper-ready results document from the benchmark result TSVs.

Reads the four tables benchmarking.py produces and writes, into this directory:
  indel-benchmarking-results.xlsx   one sheet per table, with its footnotes below it
  RESULTS.md                        all tables + footnotes as Markdown

The numbers are read from the TSVs, never hand-typed, so the document cannot drift from
the actual results. Run it wherever the TSVs live -- i.e. on the cluster once the timing
and memory jobs have finished:

    python3 benchmarking/scripts/build-results-doc.py

Missing TSVs are skipped with a warning (e.g. running before the memory job completes),
so the document always reflects whatever results exist so far.
"""
import sys
from pathlib import Path

import pandas as pd

# This script lives in benchmarking/scripts/; the result TSVs it reads and the outputs
# it writes both live in benchmarking/results/.
RESULTS = Path(__file__).resolve().parent.parent / 'results'
XLSX_OUT = RESULTS / 'indel-benchmarking-results.xlsx'
MD_OUT = RESULTS / 'RESULTS.md'

# Each table: (tsv, sheet title, source->display columns, per-column decimal places).
TIMING_COLS = {
  'Method': 'Method',
  'Proteome Preprocessing (s)': 'Proteome preprocessing (s)',
  'Searching (s)': 'Search (s)',
  'Total (s)': 'Total (s)',
  'Recall (%)': 'Recall (%)',
}
TIMING_DECIMALS = {'Proteome preprocessing (s)': 1, 'Search (s)': 1, 'Total (s)': 1, 'Recall (%)': 1}

MEMORY_COLS = {
  'Method': 'Method',
  'Preprocessing (MB)': 'Preprocessing (MB)',
  'Search (MB)': 'Search (MB)',
  'Peak total (MB)': 'Peak total (MB)',
}
MEMORY_DECIMALS = {'Preprocessing (MB)': 0, 'Search (MB)': 0, 'Peak total (MB)': 0}

TABLES = [
  ('cosmic_indel_benchmarking.tsv', 'Timing – COSMIC', TIMING_COLS, TIMING_DECIMALS, 'timing',
   'Timing — COSMIC (3,510 queries, 9-mers)'),
  ('cedar_indel_benchmarking.tsv', 'Timing – CEDAR', TIMING_COLS, TIMING_DECIMALS, 'timing',
   'Timing — CEDAR (50 queries, 8–26 aa)'),
  ('cosmic_indel_memory.tsv', 'Memory – COSMIC', MEMORY_COLS, MEMORY_DECIMALS, 'memory',
   'Memory — COSMIC (MB)'),
  ('cedar_indel_memory.tsv', 'Memory – CEDAR', MEMORY_COLS, MEMORY_DECIMALS, 'memory',
   'Memory — CEDAR (MB)'),
]

TIMING_NOTES = [
  'Recall is measured against an exhaustive brute-force oracle, independently validated.',
  'All alignment tools (BLAST, DIAMOND, MMseqs2) were run at maximum sensitivity with '
  'low-complexity masking disabled — their most accurate, and slowest, configuration.',
  'BLAST is reported in two modes: "short" (blastp-short: PAM30 matrix and a short word '
  'size, appropriate for peptide-length queries) and "default" (blastp).',
  'COSMIC queries are uniform 9-mers (~88% deletions, ~12% insertions); CEDAR queries '
  'range 8–26 aa.',
  'Query preprocessing is omitted: no method preprocesses queries for indel search.',
  'All tools were pinned to the same thread count on one compute node.',
]
MEMORY_NOTES = [
  'Values are peak resident set size (RSS) in MB, excluding the shared Python interpreter '
  'baseline, so in-process tools (PEPMatch, Brute Force) and standalone binaries (BLAST, '
  'DIAMOND, MMseqs2) are directly comparable.',
  'Preprocessing is a one-time cost, amortised across all subsequent searches; Search is '
  'paid per query.',
  'Peak total is measured end-to-end, not derived: the true peak can exceed '
  'max(preprocessing, search) when preprocessing memory is not released before search.',
  'Brute Force reports "inseparable": its proteome resides in process memory with no '
  'on-disk artifact, so preprocessing and search cannot be measured independently — only '
  'the end-to-end peak is reportable.',
  'Alignment tools were run at maximum sensitivity with masking disabled (as in timing).',
]


def fmt(value, decimals):
  """Round a numeric cell; pass non-numeric cells (inseparable / N/A) through unchanged."""
  try:
    return f'{float(value):.{decimals}f}'
  except (ValueError, TypeError):
    return str(value)


def load_table(tsv, colmap, decimals):
  path = RESULTS / tsv
  if not path.exists():
    print(f'  skip: {tsv} not found (has that job finished?)')
    return None
  df = pd.read_csv(path, sep='\t', dtype=str)
  missing = [c for c in colmap if c not in df.columns]
  if missing:
    print(f'  skip: {tsv} missing columns {missing}')
    return None
  out = df[list(colmap)].rename(columns=colmap)
  for col, nd in decimals.items():
    out[col] = out[col].map(lambda v, _nd=nd: fmt(v, _nd))
  return out


def write_xlsx(loaded):
  try:
    import xlsxwriter  # noqa: F401
  except ImportError:
    print('  xlsxwriter not installed; skipping .xlsx (Markdown still written). '
          'pip install xlsxwriter')
    return False

  with pd.ExcelWriter(XLSX_OUT, engine='xlsxwriter') as writer:
    book = writer.book
    fmt_title = book.add_format({'bold': True, 'font_size': 13})
    fmt_head = book.add_format({'bold': True, 'bottom': 1, 'align': 'left'})
    fmt_cell = book.add_format({'align': 'left'})
    fmt_pep = book.add_format({'bold': True, 'align': 'left'})
    fmt_note = book.add_format({'italic': True, 'font_size': 9, 'text_wrap': True,
                                'valign': 'top'})

    for sheet, title, df, notes in loaded:
      ws = book.add_worksheet(sheet)
      ws.write(0, 0, title, fmt_title)
      # header
      for c, col in enumerate(df.columns):
        ws.write(2, c, col, fmt_head)
      # body (bold the PEPMatch row)
      for r, (_, row) in enumerate(df.iterrows(), start=3):
        rowfmt = fmt_pep if str(row['Method']).startswith('PEPMatch') else fmt_cell
        for c, col in enumerate(df.columns):
          ws.write(r, c, row[col], rowfmt)
      # footnotes below the table
      note_row = 3 + len(df) + 1
      for i, note in enumerate(notes):
        ws.write(note_row + i, 0, f'• {note}', fmt_note)
      # column widths
      ws.set_column(0, 0, 18)
      ws.set_column(1, len(df.columns) - 1, 22)
  print(f'  wrote {XLSX_OUT}')
  return True


def write_markdown(loaded):
  lines = ['# PEPMatch 1-indel benchmark — results', '']
  for _sheet, title, df, notes in loaded:
    lines.append(f'## {title}')
    lines.append('')
    lines.append('| ' + ' | '.join(df.columns) + ' |')
    lines.append('|' + '|'.join('---' for _ in df.columns) + '|')
    for _, row in df.iterrows():
      cells = [f'**{row[c]}**' if str(row['Method']).startswith('PEPMatch') else str(row[c])
               for c in df.columns]
      lines.append('| ' + ' | '.join(cells) + ' |')
    lines.append('')
    for note in notes:
      lines.append(f'- {note}')
    lines.append('')
  MD_OUT.write_text('\n'.join(lines))
  print(f'  wrote {MD_OUT}')
  return '\n'.join(lines)


def main():
  loaded = []
  for tsv, sheet, colmap, decimals, kind, title in TABLES:
    df = load_table(tsv, colmap, decimals)
    if df is None:
      continue
    notes = TIMING_NOTES if kind == 'timing' else MEMORY_NOTES
    loaded.append((sheet, title, df, notes))

  if not loaded:
    print('No result TSVs found — nothing to build.')
    return 1

  write_xlsx(loaded)
  md = write_markdown(loaded)
  print('\n' + md)
  return 0


if __name__ == '__main__':
  sys.exit(main())
