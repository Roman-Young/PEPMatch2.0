#!/usr/bin/env python3
"""Gate a synthetic indel dataset before it is trusted or shipped to the cluster.

WHY: the expensive failures in this benchmark are all SILENT. An off-by-one in
`Index start`, an accession carrying a `.2` sequence version, or an exact-match row left
in the expected file do not crash anything -- they produce a clean-looking table where
every tool scores 0%, or where Brute Force mysteriously sits below 100%. Those are
indistinguishable from a real result until someone digs. This script turns each of them
into a loud, specific failure, cheaply, before any cluster time is spent.

Gates, in order of cost:

  FAST (seconds, every row of every size)
    1  coordinate literalism -- Matched Sequence IS proteome[start-1 : start-1+len]
    2  schema -- the 4 join columns present, named exactly, Index start a positive int
    3  Protein ID carries no .SV sequence version
    4  NO exact-match rows (methods/brute_force.py drops exact matches, so an exact row
       in expected is unreachable ground truth and silently caps Brute Force below 100%)
    5  no duplicate join-tuples (recall() joins on exactly those 4 columns)
    6  strict subset nesting, and each subset's expected == the largest filtered to it
    7  50/50 insertion/deletion balance in EVERY nested subset
    8  every edit interior, per the engine's own _indel_placements
    9  FASTA hygiene: consecutive integer ids from 1, unwrapped, distinct sequences,
       every sequence present in expected
   10  manifest sha256s still match the files on disk

  REAL (minutes, --run-tools; on the smallest subset by default)
   11  PEPMatch scores exactly 100.0
   12  Brute Force scores exactly 100.0
       -- both via the SAME benchmarking.recall() and RESULT_COLUMNS the published table
          uses, so a passing gate here means the published number will agree.

Exit code is non-zero on any failure, so it can guard an sbatch preflight.

Usage
  verify-synthetic-dataset.py --prefix synth-1indel
  verify-synthetic-dataset.py --prefix synth-1indel --fast          # skip 11-12
  verify-synthetic-dataset.py --prefix synth-1indel --run-tools-on 100
"""
import argparse
import csv
import hashlib
import json
import random
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BENCH = REPO / 'benchmarking'

sys.path.insert(0, str(BENCH))
csv.field_size_limit(1 << 30)

JOIN_COLUMNS = ['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start']


class Failures:
  def __init__(self):
    self.items = []

  def check(self, ok, gate, detail=''):
    mark = 'PASS' if ok else 'FAIL'
    print(f'  [{mark}] {gate}' + (f' -- {detail}' if detail and not ok else
                                  (f' ({detail})' if detail else '')))
    if not ok:
      self.items.append(f'{gate}: {detail}')
    return ok


def load_proteome(path):
  """Same accession parse as methods/brute_force.py (`line[1:].split('|')[1]`). Using a
  different one here would let a real mismatch pass verification."""
  seqs, header, chunks = {}, None, []
  with open(path) as f:
    for line in f:
      line = line.rstrip('\n')
      if line.startswith('>'):
        if header is not None:
          seqs[header] = ''.join(chunks)
        parts = line[1:].split('|')
        header = parts[1] if len(parts) >= 2 else line[1:].split()[0]
        chunks = []
      else:
        chunks.append(line)
    if header is not None:
      seqs[header] = ''.join(chunks)
  return seqs


def read_tsv(path):
  with open(path, newline='') as f:
    return list(csv.DictReader(f, delimiter='\t'))


def read_fasta(path):
  ids, seqs = [], []
  with open(path) as f:
    for line in f:
      line = line.strip()
      if not line:
        continue
      (ids if line.startswith('>') else seqs).append(line.lstrip('>'))
  return ids, seqs


def sha256_file(path, _buf=1 << 20):
  h = hashlib.sha256()
  with open(path, 'rb') as f:
    while chunk := f.read(_buf):
      h.update(chunk)
  return h.hexdigest()


def main():
  p = argparse.ArgumentParser()
  p.add_argument('--prefix', required=True, help='e.g. synth-1indel')
  p.add_argument('--proteome', default=str(BENCH / 'proteomes' / 'human.fasta'))
  p.add_argument('--fast', action='store_true', help='skip the PEPMatch/Brute Force gates')
  p.add_argument('--run-tools-on', type=int, default=None,
                 help='subset size for gates 11-12 (default: the smallest emitted)')
  p.add_argument('--threads', type=int, default=4)
  args = p.parse_args()

  manifest_path = BENCH / 'expected' / f'{args.prefix}-manifest.json'
  if not manifest_path.exists():
    raise SystemExit(f'no manifest at {manifest_path} -- generate the dataset first')
  manifest = json.loads(manifest_path.read_text())
  indels = manifest['indels']
  files = manifest['files']
  sizes = sorted(files.values(), key=lambda i: i['size'])

  print(f'{"=" * 72}\n  VERIFY {args.prefix}  ({indels}-indel, {len(sizes)} sizes)\n{"=" * 72}')
  f = Failures()

  print('\nloading proteome...')
  proteome = load_proteome(args.proteome)
  print(f'  {len(proteome):,} proteins')

  actual_sha = sha256_file(args.proteome)
  f.check(actual_sha == manifest['proteome']['sha256'],
          'proteome matches the one the dataset was built against',
          'sha256 differs -- coordinates may not correspond' if
          actual_sha != manifest['proteome']['sha256'] else actual_sha[:12])

  largest = sizes[-1]
  largest_rows = read_tsv(BENCH / largest['expected_tsv'])
  largest_by_query = {}
  for r in largest_rows:
    largest_by_query.setdefault(r['Query Sequence'], []).append(r)

  prev_queries = None
  for info in sizes:
    lab = [k for k, v in files.items() if v is info][0]
    print(f'\n--- subset {lab} ({info["size"]:,} queries) ---')
    rows = read_tsv(BENCH / info['expected_tsv'])
    ids, seqs = read_fasta(BENCH / info['queries_fasta'])

    # 2 schema
    missing_cols = [c for c in JOIN_COLUMNS if not rows or c not in rows[0]]
    f.check(not missing_cols, 'join columns present', ','.join(missing_cols))

    bad_start = [r for r in rows if not r['Index start'].isdigit() or int(r['Index start']) < 1]
    f.check(not bad_start, 'Index start is a positive integer',
            f'{len(bad_start)} bad, e.g. {bad_start[0]["Index start"]!r}' if bad_start else '')

    # 3 no sequence version
    versioned = [r for r in rows if '.' in r['Protein ID']]
    f.check(not versioned, 'Protein ID has no .SV suffix',
            f'{len(versioned)} versioned, e.g. {versioned[0]["Protein ID"]}' if versioned else '')

    # 1 coordinate literalism -- the single highest-value check
    bad_coord, missing_prot = [], []
    for r in rows:
      seq = proteome.get(r['Protein ID'])
      if seq is None:
        missing_prot.append(r['Protein ID'])
        continue
      s = int(r['Index start']) - 1
      if seq[s:s + len(r['Matched Sequence'])] != r['Matched Sequence']:
        bad_coord.append(r)
    f.check(not missing_prot, 'every Protein ID exists in the proteome',
            f'{len(missing_prot)} missing, e.g. {missing_prot[0]}' if missing_prot else '')
    f.check(not bad_coord, 'Matched Sequence is the literal proteome slice at Index start',
            (f'{len(bad_coord)} mismatched, e.g. {bad_coord[0]["Query Sequence"]} '
             f'@ {bad_coord[0]["Protein ID"]}:{bad_coord[0]["Index start"]}')
            if bad_coord else f'{len(rows):,} rows')

    # 4 no exact matches
    exact = [r for r in rows if len(r['Matched Sequence']) == len(r['Query Sequence'])]
    f.check(not exact, 'no exact-match rows (brute_force.py drops these)',
            f'{len(exact)} exact rows, e.g. {exact[0]["Query Sequence"]}' if exact else '')

    # length delta must be exactly the indel budget
    bad_delta = [r for r in rows
                 if abs(len(r['Matched Sequence']) - len(r['Query Sequence'])) != indels]
    f.check(not bad_delta, f'every row differs by exactly {indels} residue(s)',
            f'{len(bad_delta)} rows with the wrong delta' if bad_delta else '')

    # Match Type must agree with the length delta. matcher.py derives the engine's own
    # label purely from sign(len(matched) - len(query)), and the natural way to build a
    # deletion-match is to INSERT into the source window -- so a generator that reasons
    # from the edit it applied rather than the resulting delta inverts every label.
    mislabelled = [
      r for r in rows
      if r.get('Match Type') != ('insertion_match'
                                 if len(r['Matched Sequence']) > len(r['Query Sequence'])
                                 else 'deletion_match')
    ]
    f.check(not mislabelled, 'Match Type agrees with sign(len(matched) - len(query))',
            f'{len(mislabelled)} inverted, e.g. {mislabelled[0]["Query Sequence"]}'
            if mislabelled else '')

    # 5 duplicates on the join key
    keys = [tuple(r[c] for c in JOIN_COLUMNS) for r in rows]
    dupes = [k for k, c in Counter(keys).items() if c > 1]
    f.check(not dupes, 'no duplicate join-tuples',
            f'{len(dupes)} duplicated, e.g. {dupes[0]}' if dupes else '')

    # 9 FASTA hygiene
    f.check(ids == [str(i) for i in range(1, len(seqs) + 1)],
            'FASTA ids are consecutive integers from 1')
    f.check(len(seqs) == info['size'], 'FASTA record count matches the manifest',
            f'{len(seqs)} vs {info["size"]}' if len(seqs) != info['size'] else '')
    f.check(len(set(seqs)) == len(seqs), 'FASTA sequences are distinct',
            f'{len(seqs) - len(set(seqs))} duplicates' if len(set(seqs)) != len(seqs) else '')
    covered = set(r['Query Sequence'] for r in rows)
    uncovered = [s for s in seqs if s not in covered]
    f.check(not uncovered, 'every query has >=1 expected row',
            f'{len(uncovered)} with none, e.g. {uncovered[0]}' if uncovered else '')

    # 7 balance
    balance = Counter(r['Match Type'] for r in rows if r.get('Planted') == 'yes')
    ins, dele = balance['insertion_match'], balance['deletion_match']
    f.check(abs(ins - dele) <= 1, '50/50 insertion/deletion balance among planted rows',
            f'{ins} ins vs {dele} del')

    # 6 nesting
    if prev_queries is not None:
      f.check(prev_queries <= set(seqs), 'strictly nests the previous (smaller) subset',
              f'{len(prev_queries - set(seqs))} queries lost')
    prev_queries = set(seqs)

    # 6b subset expected == largest filtered to this subset
    if info is not largest:
      expect = sorted(
        tuple(r[c] for c in JOIN_COLUMNS)
        for q in seqs for r in largest_by_query.get(q, [])
      )
      got = sorted(tuple(r[c] for c in JOIN_COLUMNS) for r in rows)
      f.check(expect == got, "expected rows == the largest subset's rows filtered to these queries",
              f'{len(expect)} vs {len(got)}' if expect != got else f'{len(got):,} rows')

    # 10 manifest hashes
    f.check(sha256_file(BENCH / info['expected_tsv']) == info['expected_sha256'],
            'expected TSV matches its manifest sha256')
    f.check(sha256_file(BENCH / info['queries_fasta']) == info['queries_sha256'],
            'queries FASTA matches its manifest sha256')

  # 8 interior edits, via the engine's own definition
  print('\n--- interior-edit check (engine definition) ---')
  try:
    from pepmatch.matcher import _indel_placements
    # RANDOM sample, never a prefix. The expected TSV is sorted by Query Sequence, so
    # largest_rows[:CAP] is every alphabetically-early query -- it systematically skips
    # the rest of the composition space, so an interiority bug confined to (say) queries
    # starting W/Y/V would pass a gate that looks like it covered 200k rows. Seeded so
    # the check is reproducible, and the coverage fraction is PRINTED rather than implied.
    CAP = 200_000
    if len(largest_rows) <= CAP:
      sample = largest_rows
    else:
      sample = random.Random(20260825).sample(largest_rows, CAP)
    pct = 100.0 * len(sample) / max(len(largest_rows), 1)
    non_interior = [r for r in sample
                    if not _indel_placements(r['Query Sequence'], r['Matched Sequence'])]
    f.check(not non_interior,
            'every row has a legal non-terminal placement per _indel_placements',
            (f'{len(non_interior)} illegal, e.g. {non_interior[0]["Query Sequence"]} vs '
             f'{non_interior[0]["Matched Sequence"]}') if non_interior
            else f'{len(sample):,}/{len(largest_rows):,} rows sampled at random '
                 f'({pct:.1f}% coverage)')
  except ImportError as exc:
    f.check(False, 'import pepmatch.matcher._indel_placements', str(exc))

  # audit provenance
  audit = manifest.get('audit', {})
  f.check(audit.get('ok') and audit.get('sampled', 0) > 0,
          'exhaustive no-prefilter audit was run at generation time',
          audit.get('note', 'missing'))

  # 11/12 the real gate
  if not args.fast:
    target = args.run_tools_on or sizes[0]['size']
    info = next(i for i in files.values() if i['size'] == target)
    print(f'\n--- tool gates on the {target}-query subset ---')
    import pandas as pd
    import benchmarking as B

    expected_df = pd.read_csv(BENCH / info['expected_tsv'], sep='\t')
    # Construct the tools through the harness's own load_method and score with its own
    # recall(), so a pass here means the published table will report the same number.
    # Rebuilding either by hand would let the gate and the benchmark drift apart.
    dataset = {
      'query': info['queries_fasta'],
      'proteome': f'proteomes/{Path(args.proteome).name}',
      'lengths': sorted({len(s) for s in read_fasta(BENCH / info['queries_fasta'])[1]}),
      'mismatches': 1,
      'indels': indels,
    }

    for method_name, pretty in (('brute_force', 'Brute Force'), ('PEPMatch', 'PEPMatch')):
      try:
        tool = B.load_method(method_name, 'synthetic', dataset, {'threads': args.threads})
        for step in ('preprocess_proteome', 'preprocess_query'):
          try:
            getattr(tool, step)()
          except (TypeError, B.StepNotApplicable):
            pass          # method genuinely has no such phase (e.g. Brute Force)
        got = B.recall(tool.search(), expected_df)
        f.check(got == 100.0, f'{pretty} scores exactly 100.0', f'got {got}')
      except Exception as exc:
        f.check(False, f'{pretty} gate ran', f'{type(exc).__name__}: {exc}')

  print(f'\n{"=" * 72}')
  if f.items:
    print(f'  {len(f.items)} GATE(S) FAILED')
    for item in f.items:
      print(f'    - {item}')
    print('=' * 72)
    return 1
  print('  ALL GATES PASSED')
  print('=' * 72)
  return 0


if __name__ == '__main__':
  sys.exit(main())
