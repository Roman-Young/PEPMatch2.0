#!/usr/bin/env python3
"""Pre-flight check: do all tools agree on the 'Index start' convention?

`recall()` merges results against the expected file on an exact 4-column key that
includes Index start. If a tool reports 0-based coordinates where the expected file is
1-based, EVERY row misses and that tool scores 0% -- indistinguishable from a genuine
"found nothing" result. BLAST, PEPMatch and Brute Force were verified 1-based by
reading their code; DIAMOND and MMseqs2 apply no adjustment at all and simply trust
their raw sstart/tstart, which has never been checked against real output.

This builds a one-protein database, searches it with a peptide whose true position is
known, and asserts each tool reports the 1-based position.

Run it before submitting the cluster job:
    python3 benchmarking/verify_coordinates.py
Exit code 0 = all tools agree. Non-zero = do NOT trust that tool's recall.
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROTEOME = HERE / 'proteomes' / 'human.fasta'


def first_protein():
  """Return (header, sequence) of the first record in the human proteome."""
  header, chunks = None, []
  with open(PROTEOME) as f:
    for line in f:
      line = line.rstrip('\n')
      if line.startswith('>'):
        if header is not None:
          break
        header = line
      else:
        chunks.append(line)
  return header, ''.join(chunks)


def run(cmd, cwd):
  return subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True)


def main():
  header, seq = first_protein()
  # Take a peptide from the interior so its position is unambiguous and non-zero.
  offset0 = 100                      # 0-based
  expected_1based = offset0 + 1      # what every tool must report
  peptide = seq[offset0:offset0 + 12]
  if len(peptide) < 12:
    print('FAIL: first protein too short for the probe')
    return 1

  print(f'probe peptide : {peptide}')
  print(f'true position : {expected_1based} (1-based)\n')

  work = Path(tempfile.mkdtemp(prefix='coordcheck-'))
  try:
    mini = work / 'mini.fasta'
    mini.write_text(f'{header}\n{seq}\n')
    query = work / 'query.fasta'
    query.write_text(f'>probe\n{peptide}\n')

    results = {}

    # ---- BLAST (known-good control: if this is wrong, the harness is misconfigured)
    if shutil.which('makeblastdb') and shutil.which('blastp'):
      run(f'makeblastdb -in {mini} -dbtype prot', work)
      run(
        f'blastp -query {query} -db {mini} -task blastp-short -evalue 10000 '
        f'-comp_based_stats 0 -outfmt "6 sstart" -out blast.tsv', work
      )
      out = (work / 'blast.tsv').read_text().strip().splitlines()
      results['BLAST'] = int(out[0]) if out else None
    else:
      results['BLAST'] = 'binary missing'

    # ---- DIAMOND (the unverified one)
    if shutil.which('diamond'):
      run(f'diamond makedb --in {mini} -d {work}/mini', work)
      run(
        f'diamond blastp -d {work}/mini -q {query} -o dmnd.tsv -e 10000 -k 0 '
        f'--ultra-sensitive --masking 0 --comp-based-stats 0 --ignore-warnings '
        f'-f 6 sstart', work
      )
      out = (work / 'dmnd.tsv').read_text().strip().splitlines()
      results['DIAMOND'] = int(out[0]) if out else None
    else:
      results['DIAMOND'] = 'binary missing'

    # ---- MMseqs2 (the other unverified one)
    if shutil.which('mmseqs'):
      run(f'mmseqs createdb {mini} {work}/mmdb', work)
      run(
        f'mmseqs easy-search {query} {work}/mmdb mm.tsv {work}/tmp '
        f'-s 7.5 -e 10000 --max-seqs 100000 --mask 0 --comp-bias-corr 0 '
        f'--format-output "tstart"', work
      )
      out = (work / 'mm.tsv').read_text().strip().splitlines()
      results['MMseqs2'] = int(out[0]) if out else None
    else:
      results['MMseqs2'] = 'binary missing'

    print(f'{"Tool":10s} {"reported":>10s}   verdict')
    print('-' * 46)
    bad = 0
    for tool, got in results.items():
      if isinstance(got, str):
        verdict = f'SKIPPED ({got})'
      elif got is None:
        verdict = 'NO HITS -- cannot verify'
        bad += 1
      elif got == expected_1based:
        verdict = 'OK (1-based, matches expected)'
      elif got == offset0:
        verdict = 'BUG: 0-BASED -> recall would be a false 0%'
        bad += 1
      else:
        verdict = f'UNEXPECTED (wanted {expected_1based})'
        bad += 1
      print(f'{tool:10s} {str(got):>10s}   {verdict}')

    print()
    if bad:
      print(f'{bad} tool(s) failed the coordinate check -- DO NOT trust their recall.')
      return 1
    print('All available tools agree on the 1-based convention.')
    return 0
  finally:
    shutil.rmtree(work, ignore_errors=True)


if __name__ == '__main__':
  sys.exit(main())
