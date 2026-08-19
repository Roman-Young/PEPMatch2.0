#!/usr/bin/env python3

import os
import glob
import pandas as pd
from Bio import SeqIO


directory = os.path.dirname(os.path.abspath(__file__))


def parse_fasta(file):
  return SeqIO.parse(file, 'fasta')


class DIAMOND(object):
  def __init__(self, query, proteome, max_mismatches, method_parameters):
    if max_mismatches == -1:
      max_mismatches = 7

    self.query = query
    self.proteome = proteome
    self.proteome_name = str(proteome).replace('.fasta', '')

    self.max_mismatches = max_mismatches
    self.threads = int(method_parameters.get('threads', 1))

    bin_directory = method_parameters['bin_directory']
    self.bin_file = os.path.join(bin_directory, 'diamond')

  def __str__(self):
    return 'DIAMOND'

  def preprocesss(self):
    os.system(
      f"{self.bin_file} makedb --in {self.proteome} -d {self.proteome_name} "
      f"--threads {self.threads}"
    )

  def diamond_search(self):
    # Tuned for maximum recall: report every alignment (-k 0), most sensitive mode,
    # masking off, no composition-based score adjustment.
    # --ignore-warnings: low-complexity peptides (e.g. poly-A/poly-C) are all valid DNA
    # letters, so DIAMOND misdetects the input as nucleotide and aborts without it.
    os.system(
      f"{self.bin_file} blastp -d {self.proteome_name} -q {self.query} -o matches.m8 "
      f"-e 10000 -k 0 --ultra-sensitive --masking 0 --comp-based-stats 0 --ignore-warnings "
      f"--threads {self.threads} -f 6 "
       "full_qseq sseq sseqid mismatch sstart"
    )

    all_matches = []
    with open('matches.m8', 'r') as file:
      lines = file.readlines()
      for line in lines:
        match = []
        result = line.split('\t')
        for i in range(len(result)):
          if i == 1:
            # aligned subject sequence: carries '-' gap chars on indel alignments;
            # stripping them yields the true contiguous match (a no-op when ungapped).
            match.append(result[i].replace('-', ''))
          elif i == 3:
            continue  # skip the max_mismatches column
          elif i == 4:
            match.append(int(result[i].replace('\n', '')))
          else:
            match.append(result[i])

        all_matches.append(match)

    return all_matches


class Benchmarker(DIAMOND):
  def __init__(
    self, benchmark: str, query: str, proteome: str, lengths: list, max_mismatches: int,
    method_parameters: dict
  ):
    self.benchmark = benchmark
    self.query = query
    self.proteome = proteome
    self.lengths = lengths
    self.max_mismatches = max_mismatches
    self.method_parameters = method_parameters

    super().__init__(query, proteome, max_mismatches, method_parameters)

  def __str__(self):
    return 'DIAMOND'

  def preprocess_proteome(self):
    return self.preprocesss()

  def preprocess_query(self):
    raise TypeError(self.__str__() + ' does not preprocess queries.\n')

  def search(self):
    matches = self.diamond_search()

    all_matches = []
    for match in matches:
      match = list(match)
      try:  # get the UniProt ID or do nothing
        match[2] = match[2].split('|')[1]
      except IndexError:
        pass
      all_matches.append([str(i) for i in match])

    columns = ['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start']

    for f in ['matches.m8'] + glob.glob(f"{self.proteome_name}.dmnd"):
      try:
        os.remove(f)
      except OSError:
        pass

    return pd.DataFrame(all_matches, columns=columns)
