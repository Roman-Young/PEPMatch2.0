#!/usr/bin/env python3

import os
import glob
import pandas as pd

from Bio import SeqIO


directory = os.path.dirname(os.path.abspath(__file__))


def parse_fasta(file):
  return SeqIO.parse(file, 'fasta')


class BLAST(object):
  def __init__(self, query, proteome, max_mismatches, method_parameters, indels=0):
    if max_mismatches == -1:
      max_mismatches = 7

    self.query = query
    self.proteome = proteome
    self.max_mismatches = max_mismatches
    self.indels = indels
    self.threads = int(method_parameters.get('threads', 1))

    bin_directory = method_parameters['bin_directory']
    self.makeblastdb_path = os.path.join(bin_directory, 'makeblastdb')
    self.blastp_path = os.path.join(bin_directory, 'blastp')

  def __str__(self):
    return 'BLAST'

  def preprocess(self):
    os.system(f"{self.makeblastdb_path} -in {self.proteome} -dbtype prot")

  def blast_search(self):
    peptides = parse_fasta(self.query)
    proteins = parse_fasta(self.proteome)

    peptide_dict = {}
    for peptide in peptides:
      peptide_dict[str(peptide.id)] = str(peptide.seq)

    protein_dict = {}
    for protein in proteins:
      try:
        protein_dict[str(protein.id).split('|')[1]] = str(protein.seq)
      except IndexError:
        protein_dict[str(protein.id)] = str(protein.seq)

    # NcbiblastpCommandline was removed in Biopython 1.80, so call blastp directly.
    # Tuned for maximum recall on short peptides: blastp-short (the correct task for
    # short queries -- default blastp is built for long proteins), uncapped targets,
    # no composition-based score adjustment.
    evalue = 100 if self.max_mismatches == 0 else 10000
    os.system(
      f"{self.blastp_path} -query {self.query} -db {self.proteome} "
      f"-task blastp-short -evalue {evalue} -max_target_seqs 100000 "
      f"-comp_based_stats 0 -num_threads {self.threads} -outfmt 10 -out output.csv"
    )

    df = pd.read_csv(
      'output.csv',
      names=[
        'Peptide Sequence', 'Protein ID', 'Sequence Identity',
        'Length', 'Mismatches', 'Gap Openings', 'Query start',
        'Query end', 'Index start', 'Index end', 'e value', 'bit score'
      ]
    )
    df['Peptide Sequence'] = df['Peptide Sequence'].apply(str)
    df = df.replace({'Peptide Sequence': peptide_dict})

    all_matches = []
    for i, row in df.iterrows():
      try:
        # for UniProt IDs - betacoronaviruses have different NCBI IDs
        row['Protein ID'] = row['Protein ID'].split('|')[1]
      except IndexError:
        pass

      peptide_sequence = row['Peptide Sequence']
      protein_id = row['Protein ID']
      index_start = int(row['Index start']) - 1  # BLAST is 1-indexed

      if self.indels > 0:
        # An indel match is not the query's length; use BLAST's own subject-end
        # coordinate so the reported match is the true aligned substring.
        matched = protein_dict[protein_id][index_start:int(row['Index end'])]
      else:
        matched = protein_dict[protein_id][index_start:index_start + len(peptide_sequence)]

      all_matches.append((peptide_sequence, matched, protein_id, index_start + 1))

    return all_matches


class Benchmarker(BLAST):
  def __init__(
    self, benchmark: str, query: str, proteome: str, lengths: list, max_mismatches: int,
    method_parameters: dict, indels: int = 0
  ):
    self.benchmark = benchmark
    self.query = query
    self.proteome = proteome
    self.lengths = lengths
    self.max_mismatches = max_mismatches
    self.method_parameters = method_parameters

    super().__init__(query, proteome, max_mismatches, method_parameters, indels)

  def __str__(self):
    return 'BLAST'

  def preprocess_proteome(self):
    return self.preprocess()

  def preprocess_query(self):
    raise TypeError(self.__str__() + ' does not preprocess queries.\n')

  def search(self):
    matches = self.blast_search()

    all_matches = []
    for match in matches:
      match = list(match)
      try:  # try taking the UniProt ID - else do nothing
        match[2] = match[2].split('|')[1]
      except IndexError:
        pass
      all_matches.append([str(i) for i in match])

    columns = ['Query Sequence', 'Matched Sequence', 'Protein ID', 'Index start']

    # clean up BLAST db files + output, tolerating whatever extension set this
    # BLAST version wrote.
    db_dir = os.path.dirname(str(self.proteome))
    for pattern in ['*.pdb', '*.phr', '*.pin', '*.psq', '*.ptf', '*.pot', '*.pto', '*.pjs', '*.pdb']:
      for f in glob.glob(os.path.join(db_dir, pattern)):
        try:
          os.remove(f)
        except OSError:
          pass
    for f in ['output.csv']:
      try:
        os.remove(f)
      except OSError:
        pass

    return pd.DataFrame(all_matches, columns=columns)
