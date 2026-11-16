import argparse
import re
from .preprocessor import Preprocessor
from .matcher import Matcher, check_search_options, index_ks


def _search_arguments(parser):
  """The options that decide which indexes a search reads. pepmatch-match and
  pepmatch-preprocess share them, so the same command line prebuilds what it later searches."""
  parser.add_argument('-p', '--proteome_file', '--proteome', dest='proteome_file', required=True)
  parser.add_argument('-m', '--max_mismatches', type=int, default=0)
  parser.add_argument('-i', '--max_indels', type=int, default=0)
  parser.add_argument('-b', '--best_match', action='store_true', default=False)
  parser.add_argument('-P', '--preprocessed_files_path', default='.')
  parser.add_argument('-n', '--proteome_name', default='')


def _lengths(text):
  """Peptide lengths from e.g. '9', '8-15' or '8-11,15'."""
  lengths = set()
  for part in text.split(','):
    span = re.fullmatch(r'\s*(\d+)\s*(?:-\s*(\d+)\s*)?', part)
    low, high = (int(span[1]), int(span[2] or span[1])) if span else (0, 0)
    if not 1 <= low <= high:
      lengths = set()
      break
    lengths.update(range(low, high + 1))
  if not lengths:
    raise argparse.ArgumentTypeError(f"{text!r} is not a list of peptide lengths such as '9', '8-15' or '8-11,15'")
  return sorted(lengths)


def run_preprocessor():
  """Build the indexes a search will read, without searching: for the peptides in a query
  file, or for a range of peptide lengths when the queries are not known yet."""
  parser = argparse.ArgumentParser()
  queries = parser.add_mutually_exclusive_group(required=True)
  queries.add_argument('-q', '--query')
  queries.add_argument('-l', '--lengths', type=_lengths)
  _search_arguments(parser)
  args = parser.parse_args()

  if args.query:
    Matcher(
      query=args.query,
      proteome_file=args.proteome_file,
      max_mismatches=args.max_mismatches,
      max_indels=args.max_indels,
      preprocessed_files_path=args.preprocessed_files_path,
      best_match=args.best_match,
      proteome_name=args.proteome_name,
    ).preprocess()
  else:
    check_search_options(args.max_mismatches, args.max_indels, args.best_match)
    Preprocessor(
      proteome=args.proteome_file,
      proteome_name=args.proteome_name,
      preprocessed_files_path=args.preprocessed_files_path,
    ).preprocess_all(index_ks(args.lengths, max(args.max_mismatches, args.max_indels), args.best_match))


def run_matcher():
  parser = argparse.ArgumentParser()
  parser.add_argument('-q', '--query', required=True)
  _search_arguments(parser)
  parser.add_argument('-f', '--output_format', default='csv')
  parser.add_argument('-o', '--output_name', default='')
  parser.add_argument('-v', '--sequence_version', action='store_false', default=True)
  args = parser.parse_args()

  Matcher(
    query=args.query,
    proteome_file=args.proteome_file,
    max_mismatches=args.max_mismatches,
    max_indels=args.max_indels,
    preprocessed_files_path=args.preprocessed_files_path,
    best_match=args.best_match,
    output_format=args.output_format,
    output_name=args.output_name,
    sequence_version=args.sequence_version,
    proteome_name=args.proteome_name,
  ).match()
