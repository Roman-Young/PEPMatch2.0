//! Builds the `.pepidx` that all searching reads. Layout is documented at `write_pepidx`;
//! match.rs mirrors it byte for byte, so a change here needs the same change there.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write, BufWriter};

use rayon::prelude::*;

// Every k-mer occurrence is stored as ONE u64: protein_number * MULTIPLIER + offset, so the
// multiplier is also the hard cap on protein length.
const PROTEIN_INDEX_MULTIPLIER: u64 = 100_000_000;

// The build holds one 24-byte record per k-mer occurrence it is sorting (0.3 GB for the human
// proteome at k=9). Up to ONE_PASS_MAX_OCCURRENCES everything is sorted at once; larger inputs
// are sorted and written one slice of at most SLICE_OCCURRENCES at a time (1.2 GB), which stays
// below the 8+ bytes per occurrence the earlier chunk-and-merge builder held at those sizes.
const ONE_PASS_MAX_OCCURRENCES: usize = 200_000_000;
const SLICE_OCCURRENCES: usize = 50_000_000;

struct Protein {
    header: String,
    sequence: String,
}

struct Metadata {
    protein_id: String,
    protein_name: String,
    species: String,
    taxon_id: String,
    gene: String,
    pe_level: String,
    sequence_version: String,
    gene_priority: String,
    swissprot: String,
}

/// Every `>` line starts a record, even one with an empty header, so records never merge.
/// Lines before the first header belong to no record and are left out (the Python side warns).
/// Sequences are uppercased because queries are.
fn parse_fasta(path: &str) -> Vec<Protein> {
    let file = File::open(path).expect("Cannot open FASTA file");
    let reader = BufReader::new(file);
    let mut proteins: Vec<Protein> = Vec::new();
    let mut current: Option<Protein> = None;

    for line in reader.lines() {
        let line = line.expect("Cannot read line");
        if let Some(header) = line.strip_prefix('>') {
            proteins.extend(current.take());
            current = Some(Protein { header: header.to_string(), sequence: String::new() });
        } else if let Some(protein) = current.as_mut() {
            protein.sequence.push_str(&line.trim().to_ascii_uppercase());
        }
    }
    proteins.extend(current);
    proteins
}

fn extract_between<'a>(header: &'a str, start_tag: &str, end_tag: &str) -> &'a str {
    if let Some(start) = header.find(start_tag) {
        let value_start = start + start_tag.len();
        let rest = &header[value_start..];
        let end = rest.find(end_tag).unwrap_or(rest.len());
        rest[..end].trim()
    } else {
        ""
    }
}

fn extract_field<'a>(header: &'a str, prefix: &str) -> &'a str {
    if let Some(start) = header.find(prefix) {
        let value_start = start + prefix.len();
        let rest = &header[value_start..];
        let end = rest.find(|c: char| c.is_whitespace()).unwrap_or(rest.len());
        rest[..end].trim()
    } else {
        ""
    }
}

fn extract_protein_id(header: &str) -> String {
    if let Some(first_pipe) = header.find('|') {
        let rest = &header[first_pipe + 1..];
        if let Some(second_pipe) = rest.find('|') {
            return rest[..second_pipe].to_string();
        }
    }
    header.split_whitespace().next().unwrap_or("").to_string()
}

fn extract_protein_name(header: &str) -> String {
    if let Some(space_pos) = header.find(' ') {
        let after_id = &header[space_pos + 1..];
        if let Some(os_pos) = after_id.find(" OS=") {
            return after_id[..os_pos].trim().to_string();
        }
        return after_id.trim().to_string();
    }
    String::new()
}

fn is_swissprot(header: &str) -> bool {
    header.starts_with("sp|")
}

/// UniProt marks isoforms with a dash suffix (`P12345-2`).
fn has_isoform_dash(protein_id: &str) -> bool {
    protein_id.contains('-')
}

fn base_accession(protein_id: &str) -> &str {
    if let Some(dash_pos) = protein_id.find('-') {
        &protein_id[..dash_pos]
    } else {
        protein_id
    }
}

fn extract_all_metadata(proteins: &[Protein]) -> Vec<Metadata> {
    let mut canonical_pe: HashMap<String, String> = HashMap::new();

    // Pass 1: isoform records often omit PE (protein existence) level and inherit their
    // canonical entry's, which may appear later in the file.
    for protein in proteins {
        let pid = extract_protein_id(&protein.header);
        if !has_isoform_dash(&pid) {
            let pe = extract_field(&protein.header, "PE=");
            if !pe.is_empty() {
                canonical_pe.insert(pid.clone(), pe.to_string());
            }
        }
    }

    let mut seen_genes: HashSet<String> = HashSet::new();
    let mut metadata_list: Vec<Metadata> = Vec::with_capacity(proteins.len());

    // `seen_genes` makes gene_priority "first entry wins": one preferred protein per gene.
    for protein in proteins {
        let h = &protein.header;
        let protein_id = extract_protein_id(h);
        let protein_name = extract_protein_name(h);
        let species = extract_between(h, "OS=", " OX=").to_string();
        let taxon_id = extract_field(h, "OX=").to_string();
        let gene = extract_field(h, "GN=").to_string();
        let sv = extract_field(h, "SV=");
        let sequence_version = if sv.is_empty() { "1".to_string() } else { sv.to_string() };
        let swissprot = if is_swissprot(h) { "1".to_string() } else { "0".to_string() };

        let mut pe_level = extract_field(h, "PE=").to_string();
        if has_isoform_dash(&protein_id) && (pe_level.is_empty() || pe_level == "0") {
            let base = base_accession(&protein_id);
            if let Some(canonical) = canonical_pe.get(base) {
                pe_level = canonical.clone();
            }
        }
        if pe_level.is_empty() {
            pe_level = "0".to_string();
        }

        let gene_priority = if is_swissprot(h)
            && !has_isoform_dash(&protein_id)
            && !gene.is_empty()
            && !seen_genes.contains(&gene)
        {
            seen_genes.insert(gene.clone());
            "1".to_string()
        } else {
            "0".to_string()
        };

        metadata_list.push(Metadata {
            protein_id,
            protein_name,
            species,
            taxon_id,
            gene,
            pe_level,
            sequence_version,
            gene_priority,
            swissprot,
        });
    }

    metadata_list
}

fn write_length_prefixed(buf: &mut Vec<u8>, s: &str) {
    let bytes = s.as_bytes();
    buf.extend_from_slice(&(bytes.len() as u16).to_le_bytes());
    buf.extend_from_slice(bytes);
}


/// One k-mer occurrence. The k-mer is packed big-endian into (hi, lo), so the derived order is
/// the k-mers' byte order and then position order: the order the index stores them in.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Occurrence {
    hi: u64,
    lo: u64,
    position: u64,
}

/// End of the run of occurrences that share the k-mer at `start`.
fn kmer_run_end(occurrences: &[Occurrence], start: usize) -> usize {
    let first = occurrences[start];
    let mut end = start + 1;
    while end < occurrences.len() && occurrences[end].hi == first.hi && occurrences[end].lo == first.lo {
        end += 1;
    }
    end
}

/// The first two residues of every k-mer as a u16 (k is at least 2). K-mers sort by these
/// first, so a contiguous range of prefixes is a contiguous stretch of the index.
fn for_each_prefix(proteins: &[Protein], k: usize, mut visit: impl FnMut(usize)) {
    for protein in proteins {
        let seq = protein.sequence.as_bytes();
        for j in 0..(seq.len() + 1).saturating_sub(k) {
            visit((seq[j] as usize) << 8 | seq[j + 1] as usize);
        }
    }
}

/// Splits the k-mers into contiguous prefix ranges of at most `slice` occurrences each (one
/// range when the total fits in `one_pass_max`), returned as (first prefix, end prefix,
/// occurrences). A single prefix with more occurrences than `slice` gets a range of its own.
fn slice_ranges(proteins: &[Protein], k: usize, total: usize, one_pass_max: usize, slice: usize) -> Vec<(usize, usize, usize)> {
    if total <= one_pass_max {
        return vec![(0, 1 << 16, total)];
    }
    let mut counts = vec![0usize; 1 << 16];
    for_each_prefix(proteins, k, |prefix| counts[prefix] += 1);
    let mut ranges = Vec::new();
    let (mut start, mut held) = (0, 0);
    for (prefix, &count) in counts.iter().enumerate() {
        if held > 0 && held + count > slice {
            ranges.push((start, prefix, held));
            start = prefix;
            held = 0;
        }
        held += count;
    }
    ranges.push((start, 1 << 16, held));
    ranges
}

/// Every occurrence of a k-mer whose prefix lies in `range`, sorted on all cores.
fn sorted_occurrences(proteins: &[Protein], k: usize, range: (usize, usize, usize)) -> Vec<Occurrence> {
    let mask: u128 = if k == 16 { u128::MAX } else { (1u128 << (8 * k)) - 1 };
    let shift = 8 * (k - 2);
    let mut occurrences = Vec::with_capacity(range.2);
    for (protein_count, protein) in proteins.iter().enumerate() {
        let protein_number = (protein_count + 1) as u64;
        // Rolling pack: each residue shifts in and the mask drops the one that left the k-mer.
        let mut kmer_val: u128 = 0;
        for (j, &b) in protein.sequence.as_bytes().iter().enumerate() {
            kmer_val = ((kmer_val << 8) | b as u128) & mask;
            if j + 1 >= k {
                let prefix = (kmer_val >> shift) as usize;
                if range.0 <= prefix && prefix < range.1 {
                    occurrences.push(Occurrence {
                        hi: (kmer_val >> 64) as u64,
                        lo: kmer_val as u64,
                        position: protein_number * PROTEIN_INDEX_MULTIPLIER + (j + 1 - k) as u64,
                    });
                }
            }
        }
    }
    occurrences.par_sort_unstable();
    occurrences
}

/// Layout, in order: 8-byte magic, k (u8), protein count (u32), sequence length (u64),
/// unique k-mer count (u64), total positions (u64), metadata length (u64); then the
/// concatenated sequences, protein offsets, metadata offsets, metadata blob, k-mer offsets,
/// k-mer data. Every section start is derivable from the fixed-size header, which is what
/// lets the query side mmap the file and binary-search it without parsing.
///
/// Each k-mer record is `[k bytes of k-mer][u32 count][count * u64 position]`, k-mers in byte
/// order and positions ascending. The k-mer offsets come before the records, so the slices
/// are sorted once to measure the records and again to write them; with one slice (any input
/// up to `one_pass_max`) the sorted occurrences are kept and sorted only once. The index is
/// written to a temporary name and renamed into place when complete, so an interrupted build
/// never leaves a partial .pepidx behind.
fn write_pepidx(path: &str, k: usize, proteins: &[Protein], metadata_list: &[Metadata], one_pass_max: usize, slice: usize) {
    let mut concat_seq = Vec::new();
    let mut protein_offsets: Vec<u64> = Vec::with_capacity(proteins.len());
    for protein in proteins {
        protein_offsets.push(concat_seq.len() as u64);
        concat_seq.extend_from_slice(protein.sequence.as_bytes());
    }

    let mut metadata_buf: Vec<u8> = Vec::new();
    let mut metadata_offsets: Vec<u64> = Vec::with_capacity(proteins.len());
    for meta in metadata_list {
        metadata_offsets.push(metadata_buf.len() as u64);
        write_length_prefixed(&mut metadata_buf, &meta.protein_id);
        write_length_prefixed(&mut metadata_buf, &meta.protein_name);
        write_length_prefixed(&mut metadata_buf, &meta.species);
        write_length_prefixed(&mut metadata_buf, &meta.taxon_id);
        write_length_prefixed(&mut metadata_buf, &meta.gene);
        write_length_prefixed(&mut metadata_buf, &meta.pe_level);
        write_length_prefixed(&mut metadata_buf, &meta.sequence_version);
        write_length_prefixed(&mut metadata_buf, &meta.gene_priority);
        write_length_prefixed(&mut metadata_buf, &meta.swissprot);
    }

    let total: usize = proteins.iter().map(|p| (p.sequence.len() + 1).saturating_sub(k)).sum();
    let ranges = slice_ranges(proteins, k, total, one_pass_max, slice);

    // Sorted, each k-mer is one contiguous run, so its record's offset is the bytes before it.
    let mut kept: Option<Vec<Occurrence>> = None;
    let mut kmer_offsets: Vec<u64> = Vec::new();
    let mut kmer_data_len: u64 = 0;
    for &range in &ranges {
        let occurrences = sorted_occurrences(proteins, k, range);
        let mut start = 0;
        while start < occurrences.len() {
            let end = kmer_run_end(&occurrences, start);
            kmer_offsets.push(kmer_data_len);
            kmer_data_len += (k + 4 + 8 * (end - start)) as u64;
            start = end;
        }
        if ranges.len() == 1 {
            kept = Some(occurrences);
        }
    }

    let tmp_path = format!("{}.tmp.{}", path, std::process::id());
    let file = File::create(&tmp_path).expect("Cannot create output file");
    let mut writer = BufWriter::new(file);
    writer.write_all(b"PEPIDX05").unwrap();
    writer.write_all(&(k as u8).to_le_bytes()).unwrap();
    writer.write_all(&(proteins.len() as u32).to_le_bytes()).unwrap();
    writer.write_all(&(concat_seq.len() as u64).to_le_bytes()).unwrap();
    writer.write_all(&(kmer_offsets.len() as u64).to_le_bytes()).unwrap();
    writer.write_all(&(total as u64).to_le_bytes()).unwrap();
    writer.write_all(&(metadata_buf.len() as u64).to_le_bytes()).unwrap();
    writer.write_all(&concat_seq).unwrap();
    for &offset in &protein_offsets {
        writer.write_all(&offset.to_le_bytes()).unwrap();
    }
    for &offset in &metadata_offsets {
        writer.write_all(&offset.to_le_bytes()).unwrap();
    }
    writer.write_all(&metadata_buf).unwrap();
    for &offset in &kmer_offsets {
        writer.write_all(&offset.to_le_bytes()).unwrap();
    }
    for &range in &ranges {
        let occurrences = kept.take().unwrap_or_else(|| sorted_occurrences(proteins, k, range));
        let mut start = 0;
        while start < occurrences.len() {
            let end = kmer_run_end(&occurrences, start);
            let kmer = ((occurrences[start].hi as u128) << 64) | occurrences[start].lo as u128;
            writer.write_all(&kmer.to_be_bytes()[16 - k..]).unwrap();
            writer.write_all(&((end - start) as u32).to_le_bytes()).unwrap();
            for occurrence in &occurrences[start..end] {
                writer.write_all(&occurrence.position.to_le_bytes()).unwrap();
            }
            start = end;
        }
    }
    writer.flush().unwrap();
    drop(writer);
    std::fs::rename(&tmp_path, path).expect("Cannot move the finished index into place");
}

/// `max_slice_occurrences` is None in normal use; tests set it to build small inputs in several
/// slices of at most that many occurrences.
pub(crate) fn run(fasta_path: &str, k: usize, output_path: &str, max_slice_occurrences: Option<usize>) {
    let proteins = parse_fasta(fasta_path);
    let metadata_list = extract_all_metadata(&proteins);
    let (one_pass_max, slice) = match max_slice_occurrences {
        Some(bound) => (bound.max(1), bound.max(1)),
        None => (ONE_PASS_MAX_OCCURRENCES, SLICE_OCCURRENCES),
    };
    write_pepidx(output_path, k, &proteins, &metadata_list, one_pass_max, slice);
}
