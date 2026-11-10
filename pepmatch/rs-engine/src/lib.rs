//! PyO3 boundary: thin wrappers over `preprocess` (build a .pepidx) and `matching` (query
//! one). Results are returned as tuples of parallel Vecs so matcher.py can feed each one
//! straight into a polars Series instead of building per-row objects.

use pyo3::prelude::*;

mod preprocess;
// `match` is a Rust keyword, so the file cannot be a module name; import it under an alias.
#[path = "match.rs"]
mod matching;

#[pyfunction]
fn rs_version() -> &'static str {
    "0.1.0"
}

#[pyfunction]
fn rs_preprocess(fasta_path: &str, k: usize, output_path: &str) {
    preprocess::run(fasta_path, k, output_path);
}

#[pyfunction]
fn rs_match(pepidx_path: &str, peptides: Vec<(String, String)>, k: usize, max_mismatches: usize) -> matching::Columns {
    matching::run(pepidx_path, peptides, k, max_mismatches)
}

#[pyfunction]
fn rs_discontinuous(pepidx_path: &str, epitopes: Vec<(String, Vec<(char, usize)>)>, max_mismatches: usize) -> matching::Columns {
    matching::run_discontinuous(pepidx_path, epitopes, max_mismatches)
}

/// Metadata for `protein_nums` (1-based FASTA ordinals in this index), or for every protein
/// when it is omitted. A number outside the index raises `ValueError`.
#[pyfunction]
#[pyo3(signature = (pepidx_path, protein_nums=None))]
fn rs_metadata(pepidx_path: &str, protein_nums: Option<Vec<u32>>) -> PyResult<matching::MetaColumns> {
    matching::run_metadata(pepidx_path, protein_nums).map_err(pyo3::exceptions::PyValueError::new_err)
}

/// Proteins with `min_len <= length < max_len` in this index, as (protein numbers, sequences).
#[pyfunction]
fn rs_short_proteins(pepidx_path: &str, min_len: usize, max_len: usize) -> (Vec<u32>, Vec<String>) {
    matching::run_short_proteins(pepidx_path, min_len, max_len)
}

/// (length, count) for proteins with `min_len <= length < max_len` in this index, ascending.
#[pyfunction]
fn rs_short_lengths(pepidx_path: &str, min_len: usize, max_len: usize) -> Vec<(usize, u64)> {
    matching::run_short_lengths(pepidx_path, min_len, max_len)
}

/// Searches the proteins of `proteome_path` with `min_len <= length < max_len` against the
/// queries index and returns only the hits (see matching::run_short_target_match).
#[pyfunction]
fn rs_short_target_match(
    proteome_path: &str, queries_path: &str, min_len: usize, max_len: usize, k: usize, edits: usize, indels: bool,
) -> matching::ShortHitColumns {
    matching::run_short_target_match(proteome_path, queries_path, min_len, max_len, k, edits, indels)
}

#[pyfunction]
fn rs_match_counts(pepidx_path: &str, peptides: Vec<(String, String)>, k: usize, max_mismatches: usize) -> matching::CountColumns {
    matching::run_counts(pepidx_path, peptides, k, max_mismatches)
}

#[pyfunction]
fn rs_indel_match(pepidx_path: &str, peptides: Vec<(String, String)>, indels_allowed: usize) -> matching::Columns {
    matching::run_indel(pepidx_path, peptides, indels_allowed)
}

#[pymodule]
fn _rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(rs_version, m)?)?;
    m.add_function(wrap_pyfunction!(rs_preprocess, m)?)?;
    m.add_function(wrap_pyfunction!(rs_match, m)?)?;
    m.add_function(wrap_pyfunction!(rs_discontinuous, m)?)?;
    m.add_function(wrap_pyfunction!(rs_metadata, m)?)?;
    m.add_function(wrap_pyfunction!(rs_short_proteins, m)?)?;
    m.add_function(wrap_pyfunction!(rs_short_lengths, m)?)?;
    m.add_function(wrap_pyfunction!(rs_short_target_match, m)?)?;
    m.add_function(wrap_pyfunction!(rs_match_counts, m)?)?;
    m.add_function(wrap_pyfunction!(rs_indel_match, m)?)?;
    Ok(())
}
