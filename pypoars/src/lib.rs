use pyo3::prelude::*;
use bio::alignment::{pairwise::Scoring, poa::Aligner};


/// Formats the sum of two numbers as string.
#[pyfunction]
fn poa_consensus(seqs: Vec<String>) -> PyResult<String> {
    let seqs_bytes = seqs.iter().map(|seq| seq.bytes().collect::<Vec<u8>>()).collect::<Vec<Vec<u8>>>();
    let scoring = Scoring::new(-12, -6, |a: u8, b: u8| if a == b { 3 } else { -4 });
    let mut aligner = Aligner::new(scoring, &seqs_bytes[0]);
    for seq in seqs_bytes.iter().skip(1) {
        aligner.global(seq).add_to_graph();
    }
    let consensus = aligner.consensus();
    Ok(std::str::from_utf8(&consensus).unwrap().to_string())
}

/// A Python module implemented in Rust.
#[pymodule]
fn pypoars(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(poa_consensus, m)?)?;
    Ok(())
}
