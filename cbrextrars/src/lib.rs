use pyo3::prelude::*;

pub mod histon;
pub mod protdeflator;

use protdeflator::screen::sum;

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok(sum(a,b).to_string())
}

/// A Python module implemented in Rust.
#[pymodule]
fn cbrextrars(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    Ok(())
}
