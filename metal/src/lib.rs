use pyo3::prelude::*;

pub mod algorithms;
pub mod core;
pub mod py_api;

use py_api::pdb_context::PyPdbContext;

/// A Python module implemented in Rust.
#[pymodule]
fn cbrextra_metal(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyPdbContext>().unwrap();
    Ok(())
}
