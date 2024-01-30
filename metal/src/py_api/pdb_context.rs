use pyo3::prelude::*;

use crate::core::pdb_context::PdbContext;

#[pyclass]
pub struct PyPdbContext {
    pdb_context : PdbContext
}

#[pymethods]
impl PyPdbContext {

    #[staticmethod]
    pub fn open_file(file: &str) -> PyResult<Self> {

        PdbContext::open_file(file)
            .map(|pdb_context| { PyPdbContext { pdb_context }})
            .map_err(|e| { e.into() })
            
    }
}