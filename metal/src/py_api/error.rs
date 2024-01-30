use pyo3::{PyErr, PyResult};
use pyo3::exceptions::PyException;

use crate::core::error::{CbrExtraError};

impl From<CbrExtraError> for PyErr {

    fn from(error: CbrExtraError) -> Self {
        PyException::new_err("")
    }
}