use pdbtbx::*;
use polars::error::{PolarsError, PolarsWarning};

use crate::py_api::error;

pub type CbrResult<T> = Result<T, CbrExtraError>;

#[derive(Clone)]
pub enum CbrExtraError {
    Unknown(),
    Many(Vec<Self>)
}

impl From<&PDBError> for CbrExtraError {

    fn from(error: &PDBError) -> CbrExtraError {
        CbrExtraError::Unknown()
    }
}

impl From<PDBError> for CbrExtraError {

    fn from(error: PDBError) -> CbrExtraError {
        error.into()
    }
}

impl<'a, T: 'a> From<&'a Vec<T>> for CbrExtraError where CbrExtraError : From<&'a T> {

    fn from(items: &'a Vec<T>) -> CbrExtraError {
        CbrExtraError::Many(
            items
                .iter()
                .map(|i| { i.into() })
                .collect()
            )
    }
}

impl From<&PolarsError> for CbrExtraError {

    fn from(error: &PolarsError) -> CbrExtraError {
        CbrExtraError::Unknown()
    }
}

impl From<PolarsError> for CbrExtraError {

    fn from(error: PolarsError) -> CbrExtraError {
        error.into()
    }
}

pub trait AsResult<T> {

    fn as_result(self) -> CbrResult<T>;
}

impl<TValue, TError> AsResult<TValue> for Result<TValue, TError>
    where CbrExtraError : From<TError> {

    fn as_result(self) -> CbrResult<TValue> {
        self.map_err(|e| { e.into() })
    }
}