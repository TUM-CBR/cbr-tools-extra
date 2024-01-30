use pdbtbx::*;
use polars::error::PolarsError;

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

impl<'a, T: 'a> From<Vec<T>> for CbrExtraError where CbrExtraError : From<&'a T> {

    fn from(items: Vec<T>) -> CbrExtraError {
        CbrExtraError::Many(items.iter().map(|i| { i.into() }).collect())
    }
}