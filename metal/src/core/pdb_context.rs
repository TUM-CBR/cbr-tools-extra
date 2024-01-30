use pdbtbx::*;

use super::error::{CbrExtraError};

pub struct PdbContext {
    pdb : PDB
}

impl PdbContext {

    pub fn open_file(path: impl AsRef<str>) -> Result<Self, CbrExtraError> {

        let result =
            open_pdb(path, StrictnessLevel::Loose);

        match result {
            Ok((pdb, errors)) => Ok(PdbContext { pdb }),
            Err(errors) => Err(errors.into())
        }
    }
}