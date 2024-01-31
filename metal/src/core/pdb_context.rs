use pdbtbx::*;

use super::data::*;
use super::error::*;

pub struct PdbContext {
    pdb : PDB
}

impl PdbContext {

    /// Open a file containing a PDB structure and load it into memory.
    pub fn open_file(path: impl AsRef<str>) -> Result<Self, CbrExtraError> {

        let result =
            open_pdb(path, StrictnessLevel::Loose);

        match result {
            Ok((pdb, errors)) => Ok(PdbContext { pdb }),
            Err(errors) => Err((&errors).into())
        }
    }

    pub fn pdb(&self) -> &PDB {
        &self.pdb
    }

    pub fn bounding_box(&self) -> (Point3d, Point3d) {
        self.pdb.bounding_box()
    }
}