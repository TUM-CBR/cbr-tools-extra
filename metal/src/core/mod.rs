pub mod data;
pub mod error;
pub mod pdb_context;
pub mod pdb;

pub use pdb::*;
pub use pdb_context::*;

use polars::frame::DataFrame;

use self::error::CbrResult;

impl PdbContext {

    pub fn all_atoms(&self) -> CbrResult<DataFrame> {
        pdb::from_pdb(
            pdb::PdbReadArgs::Default,
            self.pdb()
        )
    }
}