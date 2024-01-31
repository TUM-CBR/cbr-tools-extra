use pdbtbx::{ContainsAtomConformer, ContainsAtomConformerResidue, PDB};
use polars::prelude::*;

use crate::core::error::*;

pub enum PdbReadArgs {
    Default,
}

pub enum PdbColumn {
    XCoord,
    YCoord,
    ZCoord
}

impl From<PdbColumn> for &str {

    fn from(column: PdbColumn) -> &'static str {

        match column {
            PdbColumn::XCoord => "x_coord",
            PdbColumn::YCoord => "y_coord",
            PdbColumn::ZCoord => "z_coord"
        }
    }
}

/// Construct a `DataFrame` from the atoms in the given pdb object.
/// An instance of `PdbReadArgs` can be used to control what atoms
/// and fields to include in the final result. The fields in the
/// resulting `DataFrame` will be named using the string representations
/// of `PdbColumn`.
pub fn from_pdb(
    _args: PdbReadArgs,
    pdb: &PDB
) -> CbrResult<DataFrame> {

    //let atoms = Series::new("atoms", pdb.atoms_with_hierarchy_mut().into());
    let atom_count = pdb.atom_count();
    let mut x_coords : Vec<f64> = Vec::with_capacity(atom_count);
    let mut y_coords : Vec<f64> = Vec::with_capacity(atom_count);
    let mut z_coords : Vec<f64> = Vec::with_capacity(atom_count);

    for entry in pdb.atoms_with_hierarchy() {
        let atom = entry.atom();
        x_coords.push(atom.x());
        y_coords.push(atom.y());
        z_coords.push(atom.z());
    }

    DataFrame::new(vec![
        Series::new(
            PdbColumn::XCoord.into(),
            x_coords
        )
    ]).as_result()
}