use pdbtbx::{AtomConformerResidueChainModelMut, ContainsAtomConformer, PDB};
use polars::prelude::*;
use std::collections::LinkedList;

enum PdbReadArgs {
    Default,
}

enum PdbColumn {
    X_COORD,
    Y_COORD,
    Z_COORD
}

impl From<PdbColumn> for &str {

    fn from(column: PdbColumn) -> &'static str {

        match column {
            PdbColumn::X_COORD => "x_coord",
            PdbColumn::Y_COORD => "y_coord",
            PdbColumn::Z_COORD => "z_coord"
        }
    }
}

pub fn from_pdb(
    args: PdbReadArgs,
    pdb: &PDB
) -> DataFrame {

    //let atoms = Series::new("atoms", pdb.atoms_with_hierarchy_mut().into());
    let atom_count = pdb.atom_count();
    let mut x_coords : Vec<f64> = Vec::with_capacity(atom_count);
    let mut y_coords : Vec<f64> = Vec::with_capacity(atom_count);
    let mut z_coords : Vec<f64> = Vec::with_capacity(atom_count);

    for entry in pdb.atoms_with_hierarchy() {
        let atom = entry.atom();
        x_coords.push(atom.x());
        y_coords.push(*atom.y());
        z_coords.push(*atom.z());
    }

    DataFrame::new(vec![
        Series::new(
            PdbColumn::X_COORD.into(),
            x_coords
        )
    ])
}