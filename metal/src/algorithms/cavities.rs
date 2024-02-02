use std::{collections::btree_map::Iter, ops::AddAssign};

use polars::prelude::*;

use crate::core::*;

pub struct PdbCavities {
    cavities : Vec<Cavity>
}

pub struct Cavity {
    points : DataFrame
}



impl PdbCavities {

    pub fn find_cavities(pdb_context: &PdbContext) -> PdbCavities {

        let atoms = pdb_context.all_atoms();
        let bounds = pdb_context.bounding_box();

        panic!("Nooo!!!")
    }
}