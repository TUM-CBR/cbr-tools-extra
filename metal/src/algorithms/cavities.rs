use std::{collections::btree_map::Iter, ops::AddAssign};

use polars::prelude::*;

use crate::core::*;

use self::data::*;

pub struct PdbCavities {
    cavities : Vec<Cavity>
}

pub struct Cavity {
    points : DataFrame
}

/// Value used to iterate through all the points inside a box which
/// is delimited by two points.
struct BoxIterator {
    start: Point3d,
    end: Point3d,
    resolution: f64,
    rot_matrix: Matrix3d,
    xn: i64,
    yn: i64,
    zn: i64,
    xi: i64,
    yi: i64,
    zi: i64
}

impl BoxIterator {

    /// Construct the iterator. The approach looks convoluted, but it is rather
    /// intuitive once explained.
    /// 1. We determine how many units will be traversed along each axis. To
    ///    to this, we simply compute the difference between end and start
    ///    for each axis and divide by the resolution. We take the integer value.
    /// 2. We den compute a new box, which start at the origin (0,0,0) and goes
    ///    to the values computed in (1). The idea is that we will traverse
    ///    all the points in this box.
    /// 3. Convert the positions computed in (2) to actual positions in our box,
    ///    we must rotate each of the boints to match the potential rotation
    ///    relative to the original box and traslate them by adding the
    ///    lower corner of the original box.
    pub fn new(start: Point3d, end: Point3d, resolution: f64) -> BoxIterator {

        let (x_units, y_units, z_units) =
            end.minus(&start)
                .abs()
                .mult(1.0/resolution);
        let (xn, yn, zn) = end;
        let diagonal = end.minus(&start);
        let base_diagonal = (x_units, y_units, z_units);
        let rot_axis = base_diagonal.cross(&diagonal);
        let rot_matrix = rot_axis.rot_matrix(
            base_diagonal.angle(diagonal)
        );

        BoxIterator {
            start,
            end,
            resolution,
            rot_matrix,
            xn: x_units as i64,
            yn: y_units as i64,
            zn: z_units as i64,
            xi: 0,
            yi: 0,
            zi: 0
        }
    }

    pub fn current(&self) -> Point3d {

        let scale = self.resolution;
        let vector = (
            self.xi as f64*scale,
            self.yi as f64*scale,
            self.zi as f64*scale
        );
        vector.apply_rot(self.rot_matrix).plus(&self.start)
    }

    pub fn done(&self) -> bool {

        self.xi == self.xn && self.yi == self.yn && self.zi == self.zn
    }
}

fn inc_by(value: &mut i64, carry: i64, top: i64) -> i64 {
    let result = *value + carry;

    if carry > result {
        panic!("Carry should not be greater than result.")
    }

    if result <= top {
        *value = result;
        0
    }
    else {
        *value = 0;
        result - top
    }
}

impl Iterator for BoxIterator {
    type Item = Point3d;

    fn next(&mut self) -> Option<Point3d> {

        if self.done() {
            return None;
        }

        let result = self.current();
        let x_carry = inc_by(&mut self.xi, 1, self.xn);
        let y_carry = inc_by(&mut self.yi, x_carry, self.yn);
        let _ = inc_by(&mut self.zi, y_carry, self.zn);

        Some(result)
    }
}

/// Iterate over all of the points that exists inside the box delimited by the
/// two `Point3d` provided as argument. The resolution parameter indicates what
/// spacing will exists between the closest points.
pub fn traverse_box(
    start: Point3d,
    end: Point3d,
    resolution: f64
) -> impl Iterator<Item = Point3d> {

    BoxIterator::new(start, end, resolution)
}

impl PdbCavities {

    pub fn find_cavities(pdb_context: &PdbContext) -> PdbCavities {

        let atoms = pdb_context.all_atoms();
        let bounds = pdb_context.bounding_box();
    }
}
