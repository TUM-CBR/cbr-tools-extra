use std::fmt::Pointer;

use crate::core::data::*;

/// Value used to iterate through all the points inside a box which
/// is delimited by two points.
pub struct BoxIterator {
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

pub trait ClosedShape3D {

    type TVolumeIterator : Iterator<Item = Point3d>;

    fn iterate_volume(&self, resolution: f64) -> Self::TVolumeIterator;

    fn inside(&self, point: Point3d) -> bool;
}

pub trait HasTraslate {

    fn traslate(&self, direction: &Point3d) -> Self;
}

impl HasTraslate for Box3D {

    fn traslate(&self, direction: &Point3d) -> Self {

        self.fmap(|p| { p.plus(direction)} )
    }
}

pub trait HasScale {

    fn scale(&self, how_much: f64) -> Self;
}

impl HasScale for Box3D {
    
    fn scale(&self, how_much: f64) -> Self {
        self.fmap(|p| { p.mult(how_much)})
    }
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
        let rot_matrix = <Point3d as HasRotation>::rot_matrix(
            &rot_axis,
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
) -> BoxIterator {

    BoxIterator::new(start, end, resolution)
}

impl ClosedShape3D for Box3D {

    type TVolumeIterator = BoxIterator;

    fn iterate_volume(&self, resolution: f64) -> Self::TVolumeIterator {
        traverse_box(self.bottom, self.top, resolution)
    }

    fn inside(&self, point: Point3d) -> bool {

        let (x0,y0,z0) = self.bottom;
        let (xn, yn, zn) = self.top;
        let (x,y,z) = point;

        x >= x0 && y >= x0 && z >= z0 && x <= xn && y <= yn && z <= zn
    }
}

#[cfg(test)]
mod tests {

    use crate::{algorithms::spatial::{ClosedShape3D, HasScale, HasTraslate}, core::data::*};

    struct ArbitraryValues {
        number1 : i64,
        number2 : i64
    }

    impl ArbitraryValues {
        pub fn new() -> ArbitraryValues {

            ArbitraryValues {
                number1: 503,
                number2: 503*51 + 2
            }
        }
        pub fn next(&mut self, count: usize) -> Vec<f64> {
            let big_prime: i64 = 10007;
            let mut result = Vec::with_capacity(count);
            let mut number1: i64 = self.number1;
            let mut number2: i64 = self.number2;
    
            while result.len() < count {
                let next = (number1 * number2 + 2) % big_prime;
                number1 = number2;
                number2 = next;
                result.push(next as f64);
            }
            self.number1 = number1;
            self.number2 = number2;
            result
        }
    }

    fn apply_transforms(
        value: &Point3d,
        x_rot: f64,
        y_rot: f64,
        scale: f64,
        offset: &Point3d
    ) -> Point3d {

        value.rotate(x_rot, &(1.0,0.0,0.0))
            .rotate(y_rot, &(0.0, 1.0, 0.0))
            .mult(scale)
            .plus(&offset)
    }

    #[test]
    fn enumerates_all_points() {

        let mut arbitrary = ArbitraryValues::new();
        let box_base = Box3D::create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let angles : Vec<f64> = (2..10)
            .map(|den| { 3.14159 / den as f64})
            .collect();


        let scales : Vec<f64> = arbitrary.next(10);
        let x_offset = arbitrary.next(10);
        let y_offset = arbitrary.next(10);
        let z_offset = arbitrary.next(10);
        let offsets : Vec<(f64, f64, f64)> =
            x_offset.into_iter()
                .zip(y_offset)
                .zip_n(z_offset)
                .collect();

        let mut canary = 0;

        for angle_x in angles.clone() {
            for angle_z in angles.clone(){
                for scale in scales.clone() {
                    for offset in offsets.clone() {

                        let t_box = box_base
                            .rotate(angle_x,&(1.0,0.0,0.0))
                            .rotate(angle_z, &(0.0, 0.0, 1.0))
                            .scale(scale)
                            .traslate(&offset);

                        for point in t_box.iterate_volume(1.0){
                            canary += 1;

                            assert!(t_box.inside(point))
                        }

                    }
                }
            }
        }

        assert_eq!(8, 8)
    }
}