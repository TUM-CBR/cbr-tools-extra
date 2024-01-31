pub type Point3d = (f64, f64, f64);

pub type Matrix3d = (Point3d, Point3d, Point3d);

pub const IDENTITY_3D : Matrix3d = (
    (1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0)
);

pub trait Functor<T> {

    fn fmap<F>(&self, op: F) -> Self
        where F : Fn(T) -> T;
}

impl<T> Functor<T> for (T, T, T) {

    fn fmap<F>(&self, op : F) -> (T, T, T)
        where F : Fn(T) -> T {
        let (x,y,z) = self;
        (
            op(*x),
            op(*y),
            op(*z)
        )
    }
}

pub trait HasClosedAdd where Self : Sized {
    fn plus(&self, other: &Self) -> Self;

    fn neg(&self) -> Self;

    fn minus(&self, other: &Self) -> Self {
        self.plus(&other.neg())
    }
}

impl HasClosedAdd for Point3d {

    fn plus(&self, other: &Point3d) -> Point3d {
        let (x1, y1, z1) = self;
        let (x2, y2, z2) = other;
    
        (x1 + x2, y1 + y2, z1 + z2)
    }

    fn neg(&self) -> Point3d {
        self.fmap(|v| { -v })
    }
}

impl HasClosedAdd for Matrix3d {

    fn plus(&self, other: &Self) -> Self {
        let (a1, a2, a3) = self;
        let (b1, b2, b3) = other;

        (a1.plus(b1), a2.plus(b2), a3.plus(b3))
    }

    fn neg(&self) -> Self {
        self.fmap(|v| { v.neg()})
    }
}

pub trait HasMultWith<T> {

    fn mult(&self, value: T) -> Self;
}

impl HasMultWith<f64> for Point3d {
    fn mult(&self, scalar: f64) -> Point3d {
        let (x,y,z) = self;
    
        (x*scalar, y*scalar, z*scalar)
    }
}

impl HasMultWith<f64> for Matrix3d {
    fn mult(&self, scalar: f64) -> Self {
        let (c1, c2, c3) = self;

        (c1.mult(scalar), c2.mult(scalar), c3.mult(scalar))
    }
}

impl HasMultWith<&Matrix3d> for Point3d {

    fn mult(&self, value: &Matrix3d) -> Self {
        let (c1, c2, c3) = value;

        (
            self.inner(c1),
            self.inner(c2),
            self.inner(c3)
        )
    }
}

pub trait HasInnerProduct where Self : Sized {

    fn inner(&self, other: &Self) -> f64;

    fn norm(&self) -> f64 {
        self.inner(self).sqrt()
    }
}

impl HasInnerProduct for Point3d {

    fn inner(&self, other: &Point3d) -> f64 {
        let (x1, y1, z1) = self;
        let (x2, y2, z2) = other;

        return x1*x2 + y1*y2 + z1*z2
    }
}

pub trait HasOutterProduct where Self : Sized {
    type TOut : Sized;

    fn outter(&self, other: &Self) -> Self::TOut;
}

impl HasOutterProduct for Point3d {

    type TOut = Matrix3d;

    fn outter(&self, other: &Self) -> Self::TOut {
        let (x1, y1, z1) = self;
        let (x2, y2, z2) = self;

        let c1 = (x1*x2, y1*x2, z1*x2);
        let c2 = (x1*y2, y1*y2, z1*y2);
        let c3 = (x1*z2, y1*z2, z1*z2);

        (c1,c2,c3)
    }
}

pub trait HasCrossProduct where Self : Sized {

    fn cross(&self, other: &Self) -> Self;
}

impl HasCrossProduct for Point3d {

    fn cross(&self, other: &Self) -> Self {
        let (a1, a2, a3) = self;
        let (b1, b2, b3) = other;

        (
            a2*b3 - a3*b2,
            a3*b1 - a1*b3,
            a1*b2 - a2*b1
        )
    }
}

pub trait ScalarVector
    where
        Self : HasClosedAdd,
        Self : HasMultWith<f64>,
        Self : HasCrossProduct,
        Self : HasInnerProduct {

    fn unit(&self) -> Self {
        self.mult(1.0 / self.norm())
    }

    fn angle(&self, other: Self) -> f64 {

        let num = self.inner(&other);
        let den = self.norm()*other.norm();

        (num/den).acos()
    }
}

impl ScalarVector for Point3d {}

pub trait HasRotation {

    type RotationMatrix;

    /// Computes the rotation matrix rotate points using
    /// this vector as the axis.
    fn rot_matrix(&self, angle: f64) -> Self::RotationMatrix;

    /// Rotate the vector using the given rotation matrix
    fn apply_rot(&self, matrix: Self::RotationMatrix) -> Self;
}

pub fn skew_symmetric_matrix_3d(vec: &Point3d) -> Matrix3d {
    let (x, y, z) = vec;
    (
        (0.0, -z, *y),
        (*z, 0.0, -x),
        (-y, *x, 0.0)
    )
}

impl HasRotation for Point3d {

    type RotationMatrix = Matrix3d;

    fn rot_matrix(&self, angle: f64) -> Self::RotationMatrix {
        let unit_self = self.unit();
        let r1 = IDENTITY_3D.mult(angle.cos());
        let r2 = unit_self.outter(&unit_self).mult(1.0 - angle.cos());
        let r3 = skew_symmetric_matrix_3d(&unit_self).mult(angle.sin());

        r1.plus(&r2).plus(&r3)
    }

    fn apply_rot(&self, matrix: Self::RotationMatrix) -> Self {
        self.mult(&matrix)
    }
}

pub trait HasAbsolute {

    fn abs(&self) -> Self;
}

impl HasAbsolute for Point3d {

    fn abs(&self) -> Self {
        self.fmap(|v|{ v.abs() })

    }
}