#[derive(Debug, Copy, Clone)]
pub struct  Point {
    x : f32,
    y : f32,
    z : f32
}

impl Point {
    pub fn scale(&self, factor: f32) -> Point {
        return Point {
            x: self.x * factor,
            y: self.y * factor,
            z: self.z * factor
        }
    }
}

pub struct Sphere {
    center: Point,
    radius : f32
}

impl Sphere {
    const ZERO : Sphere = Sphere { center: Point { x: 0.0, y: 0.0, z: 0.0 }, radius: 0.0};
}

pub fn distance(
    Point { x : x1, y : y1, z: z1 }: &Point,
    Point { x : x2, y : y2, z: z2 }: &Point) -> f32 {

    return (x2 - x1).powf(2.0) + (y2 - y1).powf(2.0) + (z2 - z1).powf(2.0)
}

pub struct Range3d {
    range_x : (f32, f32),
    range_y : (f32, f32),
    range_z : (f32, f32)
}

pub fn max_f32(v1 : f32, v2 : f32) -> f32 {
    if v2 > v1 {
        return v2
    }

    return v1
}

pub fn min_f32(v1 : f32, v2 : f32) -> f32 {
    if v2 < v1 {
        return v2
    }
    return v1
}

impl Range3d {

    pub fn far_end(&self) -> Point {
        return Point {
            x : self.range_x.1,
            y : self.range_y.1,
            z : self.range_z.1
        }
    }

    pub fn widen(&mut self, x: f32, y: f32, z:f32) {

        self.range_x = (min_f32(self.range_x.0, x), max_f32(self.range_x.1, x));
        self.range_y = (min_f32(self.range_y.0, y), max_f32(self.range_y.1, y));
        self.range_z = (min_f32(self.range_z.0, z), max_f32(self.range_z.1, z))
    }

    pub fn from_point(Point{ x,y,z }: Point) -> Range3d {
        return Range3d {
            range_x: (x,x), range_y: (y,y), range_z: (z,z) }
    }

    pub fn get_range(points : &Vec<Point>) -> Range3d {

        assert!(points.len() > 0);

        let mut result = Range3d::from_point(points[0]);

        for &Point{ x,y,z } in points {
            result.widen(x,y,z);
        }

        return result
    }

    pub fn spread(&self) -> f32 {
        let squares = (self.range_x.1 - self.range_x.0).powf(2.0)
            + (self.range_y.1 - self.range_y.0).powf(2.0)
            + (self.range_z.1 - self.range_z.0).powf(2.0);

        return squares.sqrt()
    }
}

pub fn bounding_sphere(points : Vec<Point>) -> Sphere {

    if points.len() == 0 {
        return Sphere::ZERO;
    }
    else if points.len() == 1 {
        return Sphere { center : points[0], radius: 0.0 }
    }

    let range = Range3d::get_range(&points);
    let large_spread = range.spread() * 100.0;
    let far_away = range.far_end().scale(large_spread);

    return Sphere::ZERO;
}

pub fn sum(a: usize, b: usize) -> usize {
    return a + b
}