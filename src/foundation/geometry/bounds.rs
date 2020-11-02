use std::ops::Index;

use num_traits::Float;

use crate::foundation::geometry::point::Point3;
use crate::foundation::geometry::vector::Vector3;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Bounds3<T: Float> {
    p_min: Point3<T>,
    p_max: Point3<T>,
}

impl<T: Float> Index<i32> for Bounds3<T> {
    type Output = Point3<T>;

    fn index(&self, index: i32) -> &Self::Output {
        match index {
            0 => &self.p_min,
            1 => &self.p_max,
            _ => panic!("Index {} out of bounds for Bounds", index),
        }
    }
}

impl<T: Float> Bounds3<T> {
    pub fn new(p1: Point3<T>, p2: Point3<T>) -> Self {
        Self {
            p_min: Point3 {
                x: Float::min(p1.x, p2.x),
                y: Float::min(p1.y, p2.y),
                z: Float::min(p1.z, p2.z),
            },
            p_max: Point3 {
                x: Float::max(p1.x, p2.x),
                y: Float::max(p1.y, p2.y),
                z: Float::max(p1.z, p2.z),
            },
        }
    }

    /// Returns one of the eight corners Points
    pub fn corner(self, corner: i32) -> Point3<T> {
        debug_assert!(corner >= 0 && corner < 8);
        let x = corner & 1;
        let y = if (corner & 2) == 0 { 0 } else { 1 };
        let z = if (corner & 4) == 0 { 0 } else { 1 };

        Point3 {
            x: self[x].x,
            y: self[y].y,
            z: self[z].z,
        }
    }

    /// Returns a box in which both bounding boxes fit
    pub fn union(self, other: Self) -> Self {
        Self {
            p_min: Point3::min(self.p_min, other.p_min),
            p_max: Point3::max(self.p_max, other.p_max),
        }
    }

    pub fn overlaps(self, other: Self) -> bool {
        let x = (self.p_max.x >= other.p_min.x) && (self.p_min.x <= self.p_max.x);
        let y = (self.p_max.y >= other.p_min.y) && (self.p_min.y <= self.p_max.y);
        let z = (self.p_max.z >= other.p_min.z) && (self.p_min.z <= self.p_max.z);
        x && y && z
    }

    /// Returns true if the Point is inside the Bounds
    pub fn inside(self, point: Point3<T>) -> bool {
        point.x >= self.p_min.x
            && point.x <= self.p_max.x
            && point.y >= self.p_min.y
            && point.y <= self.p_max.y
            && point.z >= self.p_min.z
            && point.z <= self.p_max.z
    }

    /// Returns true if the Point is inside the Bounds, excluding points at the upper bounds
    pub fn inside_exclusive(self, point: Point3<T>) -> bool {
        point.x >= self.p_min.x
            && point.x < self.p_max.x
            && point.y >= self.p_min.y
            && point.y < self.p_max.y
            && point.z >= self.p_min.z
            && point.z < self.p_max.z
    }

    /// Returns Bounds padded by a constant factor
    pub fn expand(self, delta: T) -> Self {
        Self {
            p_min: self.p_min
                - Point3 {
                    x: delta,
                    y: delta,
                    z: delta,
                },
            p_max: self.p_max
                + Point3 {
                    x: delta,
                    y: delta,
                    z: delta,
                },
        }
    }
}

impl<T: Float> Bounds3<T> {
    pub fn diagonal(self) -> Vector3<T> {
        (self.p_max - self.p_min).into()
    }

    pub fn maximum_extent(self) -> i8 {
        let diag = self.diagonal();
        diag.max_dimension()
    }
}

#[cfg(test)]
mod tests {
    use crate::foundation::geometry::bounds::Bounds3;
    use crate::foundation::geometry::point::Point3;

    #[test]
    fn corner() {
        let p1 = Point3 {
            x: 1.0,
            y: 2.0,
            z: 3.0,
        };
        let p2 = Point3 {
            x: 3.0,
            y: 2.0,
            z: 1.0,
        };
        let b = Bounds3::new(p1, p2);

        let expected = Point3 {
            x: 3.0,
            y: 2.0,
            z: 1.0,
        };

        assert_eq!(b.corner(1), expected)
    }

    #[test]
    fn union() {
        let p1 = Point3 {
            x: 1.0,
            y: 3.0,
            z: 6.0,
        };
        let p2 = Point3 {
            x: 4.0,
            y: 2.0,
            z: 6.0,
        };
        let p3 = Point3 {
            x: 3.0,
            y: 7.0,
            z: 2.0,
        };
        let p4 = Point3 {
            x: 8.0,
            y: 2.0,
            z: 7.0,
        };
        let b1 = Bounds3::new(p1, p2);
        let b2 = Bounds3::new(p3, p4);

        let p_min = Point3 {
            x: 1.0,
            y: 2.0,
            z: 2.0,
        };
        let p_max = Point3 {
            x: 8.0,
            y: 7.0,
            z: 7.0,
        };
        let expected = Bounds3 { p_min, p_max };

        assert_eq!(expected, b1.union(b2));
    }
}
