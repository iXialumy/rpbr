use std::ops::Index;

use num_traits::{AsPrimitive, Float, FromPrimitive};

use crate::foundation::geometry::point::Point3;
use crate::foundation::geometry::ray::Ray;
use crate::foundation::geometry::vector::Vector3;
use std::mem::swap;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Bounds3<T: Float + FromPrimitive + AsPrimitive<f64>> {
    pub p_min: Point3<T>,
    pub p_max: Point3<T>,
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Bounds3<T> {
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
        debug_assert!((0..8).contains(&corner));
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

    /// Returns a box, in which the box and the point will fit
    pub fn union_point(&self, point: Point3<T>) -> Bounds3<T> {
        Self {
            p_min: Point3::min(self.p_min, point),
            p_max: Point3::max(self.p_max, point),
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
                - Vector3 {
                    x: delta,
                    y: delta,
                    z: delta,
                },
            p_max: self.p_max
                + Vector3 {
                    x: delta,
                    y: delta,
                    z: delta,
                },
        }
    }

    pub fn diagonal(self) -> Vector3<T> {
        self.p_max - self.p_min
    }

    pub fn maximum_extent(self) -> i8 {
        let diag = self.diagonal();
        diag.max_dimension()
    }

    pub fn intersect_p(self, ray: Ray<T>) -> Option<(T, T)> {
        let mut t0 = T::from_f64(0.0).unwrap();
        let mut t1 = ray.max_length;

        for i in 0..3 {
            let inv_ray_dir = T::from_f64(1.0).unwrap() / ray.direction[i];
            let mut t_near = (self.p_min[i] - ray.origin[i]) * inv_ray_dir;
            let mut t_far = (self.p_max[i] - ray.origin[i]) * inv_ray_dir;

            if t_near > t_far {
                swap(&mut t_near, &mut t_far);
            }

            // TODO comment in when gamma will be implemented
            // t_far = t_far * (T::from_f64(1.0).unwrap() + T::from_f64(2.0).unwrap() * gamma(T::from_f64(3.0).unwrap()));
            t0 = if t_near > t0 { t_near } else { t0 };
            t1 = if t_far < t1 { t_far } else { t1 };

            if t0 > t1 {
                return None;
            }
        }

        Some((t0, t1))
    }

    pub fn intersect_p_precomp(
        self,
        ray: Ray<T>,
        inv_dir: &Vector3<T>,
        dir_is_neg: Box<[i32]>,
    ) -> bool {
        let mut t_min = (self[dir_is_neg[0]].x - ray.origin.x) * inv_dir.x;
        let mut t_max = (self[1 - dir_is_neg[0]].x - ray.origin.x) * inv_dir.x;
        let mut ty_min = (self[dir_is_neg[1]].y - ray.origin.y) * inv_dir.y;
        let mut ty_max = (self[1 - dir_is_neg[1]].y - ray.origin.y) * inv_dir.y;

        // TODO comment in when gamma will be implemented
        // t_max = t_max * (T::from_f64(1.0).unwrap() + T::from_f64(2.0).unwrap() * gamma(T::from_f64(3.0).unwrap()));
        // ty_max = ty_max * (T::from_f64(1.0).unwrap() + T::from_f64(2.0).unwrap() * gamma(T::from_f64(3.0).unwrap()));

        if t_min > t_max || ty_min > t_max {
            return false;
        }
        if ty_min > t_min {
            t_min = ty_min;
        }
        if ty_max < t_max {
            t_max = ty_max;
        }

        let mut tz_min = (self[dir_is_neg[2]].z - ray.origin.z) * inv_dir.z;
        let mut tz_max = (self[1 - dir_is_neg[2]].z - ray.origin.z) * inv_dir.z;

        // TODO comment in when gamma will be implemented
        // tz_max = tz_max * (T::from_f64(1.0).unwrap() + T::from_f64(2.0).unwrap() * gamma(T::from_f64(3.0).unwrap()));

        if t_min > tz_max || tz_min > t_max {
            return false;
        }
        if tz_min > t_min {
            t_min = tz_min;
        }
        if tz_max < t_max {
            t_max = tz_max;
        }

        (t_min < ray.max_length) && (t_max > T::zero())
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Index<i32> for Bounds3<T> {
    type Output = Point3<T>;

    fn index(&self, index: i32) -> &Self::Output {
        match index {
            0 => &self.p_min,
            1 => &self.p_max,
            _ => panic!("Index {} out of bounds for Bounds", index),
        }
    }
}

impl<T: Float + AsPrimitive<f64> + FromPrimitive> From<Point3<T>> for Bounds3<T> {
    fn from(point: Point3<T>) -> Self {
        Bounds3 {
            p_min: point,
            p_max: point,
        }
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
