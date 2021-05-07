use crate::foundation::geometry::vector::Vector3;
use num_traits::{AsPrimitive, Float, FromPrimitive};

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Normal3<T: Float + FromPrimitive + AsPrimitive<f64> + Copy> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T: Float + AsPrimitive<f64> + FromPrimitive> From<Vector3<T>> for Normal3<T> {
    fn from(vector: Vector3<T>) -> Self {
        Normal3 {
            x: vector.x,
            y: vector.y,
            z: vector.z,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Normal3<T> {
    pub fn default() -> Normal3<T> {
        Normal3 {
            x: T::zero(),
            y: T::zero(),
            z: T::zero(),
        }
    }
}
