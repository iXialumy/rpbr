use crate::foundation::geometry::vector::Vector3;
use num_traits::{AsPrimitive, Float, FromPrimitive};
use std::ops::Div;

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

impl<T: Float + AsPrimitive<f64> + FromPrimitive> Div<T> for Normal3<T> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        Self {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Normal3<T> {
    pub fn zero() -> Normal3<T> {
        Normal3 {
            x: T::zero(),
            y: T::zero(),
            z: T::zero(),
        }
    }
    pub fn length_squared(self) -> T {
        self.x * self.x + self.y * self.y + self.z * self.z
    }
    pub fn length(self) -> T {
        self.length_squared().sqrt()
    }
    pub fn normalize(self) -> Self {
        self / self.length()
    }
}
