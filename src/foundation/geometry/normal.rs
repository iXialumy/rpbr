use crate::foundation::geometry::vector::Vector3;
use crate::foundation::pbr::Float;
use std::ops::Div;

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Normal3 {
    pub x: Float,
    pub y: Float,
    pub z: Float,
}

impl From<Vector3> for Normal3 {
    fn from(vector: Vector3) -> Self {
        Normal3 {
            x: vector.x,
            y: vector.y,
            z: vector.z,
        }
    }
}

impl Div<Float> for Normal3 {
    type Output = Self;

    fn div(self, rhs: Float) -> Self::Output {
        Self {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

impl Normal3 {
    pub fn zero() -> Normal3 {
        Normal3 {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }
    pub fn length_squared(self) -> Float {
        self.x * self.x + self.y * self.y + self.z * self.z
    }
    pub fn length(self) -> Float {
        self.length_squared().sqrt()
    }
    pub fn normalize(self) -> Self {
        self / self.length()
    }
}
