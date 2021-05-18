use std::ops::{Add, AddAssign, Div, Index, Mul, MulAssign, Sub, SubAssign};

use crate::foundation::geometry::vector::Vector3;
use crate::foundation::pbr::Float;

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Point3 {
    pub x: Float,
    pub y: Float,
    pub z: Float,
}

impl Point3 {
    pub fn abs(self) -> Self {
        Point3 {
            x: self.x.abs(),
            y: self.y.abs(),
            z: self.z.abs(),
        }
    }
}

impl Add for Point3 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl Add<Vector3> for Point3 {
    type Output = Self;

    fn add(self, other: Vector3) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl Sub for Point3 {
    type Output = Vector3;

    fn sub(self, other: Self) -> Self::Output {
        Vector3 {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl Sub<Vector3> for Point3 {
    type Output = Point3;

    fn sub(self, other: Vector3) -> Self::Output {
        Point3 {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl Index<i32> for Point3 {
    type Output = Float;

    fn index(&self, index: i32) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Index out of bounds: {}", index),
        }
    }
}

impl Mul<Float> for Point3 {
    type Output = Self;

    fn mul(self, rhs: Float) -> Self::Output {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl MulAssign<Float> for Point3 {
    fn mul_assign(&mut self, rhs: Float) {
        self.x = self.x * rhs;
        self.y = self.y * rhs;
        self.z = self.z * rhs;
    }
}

impl Div<Float> for Point3 {
    type Output = Self;

    fn div(self, rhs: Float) -> Self::Output {
        Self {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

impl Point3 {
    pub fn zero() -> Point3 {
        Point3 {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    pub fn new(x: Float, y: Float, z: Float) -> Self {
        Self { x, y, z }
    }

    pub fn length_squared(self) -> Float {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    pub fn length(self) -> Float {
        self.length_squared().sqrt()
    }

    pub fn distance_squared(self, other: Self) -> Float {
        (self - other).length_squared()
    }

    pub fn distance(self, other: Self) -> Float {
        (self - other).length()
    }

    pub fn min(self, other: Self) -> Self {
        Self {
            x: Float::min(self.x, other.x),
            y: Float::min(self.y, other.y),
            z: Float::min(self.z, other.z),
        }
    }

    pub fn max(self, other: Self) -> Self {
        Self {
            x: Float::max(self.x, other.x),
            y: Float::max(self.y, other.y),
            z: Float::max(self.z, other.z),
        }
    }

    pub fn normalize(&self) -> Self {
        *self / self.length()
    }

    // pub fn lerp(self) ->
}

pub struct Point2 {
    x: Float,
    y: Float,
}

impl Point2 {
    pub fn new(x: Float, y: Float) -> Self {
        Self { x, y }
    }
}

impl Add for Point2 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl AddAssign for Point2 {
    fn add_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl Sub for Point2 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl SubAssign for Point2 {
    fn sub_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl Point2 {
    pub fn length_squared(self) -> Float {
        self.x * self.x + self.y * self.y
    }

    pub fn length(self) -> Float {
        self.length_squared().sqrt()
    }

    pub fn distance_squared(self, other: Self) -> Float {
        (self - other).length_squared()
    }

    pub fn distance(self, other: Self) -> Float {
        (self - other).length()
    }

    pub fn min(self, other: Self) -> Self {
        Self {
            x: Float::min(self.x, other.x),
            y: Float::min(self.y, other.y),
        }
    }

    pub fn max(self, other: Self) -> Self {
        Self {
            x: Float::max(self.x, other.x),
            y: Float::max(self.y, other.y),
        }
    }
}

impl Index<i32> for Point2 {
    type Output = Float;

    fn index(&self, index: i32) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            _ => panic!("Index out of bounds: {}", index),
        }
    }
}

fn lerp(t: Float, v1: Float, v2: Float) -> Float {
    (1.0 - t) * v1 + t * v2
}
