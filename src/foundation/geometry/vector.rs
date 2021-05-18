use std::ops::{Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::foundation::geometry::point::Point3;
use crate::foundation::pbr::Float;

pub fn max(first: Float, second: Float) -> Float {
    if first < second {
        second
    } else {
        first
    }
}

pub fn min(first: Float, second: Float) -> Float {
    if first < second {
        first
    } else {
        second
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Vector3 {
    pub x: Float,
    pub y: Float,
    pub z: Float,
}

impl Index<i32> for Vector3 {
    type Output = Float;

    fn index(&self, index: i32) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Index out of bounds: {} for Vector3", index),
        }
    }
}

impl Add for Vector3 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl AddAssign for Vector3 {
    fn add_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl Sub for Vector3 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl SubAssign for Vector3 {
    fn sub_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl Mul<Float> for Vector3 {
    type Output = Self;

    fn mul(self, other: Float) -> Self {
        Self {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl MulAssign<Float> for Vector3 {
    fn mul_assign(&mut self, other: Float) {
        *self = Self {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl Div<Float> for Vector3 {
    type Output = Self;

    fn div(self, other: Float) -> Self {
        Self {
            x: self.x / other,
            y: self.y / other,
            z: self.z / other,
        }
    }
}

impl DivAssign<Float> for Vector3 {
    fn div_assign(&mut self, rhs: Float) {
        *self = Self {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

impl Neg for Vector3 {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl From<Point3> for Vector3 {
    fn from(point: Point3) -> Self {
        Self {
            x: point.x,
            y: point.y,
            z: point.z,
        }
    }
}

#[allow(dead_code)]
impl Vector3 {
    pub fn new(x: Float, y: Float, z: Float) -> Self {
        debug_assert!(x.is_nan());
        debug_assert!(y.is_nan());
        debug_assert!(z.is_nan());

        Self { x, y, z }
    }

    pub fn empty() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    pub fn length_squared(self) -> Float {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    pub fn length(self) -> Float {
        Float::sqrt(self.length_squared())
    }

    pub fn abs(self) -> Self {
        Self {
            x: Float::abs(self.x),
            y: Float::abs(self.y),
            z: Float::abs(self.z),
        }
    }

    pub fn dot(&self, other: Self) -> Float {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn abs_dot(&self, other: Self) -> Float {
        Float::abs(self.dot(other))
    }

    pub fn normalize(self) -> Self {
        self / self.length()
    }

    /// Returns the largest of x, y and z
    /// Example:
    ///
    /// ```
    /// use rpbr::foundation::geometry::vector::Vector3;
    ///
    /// let v = Vector3 {
    ///     x: 1.0,
    ///     y: 3.0,
    ///     z: 2.0,
    /// };
    ///
    /// assert_eq!(3.0, v.max_component());
    /// ```
    pub fn max_component(self) -> Float {
        match (self.x > self.y, self.x > self.z, self.y > self.z) {
            (true, true, _) => self.x,
            (true, false, _) => self.z,
            (false, _, true) => self.y,
            (false, _, false) => self.z,
        }
    }

    pub fn min_component(self) -> Float {
        match (self.x < self.y, self.x < self.z, self.y < self.z) {
            (true, true, _) => self.x,
            (true, false, _) => self.z,
            (false, _, true) => self.y,
            (false, _, false) => self.z,
        }
    }

    pub fn max_dimension(self) -> i8 {
        match (self.x > self.y, self.x > self.z, self.y > self.z) {
            (true, true, _) => 0,
            (true, false, _) => 1,
            (false, _, true) => 2,
            (false, _, false) => 1,
        }
    }

    pub fn min(self, other: Self) -> Self {
        Self {
            x: min(self.x, other.x),
            y: min(self.y, other.y),
            z: min(self.z, other.z),
        }
    }

    pub fn max(self, other: Self) -> Self {
        Self {
            x: max(self.x, other.x),
            y: max(self.y, other.y),
            z: max(self.z, other.z),
        }
    }

    pub fn cross(&self, other: Self) -> Self {
        let v1x: f64 = self.x as f64;
        let v2x: f64 = other.x as f64;
        let v1y: f64 = self.y as f64;
        let v2y: f64 = other.y as f64;
        let v1z: f64 = self.z as f64;
        let v2z: f64 = other.z as f64;

        Vector3 {
            x: ((v1y * v2z) - (v1z * v2y)) as Float,
            y: ((v1z * v2x) - (v1x * v2z)) as Float,
            z: ((v1x * v2y) - (v1y * v2x)) as Float,
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub(crate) struct Vector2 {
    pub x: Float,
    pub y: Float,
}

impl Add for Vector2 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl AddAssign for Vector2 {
    fn add_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl Sub for Vector2 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl SubAssign for Vector2 {
    fn sub_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl Mul<Float> for Vector2 {
    type Output = Self;

    fn mul(self, other: Float) -> Self {
        Self {
            x: self.x * other,
            y: self.y * other,
        }
    }
}

impl MulAssign<Float> for Vector2 {
    fn mul_assign(&mut self, other: Float) {
        *self = Self {
            x: self.x * other,
            y: self.y * other,
        }
    }
}

impl Div<Float> for Vector2 {
    type Output = Self;

    fn div(self, other: Float) -> Self {
        Self {
            x: self.x / other,
            y: self.y / other,
        }
    }
}

impl DivAssign<Float> for Vector2 {
    fn div_assign(&mut self, rhs: Float) {
        *self = Self {
            x: self.x / rhs,
            y: self.y / rhs,
        }
    }
}

impl Neg for Vector2 {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
        }
    }
}

#[allow(dead_code)]
impl Vector2 {
    pub fn length_squared(&self) -> Float {
        self.x * self.x + self.y * self.y
    }

    pub fn length(&self) -> Float {
        Float::sqrt(self.length_squared())
    }

    pub fn abs(self) -> Self {
        Self {
            x: Float::abs(self.x),
            y: Float::abs(self.y),
        }
    }

    pub fn dot(&self, other: Self) -> Float {
        self.x * other.x + self.y * other.y
    }

    pub fn abs_dot(&self, other: Self) -> Float {
        Float::abs(self.dot(other))
    }

    pub fn normalize(self) -> Self {
        self / self.length()
    }
}
