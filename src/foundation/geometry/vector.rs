use std::ops::{Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Neg, Sub, SubAssign};

use num_traits::{AsPrimitive, Float, FromPrimitive};

use crate::foundation::geometry::point::Point3;

pub fn max<T: Float + FromPrimitive + AsPrimitive<f64>>(first: T, second: T) -> T {
    if first < second {
        second
    } else {
        first
    }
}

pub fn min<T: Float + FromPrimitive + AsPrimitive<f64>>(first: T, second: T) -> T {
    if first < second {
        first
    } else {
        second
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Vector3<T: Float + FromPrimitive + AsPrimitive<f64> + Copy> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T: Float + FromPrimitive + AsPrimitive<f64> + Copy> Index<i32> for Vector3<T> {
    type Output = T;

    fn index(&self, index: i32) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Index out of bounds: {} for Vector3", index),
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Add for Vector3<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> AddAssign for Vector3<T> {
    fn add_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Sub for Vector3<T> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> SubAssign for Vector3<T> {
    fn sub_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Mul<T> for Vector3<T> {
    type Output = Self;

    fn mul(self, other: T) -> Self {
        Self {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> MulAssign<T> for Vector3<T> {
    fn mul_assign(&mut self, other: T) {
        *self = Self {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Div<T> for Vector3<T> {
    type Output = Self;

    fn div(self, other: T) -> Self {
        Self {
            x: self.x / other,
            y: self.y / other,
            z: self.z / other,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> DivAssign<T> for Vector3<T> {
    fn div_assign(&mut self, rhs: T) {
        *self = Self {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Neg for Vector3<T> {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> From<Point3<T>> for Vector3<T> {
    fn from(point: Point3<T>) -> Self {
        Self {
            x: point.x,
            y: point.y,
            z: point.z,
        }
    }
}

#[allow(dead_code)]
impl<T: Float + FromPrimitive + AsPrimitive<f64>> Vector3<T> {
    pub fn new(x: T, y: T, z: T) -> Self {
        debug_assert!(x.is_nan());
        debug_assert!(y.is_nan());
        debug_assert!(z.is_nan());

        Self { x, y, z }
    }

    pub fn empty() -> Self {
        Self {
            x: T::from_f64(0.0).unwrap(),
            y: T::from_f64(0.0).unwrap(),
            z: T::from_f64(0.0).unwrap(),
        }
    }

    pub fn length_squared(self) -> T {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    pub fn length(self) -> T {
        Float::sqrt(self.length_squared())
    }

    pub fn abs(self) -> Self {
        Self {
            x: Float::abs(self.x),
            y: Float::abs(self.y),
            z: Float::abs(self.z),
        }
    }

    pub fn dot(&self, other: Self) -> T {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn abs_dot(&self, other: Self) -> T {
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
    pub fn max_component(self) -> T {
        match (self.x > self.y, self.x > self.z, self.y > self.z) {
            (true, true, _) => self.x,
            (true, false, _) => self.z,
            (false, _, true) => self.y,
            (false, _, false) => self.z,
        }
    }

    pub fn min_component(self) -> T {
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
        let v1x: f64 = (self.x).as_();
        let v2x: f64 = other.x.as_();
        let v1y: f64 = self.y.as_();
        let v2y: f64 = other.y.as_();
        let v1z: f64 = self.z.as_();
        let v2z: f64 = other.z.as_();

        Vector3 {
            x: T::from_f64((v1y * v2z) - (v1z * v2y)).unwrap(),
            y: T::from_f64((v1z * v2x) - (v1x * v2z)).unwrap(),
            z: T::from_f64((v1x * v2y) - (v1y * v2x)).unwrap(),
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub(crate) struct Vector2<T: Float + FromPrimitive + AsPrimitive<f64>> {
    pub x: T,
    pub y: T,
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Add for Vector2<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> AddAssign for Vector2<T> {
    fn add_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Sub for Vector2<T> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> SubAssign for Vector2<T> {
    fn sub_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Mul<T> for Vector2<T> {
    type Output = Self;

    fn mul(self, other: T) -> Self {
        Self {
            x: self.x * other,
            y: self.y * other,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> MulAssign<T> for Vector2<T> {
    fn mul_assign(&mut self, other: T) {
        *self = Self {
            x: self.x * other,
            y: self.y * other,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Div<T> for Vector2<T> {
    type Output = Self;

    fn div(self, other: T) -> Self {
        Self {
            x: self.x / other,
            y: self.y / other,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> DivAssign<T> for Vector2<T> {
    fn div_assign(&mut self, rhs: T) {
        *self = Self {
            x: self.x / rhs,
            y: self.y / rhs,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Neg for Vector2<T> {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
        }
    }
}

#[allow(dead_code)]
impl<T: Float + FromPrimitive + AsPrimitive<f64>> Vector2<T> {
    pub fn length_squared(&self) -> T {
        self.x * self.x + self.y * self.y
    }

    pub fn length(&self) -> T {
        Float::sqrt(self.length_squared())
    }

    pub fn abs(self) -> Self {
        Self {
            x: Float::abs(self.x),
            y: Float::abs(self.y),
        }
    }

    pub fn dot(&self, other: Self) -> T {
        self.x * other.x + self.y * other.y
    }

    pub fn abs_dot(&self, other: Self) -> T {
        Float::abs(self.dot(other))
    }

    pub fn normalize(self) -> Self {
        self / self.length()
    }
}
