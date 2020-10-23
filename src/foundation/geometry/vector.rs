use core::fmt;
use std::cmp::Ordering;
use std::fmt::{Display, Formatter};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

pub trait Float:
    Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + PartialEq
    + PartialOrd
    + Copy
    + Display
{
    fn abs(self) -> Self;
    fn partial_eq(&self, other: &Self) -> bool;
    fn partial_cmp(&self, other: &Self) -> Option<Ordering>;
    fn sqrt(self) -> Self;
    fn is_nan(self) -> bool;
}

impl Float for f32 {
    fn abs(self) -> Self {
        self.abs()
    }
    fn partial_eq(&self, other: &Self) -> bool {
        self == other
    }
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        std::cmp::PartialOrd::partial_cmp(&self, &other)
    }
    fn sqrt(self) -> Self {
        self.sqrt()
    }
    fn is_nan(self) -> bool {
        self.is_nan()
    }
}

impl Float for f64 {
    fn abs(self) -> Self {
        self.abs()
    }
    fn partial_eq(&self, other: &Self) -> bool {
        self == other
    }
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        std::cmp::PartialOrd::partial_cmp(&self, &other)
    }
    fn sqrt(self) -> Self {
        self.sqrt()
    }
    fn is_nan(self) -> bool {
        self.is_nan()
    }
}

fn max<T: Float>(first: T, second: T) -> T {
    if first < second {
        second
    } else {
        first
    }
}

fn min<T: Float>(first: T, second: T) -> T {
    if first < second {
        first
    } else {
        second
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Vector3<T: Float> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T: Float> Add for Vector3<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<T: Float> AddAssign for Vector3<T> {
    fn add_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<T: Float> Sub for Vector3<T> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<T: Float> SubAssign for Vector3<T> {
    fn sub_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<T: Float> Mul<T> for Vector3<T> {
    type Output = Self;

    fn mul(self, other: T) -> Self {
        Self {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl<T: Float> MulAssign<T> for Vector3<T> {
    fn mul_assign(&mut self, other: T) {
        *self = Self {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl<T: Float> Div<T> for Vector3<T> {
    type Output = Self;

    fn div(self, other: T) -> Self {
        Self {
            x: self.x / other,
            y: self.y / other,
            z: self.z / other,
        }
    }
}

impl<T: Float> DivAssign<T> for Vector3<T> {
    fn div_assign(&mut self, rhs: T) {
        *self = Self {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

impl<T: Float> Neg for Vector3<T> {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl<T: Float> Display for Vector3<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "x: {}, y: {}, z: {}", self.x, self.y, self.z)
    }
}

#[allow(dead_code)]
impl<T: Float> Vector3<T> {
    pub fn new(x: T, y: T, z: T) -> Self {
        debug_assert!(x.is_nan());
        debug_assert!(y.is_nan());
        debug_assert!(z.is_nan());

        Self { x, y, z }
    }

    pub fn length_squared(self) -> T {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    pub fn length(self) -> T {
        self.length_squared().sqrt()
    }

    pub fn abs(self) -> Self {
        Self {
            x: self.x.abs(),
            y: self.y.abs(),
            z: self.z.abs(),
        }
    }

    pub fn dot(&self, other: Self) -> T {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn abs_dot(&self, other: Self) -> T {
        self.dot(other).abs()
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
}

impl Vector3<f32> {
    pub fn cross(&self, other: Self) -> Self {
        let v1x = self.x as f64;
        let v2x = other.x as f64;
        let v1y = self.y as f64;
        let v2y = other.y as f64;
        let v1z = self.z as f64;
        let v2z = other.z as f64;

        Vector3 {
            x: ((v1y * v2z) - (v1z * v2y)) as f32,
            y: ((v1z * v2x) - (v1x * v2z)) as f32,
            z: ((v1x * v2y) - (v1y * v2x)) as f32,
        }
    }
}

impl Vector3<f64> {
    pub fn cross(&self, other: Self) -> Self {
        Self {
            x: (self.y * other.z) - (self.z * other.y),
            y: (self.z * other.x) - (self.x * other.z),
            z: (self.x * other.y) - (self.y * other.x),
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub(crate) struct Vector2<T: Float> {
    pub x: T,
    pub y: T,
}

impl<T: Float> Add for Vector2<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl<T: Float> AddAssign for Vector2<T> {
    fn add_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl<T: Float> Sub for Vector2<T> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl<T: Float> SubAssign for Vector2<T> {
    fn sub_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl<T: Float> Mul<T> for Vector2<T> {
    type Output = Self;

    fn mul(self, other: T) -> Self {
        Self {
            x: self.x * other,
            y: self.y * other,
        }
    }
}

impl<T: Float> MulAssign<T> for Vector2<T> {
    fn mul_assign(&mut self, other: T) {
        *self = Self {
            x: self.x * other,
            y: self.y * other,
        }
    }
}

impl<T: Float> Div<T> for Vector2<T> {
    type Output = Self;

    fn div(self, other: T) -> Self {
        Self {
            x: self.x / other,
            y: self.y / other,
        }
    }
}

impl<T: Float> DivAssign<T> for Vector2<T> {
    fn div_assign(&mut self, rhs: T) {
        *self = Self {
            x: self.x / rhs,
            y: self.y / rhs,
        }
    }
}

impl<T: Float> Neg for Vector2<T> {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
        }
    }
}

impl<T: Float> Display for Vector2<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "x: {}, y: {}", self.x, self.y)
    }
}

#[allow(dead_code)]
impl<T: Float> Vector2<T> {
    pub fn length_squared(&self) -> T {
        self.x * self.x + self.y * self.y
    }

    pub fn length(&self) -> T {
        self.length_squared().sqrt()
    }

    pub fn abs(self) -> Self {
        Self {
            x: self.x.abs(),
            y: self.y.abs(),
        }
    }

    pub fn dot(&self, other: Self) -> T {
        self.x * other.x + self.y * other.y
    }

    pub fn abs_dot(&self, other: Self) -> T {
        self.dot(other).abs()
    }

    pub fn normalize(self) -> Self {
        self / self.length()
    }
}
