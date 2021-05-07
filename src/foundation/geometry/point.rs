use std::ops::{Add, AddAssign, Div, Index, Mul, MulAssign, Sub, SubAssign};

use num_traits::{AsPrimitive, Float, FromPrimitive};

use crate::foundation::geometry::vector::Vector3;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Point3<T: Float + FromPrimitive + AsPrimitive<f64>> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Add for Point3<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Add<Vector3<T>> for Point3<T> {
    type Output = Self;

    fn add(self, other: Vector3<T>) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Sub for Point3<T> {
    type Output = Vector3<T>;

    fn sub(self, other: Self) -> Self::Output {
        Vector3 {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Sub<Vector3<T>> for Point3<T> {
    type Output = Point3<T>;

    fn sub(self, other: Vector3<T>) -> Self::Output {
        Point3 {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Index<i32> for Point3<T> {
    type Output = T;

    fn index(&self, index: i32) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Index out of bounds: {}", index),
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Mul<T> for Point3<T> {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> MulAssign<T> for Point3<T> {
    fn mul_assign(&mut self, rhs: T) {
        self.x = self.x * rhs;
        self.y = self.y * rhs;
        self.z = self.z * rhs;
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Div<T> for Point3<T> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        Self {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Point3<T> {
    pub fn zero() -> Point3<T> {
        Point3 {
            x: T::zero(),
            y: T::zero(),
            z: T::zero(),
        }
    }

    pub fn new(x: T, y: T, z: T) -> Self {
        Self { x, y, z }
    }

    pub fn length_squared(self) -> T {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    pub fn length(self) -> T {
        self.length_squared().sqrt()
    }

    pub fn distance_squared(self, other: Self) -> T {
        (self - other).length_squared()
    }

    pub fn distance(self, other: Self) -> T {
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

pub struct Point2<T: Float + FromPrimitive + AsPrimitive<f64>> {
    x: T,
    y: T,
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Add for Point2<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> AddAssign for Point2<T> {
    fn add_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Sub for Point2<T> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> SubAssign for Point2<T> {
    fn sub_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Point2<T> {
    pub fn length_squared(self) -> T {
        self.x * self.x + self.y * self.y
    }

    pub fn length(self) -> T {
        self.length_squared().sqrt()
    }

    pub fn distance_squared(self, other: Self) -> T {
        (self - other).length_squared()
    }

    pub fn distance(self, other: Self) -> T {
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

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Index<i32> for Point2<T> {
    type Output = T;

    fn index(&self, index: i32) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            _ => panic!("Index out of bounds: {}", index),
        }
    }
}

fn lerp<T: Float + FromPrimitive>(t: T, v1: T, v2: T) -> T {
    (T::from_f64(1.0).unwrap() - t) * v1 + t * v2
}
