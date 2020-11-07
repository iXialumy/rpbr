use std::ops::{Add, AddAssign, Index, Sub, SubAssign};

use num_traits::{Float, FromPrimitive};

use crate::foundation::geometry::vector::Vector3;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Point3<T: Float> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T: Float> Add for Point3<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<T: Float> Add<Vector3<T>> for Point3<T> {
    type Output = Self;

    fn add(self, other: Vector3<T>) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

// impl<T: Float> AddAssign for Point3<T> {
//     fn add_assign(&mut self, other: Self) {
//         *self = Self {
//             x: self.x + other.x,
//             y: self.y + other.y,
//             z: self.z + other.z,
//         }
//     }
// }

impl<T: Float> Sub for Point3<T> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

// impl<T: Float> SubAssign for Point3<T> {
//     fn sub_assign(&mut self, other: Self) {
//         *self = Self {
//             x: self.x - other.x,
//             y: self.y - other.y,
//             z: self.z - other.z,
//         }
//     }
// }

impl<T: Float> Index<i32> for Point3<T> {
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

impl<T: Float> Point3<T> {
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

    // pub fn lerp(self) ->
}

impl Point3<f32> {
    pub fn new() -> Self {
        Point3 {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }
}

impl Point3<f64> {
    pub fn new() -> Self {
        Point3 {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }
}

pub struct Point2<T: Float> {
    x: T,
    y: T,
}

impl<T: Float> Add for Point2<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl<T: Float> AddAssign for Point2<T> {
    fn add_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl<T: Float> Sub for Point2<T> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl<T: Float> SubAssign for Point2<T> {
    fn sub_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl<T: Float> Point2<T> {
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

impl<T: Float> Index<i32> for Point2<T> {
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
