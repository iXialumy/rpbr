use std::ops::{Add, AddAssign, Sub, SubAssign};

use crate::foundation::geometry::traits::{Num,min,max};

pub struct Point3<T: Num> {
    x: T,
    y: T,
    z: T,
}

impl<T: Num> Add for Point3<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<T: Num> AddAssign for Point3<T> {
    fn add_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<T: Num> Sub for Point3<T> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<T: Num> SubAssign for Point3<T> {
    fn sub_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<T: Num> Point3<T> {
    pub fn length_squared(self) -> T {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    pub fn length(self) -> f64 {
        self.length_squared().sqrt()
    }

    pub fn distance_squared(self, other: Self) -> T {
        (self - other).length_squared()
    }

    pub fn distance(self, other: Self) -> f64 {
        (self - other).length()
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

pub struct Point2<T: Num> {
    x: T,
    y: T,
}

impl<T: Num> Add for Point2<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl<T: Num> AddAssign for Point2<T> {
    fn add_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl<T: Num> Sub for Point2<T> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl<T: Num> SubAssign for Point2<T> {
    fn sub_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl<T: Num> Point2<T> {
    pub fn length_squared(self) -> T {
        self.x * self.x + self.y * self.y
    }

    pub fn length(self) -> f64 {
        self.length_squared().sqrt()
    }

    pub fn distance_squared(self, other: Self) -> T {
        (self - other).length_squared()
    }

    pub fn distance(self, other: Self) -> f64 {
        (self - other).length()
    }

    pub fn min(self, other: Self) -> Self {
        Self {
            x: min(self.x, other.x),
            y: min(self.y, other.y),
        }
    }

    pub fn max(self, other: Self) -> Self {
        Self {
            x: max(self.x, other.x),
            y: max(self.y, other.y),
        }
    }
}
