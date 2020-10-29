use std::cmp::Ordering;
use std::fmt::Display;
use std::ops::{Add, Div, Mul, Neg, Sub};

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
    + Num
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

pub trait Num:
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
    fn ceil(self) -> Self;
    fn floor(self) -> Self;
    fn sqrt(self) -> f64;
}

pub fn max<T: Num>(first: T, second: T) -> T {
    if first < second {
        second
    } else {
        first
    }
}

pub fn min<T: Num>(first: T, second: T) -> T {
    if first < second {
        first
    } else {
        second
    }
}

impl Num for f32 {
    fn abs(self) -> Self {
        self.abs()
    }

    fn ceil(self) -> Self {
        self.ceil()
    }

    fn floor(self) -> Self {
        self.floor()
    }

    fn sqrt(self) -> f64 {
        self.sqrt() as f64
    }
}
impl Num for f64 {
    fn abs(self) -> Self {
        self.abs()
    }

    fn ceil(self) -> Self {
        self.ceil()
    }

    fn floor(self) -> Self {
        self.floor()
    }

    fn sqrt(self) -> f64 {
        self.sqrt()
    }
}

impl Num for i8 {
    fn abs(self) -> Self {
        self.abs()
    }

    fn ceil(self) -> Self {
        self
    }

    fn floor(self) -> Self {
        self
    }

    fn sqrt(self) -> f64 {
        (self as f64).sqrt()
    }
}
impl Num for i16 {
    fn abs(self) -> Self {
        self.abs()
    }

    fn ceil(self) -> Self {
        self
    }

    fn floor(self) -> Self {
        self
    }

    fn sqrt(self) -> f64 {
        (self as f64).sqrt()
    }
}
impl Num for i32 {
    fn abs(self) -> Self {
        self.abs()
    }

    fn ceil(self) -> Self {
        self
    }

    fn floor(self) -> Self {
        self
    }

    fn sqrt(self) -> f64 {
        (self as f64).sqrt()
    }
}
impl Num for i64 {
    fn abs(self) -> Self {
        self.abs()
    }

    fn ceil(self) -> Self {
        self
    }

    fn floor(self) -> Self {
        self
    }

    fn sqrt(self) -> f64 {
        (self as f64).sqrt()
    }
}
impl Num for i128 {
    fn abs(self) -> Self {
        self.abs()
    }

    fn ceil(self) -> Self {
        self
    }

    fn floor(self) -> Self {
        self
    }

    fn sqrt(self) -> f64 {
        (self as f64).sqrt()
    }
}
