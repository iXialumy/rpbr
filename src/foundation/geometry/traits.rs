use std::fmt::Display;
use std::ops::{Add, Div, Mul, Neg, Sub};

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