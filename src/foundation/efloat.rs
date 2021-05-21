use std::cmp::Ordering;
use std::ops::{Add, Div, Mul, Sub};

use crate::foundation::pbr::Float;

#[derive(Copy, Clone)]
pub(crate) struct EFloat {
    value: Float,
    low: Float,
    high: Float,
}

impl EFloat {
    pub(crate) fn new(value: Float, error: Float) -> EFloat {
        EFloat {
            value,
            low: value - error,
            high: value + error,
        }
    }
    pub(crate) fn zero() -> Self {
        EFloat {
            value: 0.0,
            low: 0.0,
            high: 0.0,
        }
    }
    pub(crate) fn value(self) -> Float {
        self.value
    }
    pub(crate) fn upper_bound(self) -> Float {
        self.high
    }
    pub(crate) fn lower_bound(self) -> Float {
        self.low
    }
}

impl From<Float> for EFloat {
    fn from(value: Float) -> Self {
        EFloat {
            value,
            ..EFloat::zero()
        }
    }
}

impl Add for EFloat {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let value = self.value + rhs.value;
        let low = next_float_down(self.low + rhs.low);
        let high = next_float_up(self.high + rhs.high);

        // TODO: Check
        Self { value, low, high }
    }
}

impl Sub for EFloat {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let value = self.value - rhs.value;
        let low = next_float_down(self.low - rhs.low);
        let high = next_float_up(self.high - rhs.high);

        // TODO: Check
        Self { value, low, high }
    }
}

impl Mul for EFloat {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let value = self.value * rhs.value;

        let prod = vec![
            self.low * rhs.low,
            self.high * rhs.low,
            self.low * rhs.high,
            self.high * rhs.high,
        ];

        let low = next_float_down(Float::min(
            Float::min(prod[0], prod[1]),
            Float::min(prod[2], prod[3]),
        ));
        let high = next_float_up(Float::max(
            Float::max(prod[0], prod[1]),
            Float::max(prod[2], prod[3]),
        ));

        // TODO: Check
        Self { value, low, high }
    }
}

impl Mul<Float> for EFloat {
    type Output = Self;

    fn mul(self, rhs: Float) -> Self::Output {
        EFloat::from(rhs) * self
    }
}

impl Div for EFloat {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let value = self.value / rhs.value;

        // TODO Precise
        let (low, high) = if rhs.low < 0.0 && rhs.high > 0.0 {
            (Float::NEG_INFINITY, Float::INFINITY)
        } else {
            let prod = vec![
                self.low / rhs.low,
                self.high / rhs.low,
                self.low / rhs.high,
                self.high / rhs.high,
            ];

            let l = next_float_down(Float::min(
                Float::min(prod[0], prod[1]),
                Float::min(prod[2], prod[3]),
            ));
            let h = next_float_up(Float::max(
                Float::max(prod[0], prod[1]),
                Float::max(prod[2], prod[3]),
            ));
            (l, h)
        };

        EFloat { value, low, high }
    }
}

impl PartialEq for EFloat {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

impl PartialOrd for EFloat {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self == other {
            return Some(Ordering::Equal);
        }

        if self.value < other.value {
            Some(Ordering::Less)
        } else {
            Some(Ordering::Greater)
        }
    }
}

pub fn next_float_up(f: Float) -> Float {
    unimplemented!("EFloat::next_float_up({})", f)
}

fn next_float_down(f: Float) -> Float {
    unimplemented!("EFloat::next_float_down({})", f)
}
