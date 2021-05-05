use std::cmp::Ordering;
use std::ops::{Add, Div, Mul, Sub};

use num_traits::{AsPrimitive, Float, FromPrimitive};

#[derive(Copy, Clone)]
pub(crate) struct EFloat<T: Float> {
    value: T,
    low: T,
    high: T,
}

impl<T: Float> EFloat<T> {
    pub(crate) fn new(value: T, error: T) -> EFloat<T> {
        EFloat {
            value,
            low: value - error,
            high: value + error,
        }
    }
    pub(crate) fn default() -> Self {
        EFloat {
            value: T::zero(),
            low: T::zero(),
            high: T::zero(),
        }
    }
}

impl<T: Float> From<T> for EFloat<T> {
    fn from(value: T) -> Self {
        EFloat {
            value,
            ..EFloat::default()
        }
    }
}

impl<T: Float> Add for EFloat<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let value = self.value + rhs.value;
        let low = next_float_down(self.low + rhs.low);
        let high = next_float_up(self.high + rhs.high);

        // TODO: Check
        Self { value, low, high }
    }
}

impl<T: Float> Sub for EFloat<T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let value = self.value - rhs.value;
        let low = next_float_down(self.low - rhs.low);
        let high = next_float_up(self.high - rhs.high);

        // TODO: Check
        Self { value, low, high }
    }
}

impl<T: Float> Mul for EFloat<T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let value = self.value * rhs.value;

        let prod = vec![
            self.low * rhs.low,
            self.high * rhs.low,
            self.low * rhs.high,
            self.high * rhs.high,
        ];

        let low = next_float_down(T::min(T::min(prod[0], prod[1]), T::min(prod[2], prod[3])));
        let high = next_float_up(T::max(T::max(prod[0], prod[1]), T::max(prod[2], prod[3])));

        // TODO: Check
        Self { value, low, high }
    }
}

impl<T: Float> Div for EFloat<T> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let value = self.value / rhs.value;

        // TODO Precise
        let (low, high) = if rhs.low < T::zero() && rhs.high > T::zero() {
            (T::neg_infinity(), T::infinity())
        } else {
            let prod = vec![
                self.low / rhs.low,
                self.high / rhs.low,
                self.low / rhs.high,
                self.high / rhs.high,
            ];

            let l = next_float_down(T::min(T::min(prod[0], prod[1]), T::min(prod[2], prod[3])));
            let h = next_float_up(T::max(T::max(prod[0], prod[1]), T::max(prod[2], prod[3])));
            (l, h)
        };

        EFloat { value, low, high }
    }
}

impl<T: Float> PartialEq for EFloat<T> {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

impl<T: Float> PartialOrd for EFloat<T> {
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

fn next_float_up<T: Float>(f: T) -> T {
    unimplemented!()
}

fn next_float_down<T: Float>(f: T) -> T {
    unimplemented!()
}

fn quadratic<T: Float + FromPrimitive + AsPrimitive<f64>>(
    a: EFloat<T>,
    b: EFloat<T>,
    c: EFloat<T>,
) -> Option<(EFloat<T>, EFloat<T>)> {
    let av = a.value.as_();
    let bv = b.value.as_();
    let cv = c.value.as_();

    let discrim = bv * bv - 4.0 * av * cv;

    if discrim < 0.0 {
        return None;
    }

    let root_discrim = discrim.sqrt();
    let float_root_discrim = EFloat::new(
        T::from(root_discrim).unwrap(),
        T::from(f64::EPSILON * root_discrim).unwrap(),
    );

    let q = EFloat::from(T::from_f64(-0.5).unwrap())
        * (if b.value < T::zero() {
            b - float_root_discrim
        } else {
            b + float_root_discrim
        });

    let t0 = q / a;
    let t1 = c / q;
    if t0 < t1 {
        Some((t0, t1))
    } else {
        Some((t1, t0))
    }
}
