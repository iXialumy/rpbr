use crate::foundation::efloat::EFloat;
use num_traits::{AsPrimitive, Float, FromPrimitive};
use std::mem::swap;

pub fn gamma<T: Float>(n: T) -> T {
    (n * T::epsilon()) / (T::one() - n * T::epsilon())
}

pub fn quadratic<T: Float + AsPrimitive<f64> + FromPrimitive>(a: T, b: T, c: T) -> Option<(T, T)> {
    let a = a.as_();
    let b = b.as_();
    let c = c.as_();

    // Find quadratic discriminant
    let discrim = b * b - 4.0 * a * c;
    if discrim < 0.0 {
        return None;
    }
    let root_discrim = discrim.sqrt();

    // Compute quadratic _t_ values
    let mut q = 0.0;
    if b < 0.0 {
        q = -0.5 * (b - root_discrim);
    } else {
        q = -0.5 * (b + root_discrim);
    }

    let mut t0 = q / a;
    let mut t1 = c / q;
    if t0 > t1 {
        swap(&mut t0, &mut t1);
    }

    Some((T::from_f64(t0)?, T::from_f64(t1)?))
}

pub(crate) fn quadratic_ef<T: Float + AsPrimitive<f64> + FromPrimitive>(
    a: EFloat<T>,
    b: EFloat<T>,
    c: EFloat<T>,
) -> Option<(EFloat<T>, EFloat<T>)> {
    let a64 = a.value().as_();
    let b64 = b.value().as_();
    let c64 = c.value().as_();

    // Find quadratic discriminant
    let discrim = b64 * b64 - 4.0 * a64 * c64;
    if discrim < 0.0 {
        return None;
    }
    let root_discrim = T::from_f64(discrim.sqrt()).unwrap();
    let float_root_discrim = EFloat::new(root_discrim, T::epsilon() * root_discrim);

    // Compute quadratic _t_ values
    let mut q = EFloat::zero();
    if b64 < 0.0 {
        q = EFloat::from(T::from_f64(-0.5).unwrap()) * (b - float_root_discrim);
    } else {
        q = EFloat::from(T::from_f64(-0.5).unwrap()) * (b + float_root_discrim);
    }

    let mut t0 = q / a;
    let mut t1 = c / q;
    if t0 > t1 {
        swap(&mut t0, &mut t1);
    }

    Some((t0, t1))
}
