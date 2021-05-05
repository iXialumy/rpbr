use num_traits::{AsPrimitive, Float, FromPrimitive};
use std::mem::swap;

pub fn gamma<T: Float>(n: T) -> T {
    (n * T::epsilon()) / (T::one() - n * T::epsilon())
}

pub fn quadratic<T: Float + AsPrimitive<f64> + FromPrimitive>(a: T, b: T, c: T) -> Option<(T, T)> {
    let a = a.as_();
    let b = b.as_();
    let c = c.as_();

    let discrim = b * b - 4.0 * a * c;
    if discrim < 0.0 {
        return None;
    }
    let root_discrim = discrim.sqrt();
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
