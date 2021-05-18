use crate::foundation::efloat::EFloat;
use crate::foundation::pbr::Float;
use std::mem::swap;

pub fn gamma(n: Float) -> Float {
    (n * Float::EPSILON) / (1.0 - n * Float::EPSILON)
}

pub fn quadratic(a: Float, b: Float, c: Float) -> Option<(Float, Float)> {
    let a = a as f64;
    let b = b as f64;
    let c = c as f64;

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

    Some((t0 as Float, t1 as Float))
}

pub(crate) fn quadratic_ef(a: EFloat, b: EFloat, c: EFloat) -> Option<(EFloat, EFloat)> {
    let a64 = a.value() as f64;
    let b64 = b.value() as f64;
    let c64 = c.value() as f64;

    // Find quadratic discriminant
    let discrim = b64 * b64 - 4.0 * a64 * c64;
    if discrim < 0.0 {
        return None;
    }
    let root_discrim = discrim.sqrt() as Float;
    let float_root_discrim = EFloat::new(root_discrim, Float::EPSILON * root_discrim);

    // Compute quadratic _t_ values
    let mut q = EFloat::zero();
    if b64 < 0.0 {
        q = EFloat::from(-0.5) * (b - float_root_discrim);
    } else {
        q = EFloat::from(-0.5) * (b + float_root_discrim);
    }

    let mut t0 = q / a;
    let mut t1 = c / q;
    if t0 > t1 {
        swap(&mut t0, &mut t1);
    }

    Some((t0, t1))
}
