use num_traits::{AsPrimitive, Float, FromPrimitive};

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Normal3<T: Float + FromPrimitive + AsPrimitive<f64> + Copy> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Normal3<T> {}
