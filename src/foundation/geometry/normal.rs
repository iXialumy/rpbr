use num_traits::{AsPrimitive, Float, FromPrimitive};

use crate::foundation::geometry::point::Point3;

type Normal3<T> = Point3<T>;

impl<T: Float + FromPrimitive + AsPrimitive<f64>> Normal3<T> {}
