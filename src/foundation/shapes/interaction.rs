use std::fmt::Debug;

use num_traits::{AsPrimitive, Float, FromPrimitive};

use crate::foundation::geometry::normal::Normal3;
use crate::foundation::geometry::point::{Point2, Point3};
use crate::foundation::geometry::vector::Vector3;

pub trait Interaction<T: Float + FromPrimitive + AsPrimitive<f64> + Debug> {
    fn point(&self) -> Point3<T>;
    fn time(&self) -> T;
    fn error(&self) -> Vector3<T>;
    fn wo(&self) -> Vector3<T>;
    fn normal(&self) -> Normal3<T>;
    fn is_surface_interaction(&self) -> bool;
}
