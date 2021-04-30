use std::fmt::Debug;

use num_traits::{AsPrimitive, Float, FromPrimitive};

use crate::foundation::geometry::bounds::Bounds3;
use crate::foundation::transform::transform::Transform;

pub trait Shape<T: Float + FromPrimitive + AsPrimitive<f64> + Copy + Debug>: Copy {
    fn new(
        object_to_world: &Transform<T>,
        world_to_object: &Transform<T>,
        reverse_orientation: bool,
    ) -> Self;
    fn object_to_world<'a>(self) -> &'a Transform<T>;
    fn object_bounds(self) -> Bounds3<T>;
    fn world_bounds(self) -> Bounds3<T> {
        (*self.object_to_world())(self.object_bounds())
    }
    fn area(&self) -> T;
}
