use std::fmt::Debug;

use num_traits::{AsPrimitive, Float, FromPrimitive};

use crate::foundation::geometry::bounds::Bounds3;
use crate::foundation::geometry::ray::Ray;
use crate::foundation::shapes::surface_interaction::SurfaceInteraction;
use crate::foundation::transform::transform::Transform;

pub struct Intersection<T: Float + FromPrimitive + AsPrimitive<f64> + Copy + Debug> {
    interaction: SurfaceInteraction<T>,
    distance: T,
}

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
    fn intersect(&self, ray: Ray<T>) -> Option<Intersection<T>>;
    fn intersect_p(&self, ray: Ray<T>, test_alpha_texture: bool) -> bool;
}
