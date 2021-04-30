use crate::foundation::geometry::normal::Normal3;
use crate::foundation::geometry::point::Point3;
use crate::foundation::geometry::vector::Vector3;
use crate::foundation::shapes::interaction::Interaction;
use num_traits::{AsPrimitive, Float, FromPrimitive};
use std::fmt::Debug;

pub struct SurfaceInteraction<T: Float + FromPrimitive + AsPrimitive<f64> + Debug> {
    point: Point3<T>,
    time: T,
    error: Vector3<T>,
    wo: Vector3<T>,
    normal: Normal3<T>,
}

impl<T: Float + FromPrimitive + AsPrimitive<f64> + Debug> Interaction<T> for SurfaceInteraction<T> {
    fn point(&self) -> Point3<T> {
        self.point
    }

    fn time(&self) -> T {
        self.time
    }

    fn error(&self) -> Vector3<T> {
        self.error
    }

    fn wo(&self) -> Vector3<T> {
        self.wo
    }

    fn normal(&self) -> Normal3<T> {
        self.normal
    }

    fn is_surface_interaction(&self) -> bool {
        self.normal != Normal3::default()
    }
}
