use std::fmt::Debug;

use crate::foundation::geometry::bounds::Bounds3;
use crate::foundation::geometry::point::Point3;
use crate::foundation::geometry::ray::Ray;
use crate::foundation::shapes::shape::{Intersection, Shape};
use crate::foundation::transform::transform::Transform;
use num_traits::{clamp, AsPrimitive, Float, FromPrimitive};
use std::cmp::{max, min};

#[derive(Copy, Clone)]
struct Sphere<T: Float + FromPrimitive + AsPrimitive<f64> + Copy + Debug> {
    radius: T,
    z_min: T,
    z_max: T,
    theta_min: T,
    theta_max: T,
    phi_max: T,
    reverse_orientation: bool,
    world_to_object: Transform<T>,
    object_to_world: Transform<T>,
}
impl<T: AsPrimitive<f64> + Copy + Float + FromPrimitive + Debug + Ord> Sphere<T> {
    fn new(
        object_to_world: Transform<T>,
        world_to_object: Transform<T>,
        reverse_orientation: bool,
        radius: T,
        z_min: T,
        z_max: T,
        phi_max: T,
    ) -> Self {
        Sphere {
            radius,
            z_min: clamp(min(z_min, z_max), -radius, radius),
            z_max: clamp(max(z_min, z_max), -radius, radius),
            theta_min: T::acos(clamp(z_min / radius, -T::one(), T::one())),
            theta_max: T::acos(clamp(z_max / radius, -T::one(), T::one())),
            phi_max: T::to_radians(clamp(phi_max, T::zero(), T::from_f64(360.0).unwrap())),
            reverse_orientation,
            world_to_object,
            object_to_world,
        }
    }
}

impl<T: AsPrimitive<f64> + Copy + Float + FromPrimitive + Debug> Shape<T> for Sphere<T> {
    fn object_to_world(self) -> Transform<T> {
        self.object_to_world
    }

    fn object_bounds(self) -> Bounds3<T> {
        Bounds3 {
            p_min: Point3 {
                x: -self.radius,
                y: -self.radius,
                z: self.z_min,
            },
            p_max: Point3 {
                x: self.radius,
                y: self.radius,
                z: self.z_max,
            },
        }
    }

    fn area(&self) -> T {
        todo!()
    }

    fn intersect(&self, ray: Ray<T>) -> Option<Intersection<T>> {
        todo!()
    }

    fn intersect_p(&self, ray: Ray<T>, test_alpha_texture: bool) -> bool {
        todo!()
    }
}
