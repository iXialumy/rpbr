use std::cmp::{max, min};
use std::fmt::Debug;

use num_traits::{clamp, AsPrimitive, Float, FromPrimitive, Signed};

use crate::foundation::efloat::EFloat;
use crate::foundation::geometry::bounds::Bounds3;
use crate::foundation::geometry::normal::Normal3;
use crate::foundation::geometry::point::{Point2, Point3};
use crate::foundation::geometry::ray::Ray;
use crate::foundation::geometry::vector::Vector3;
use crate::foundation::shapes::shape::{Intersection, Shape};
use crate::foundation::shapes::surface_interaction::SurfaceInteraction;
use crate::foundation::transform::transform::Transform;
use crate::foundation::util::{gamma, quadratic, quadratic_ef};
use std::f64::consts::PI;
use std::ops::AddAssign;

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

impl<T: AsPrimitive<f64> + Copy + Float + FromPrimitive + Debug + Signed + AddAssign> Shape<T>
    for Sphere<T>
{
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
        let (ray_obj, o_err, d_err) = self.world_to_object.transform_ray_with_error(ray);

        // Initialize _EFloat_ ray coordinate values
        let ox = EFloat::new(ray_obj.origin.x, o_err.x);
        let oy = EFloat::new(ray_obj.origin.y, o_err.y);
        let oz = EFloat::new(ray_obj.origin.z, o_err.z);
        let dx = EFloat::new(ray_obj.direction.x, d_err.x);
        let dy = EFloat::new(ray_obj.direction.y, d_err.y);
        let dz = EFloat::new(ray_obj.direction.z, d_err.z);

        let a = dx * dx + dy * dy + dz * dz;
        let b = (dx * ox + dy * oy + dz * oz) * T::from_f64(2.0).unwrap();
        let c = ox * ox + oy * oy + oz * oz - EFloat::from(self.radius) * EFloat::from(self.radius);

        // Solve quadratic equation for _t_ values
        let option = quadratic_ef(a, b, c);
        if option == None {
            return None;
        }
        let (t0, t1) = option.unwrap();

        // Check quadratic shape _t0_ and _t1_ for nearest intersection
        if t0.upper_bound() > ray_obj.max_length || t1.lower_bound() <= T::zero() {
            return None;
        }
        let mut t_shape_hit = t0;
        if t_shape_hit.lower_bound() < T::zero() {
            t_shape_hit = t1;
            if t_shape_hit.upper_bound() > ray_obj.max_length {
                return None;
            }
        }

        // Compute sphere hit position and phi
        let mut p_hit = ray_obj.at(t_shape_hit.value());

        // Refine sphere intersection point
        p_hit *= self.radius / p_hit.distance(Point3::zero());
        if p_hit.x == T::zero() && p_hit.y == T::zero() {
            p_hit.x = T::from_f64(1e-5).unwrap() * self.radius;
        }
        let mut phi = T::atan2(p_hit.y, p_hit.x);
        if phi < T::zero() {
            phi += T::from_f64(2.0).unwrap() * T::from_f64(PI).unwrap()
        }

        // Test sphere intersection against clipping parameters
        if (self.z_min > -self.radius && p_hit.z < self.z_min)
            || (self.z_max < self.radius && p_hit.z > self.z_max)
            || phi > self.phi_max
        {
            if t_shape_hit == t1 {
                return None;
            }
            if t1.upper_bound() > ray_obj.max_length {
                return None;
            }
            t_shape_hit = t1;

            //Compute sphere hit position and phi
            p_hit = ray_obj.at(t_shape_hit.value());

            // Refine sphere intersection point
            p_hit *= self.radius / p_hit.distance(Point3::zero());
            if p_hit.x == T::zero() && p_hit.y == T::zero() {
                p_hit.x = T::from_f64(1e-5).unwrap() * self.radius;
            }
            let mut phi = T::atan2(p_hit.y, p_hit.x);
            if phi < T::zero() {
                phi += T::from_f64(2.0).unwrap() * T::from_f64(PI).unwrap()
            }
            if (self.z_min > -self.radius && p_hit.z < self.z_min)
                || (self.z_max < self.radius && p_hit.z > self.z_max)
                || phi > self.phi_max
            {
                return None;
            }
        }

        // Find parametric representation of sphere hit
        let u = phi / self.phi_max;
        let theta = T::acos(clamp(p_hit.z / self.radius, -T::one(), T::one()));
        let v = (theta - self.theta_min) / (self.theta_max - self.theta_min);

        // Compute sphere dpdu and dpdv
        let z_radius = T::sqrt(p_hit.x * p_hit.x + p_hit.y * p_hit.y);
        let inv_z_radius = T::one() / z_radius;
        let cos_phi = p_hit.x * inv_z_radius;
        let sin_phi = p_hit.y * inv_z_radius;
        let dpdu = Vector3 {
            x: -self.phi_max * p_hit.y,
            y: self.phi_max * p_hit.x,
            z: T::zero(),
        };
        let dpdv = Vector3 {
            x: p_hit.z * cos_phi,
            y: T::zero(),
            z: T::zero(),
        } * (self.theta_max - self.theta_min);

        // Compute sphere dndu and dndv
        let d2pduu = Vector3 {
            x: p_hit.x,
            y: p_hit.y,
            z: T::zero(),
        } * -self.phi_max
            * self.phi_max;
        let d2pduv = Vector3 {
            x: -sin_phi,
            y: cos_phi,
            z: T::zero(),
        } * (self.theta_max - self.theta_min)
            * p_hit.z
            * self.phi_max;
        let d2pdvv = Vector3 {
            x: p_hit.x,
            y: p_hit.y,
            z: p_hit.z,
        } * -(self.theta_max - self.theta_min)
            * (self.theta_max - self.theta_min);

        // CCompute coefficients for fundamental forms
        let E = dpdu.dot(dpdu);
        let F = dpdu.dot(dpdv);
        let G = dpdv.dot(dpdv);
        let N = dpdu.cross(dpdv).normalize();
        let e = N.dot(d2pduu);
        let f = N.dot(d2pduv);
        let g = N.dot(d2pdvv);

        // Compute dndu and dndv from fundamental form coefficients
        let inv_egf2 = T::one() / (E * G - F * F);
        let dndu = Normal3::from(
            (dpdu * (f * F - e * G) * inv_egf2) + (dpdv * (e * F - f * E) * inv_egf2),
        );
        let dndv = Normal3::from(
            (dpdu * (g * F - f * G) * inv_egf2) + (dpdv * (f * F - g * E) * inv_egf2),
        );

        // Compute error bounds for sphere intersection
        let p_error: Vector3<T> = (p_hit.abs() * gamma(T::from_f64(5.0).unwrap())).into();

        // Initialize SurfaceInteraction from parametric information
        let interact = SurfaceInteraction::new(
            p_hit,
            p_error,
            Point2::new(u, v),
            -ray.direction,
            dpdu,
            dpdv,
            dndu,
            dndv,
            ray.time,
        );

        let interaction = self.object_to_world.transform_surface_interaction(interact);

        Some(Intersection {
            interaction,
            distance: t_shape_hit.value(),
        })
    }

    fn intersect_p(&self, ray: Ray<T>, test_alpha_texture: bool) -> bool {
        todo!()
    }
}
