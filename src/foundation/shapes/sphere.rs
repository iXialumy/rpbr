use std::f64::consts::PI;
use std::ops::AddAssign;

use num::clamp;
use num::traits::FloatConst;

use crate::foundation::efloat::EFloat;
use crate::foundation::geometry::bounds::Bounds3;
use crate::foundation::geometry::normal::Normal3;
use crate::foundation::geometry::point::{Point2, Point3};
use crate::foundation::geometry::ray::Ray;
use crate::foundation::geometry::vector::Vector3;
use crate::foundation::pbr::Float;
use crate::foundation::shapes::shape::{Intersection, Shape};
use crate::foundation::shapes::surface_interaction::SurfaceInteraction;
use crate::foundation::transform::transform::Transform;
use crate::foundation::util::{gamma, quadratic_ef};

#[derive(Copy, Clone)]
struct Sphere {
    radius: Float,
    z_min: Float,
    z_max: Float,
    theta_min: Float,
    theta_max: Float,
    phi_max: Float,
    reverse_orientation: bool,
    world_to_object: Transform,
    object_to_world: Transform,
}
impl Sphere {
    fn new(
        object_to_world: Transform,
        world_to_object: Transform,
        reverse_orientation: bool,
        radius: Float,
        z_min: Float,
        z_max: Float,
        phi_max: Float,
    ) -> Self {
        Sphere {
            radius,
            z_min: clamp(Float::min(z_min, z_max), -radius, radius),
            z_max: clamp(Float::max(z_min, z_max), -radius, radius),
            theta_min: Float::acos(clamp(z_min / radius, -1.0, 1.0)),
            theta_max: Float::acos(clamp(z_max / radius, -1.0, 1.0)),
            phi_max: Float::to_radians(clamp(phi_max, 0.0, 360.0)),
            reverse_orientation,
            world_to_object,
            object_to_world,
        }
    }
}

impl Shape for Sphere {
    fn object_to_world(self) -> Transform {
        self.object_to_world
    }

    fn object_bounds(self) -> Bounds3 {
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

    fn area(&self) -> Float {
        self.phi_max * self.radius * (self.z_max - self.z_min)
    }

    fn intersect(&self, ray: Ray, _test_alpha_texture: bool) -> Option<Intersection> {
        let (t_shape_hit, p_hit, phi) = self.intersect_common(ray)?;

        // Find parametric representation of sphere hit
        let u = phi / self.phi_max;
        let theta = Float::acos(clamp(p_hit.z / self.radius, -1.0, 1.0));
        let v = (theta - self.theta_min) / (self.theta_max - self.theta_min);

        // Compute sphere dpdu and dpdv
        let z_radius = Float::sqrt(p_hit.x * p_hit.x + p_hit.y * p_hit.y);
        let inv_z_radius = 1.0 / z_radius;
        let cos_phi = p_hit.x * inv_z_radius;
        let sin_phi = p_hit.y * inv_z_radius;
        let dpdu = Vector3 {
            x: -self.phi_max * p_hit.y,
            y: self.phi_max * p_hit.x,
            z: 0.0,
        };
        let dpdv = Vector3 {
            x: p_hit.z * cos_phi,
            y: 0.0,
            z: 0.0,
        } * (self.theta_max - self.theta_min);

        // Compute sphere dndu and dndv
        let d2pduu = Vector3 {
            x: p_hit.x,
            y: p_hit.y,
            z: 0.0,
        } * -self.phi_max
            * self.phi_max;
        let d2pduv = Vector3 {
            x: -sin_phi,
            y: cos_phi,
            z: 0.0,
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
        let inv_egf2 = 1.0 / (E * G - F * F);
        let dndu = Normal3::from(
            (dpdu * (f * F - e * G) * inv_egf2) + (dpdv * (e * F - f * E) * inv_egf2),
        );
        let dndv = Normal3::from(
            (dpdu * (g * F - f * G) * inv_egf2) + (dpdv * (f * F - g * E) * inv_egf2),
        );

        // Compute error bounds for sphere intersection
        let p_error: Vector3 = (p_hit.abs() * gamma(5.0)).into();

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

    fn intersect_p(&self, ray: Ray, _test_alpha_texture: bool) -> bool {
        self.intersect_common(ray).is_some()
    }
}

impl Sphere {
    fn intersect_common(&self, ray: Ray) -> Option<(EFloat, Point3, f32)> {
        let (ray_obj, o_err, d_err) = self.world_to_object.transform_ray_with_error(ray);

        // Initialize _EFloat_ ray coordinate values
        let ox = EFloat::new(ray_obj.origin.x, o_err.x);
        let oy = EFloat::new(ray_obj.origin.y, o_err.y);
        let oz = EFloat::new(ray_obj.origin.z, o_err.z);
        let dx = EFloat::new(ray_obj.direction.x, d_err.x);
        let dy = EFloat::new(ray_obj.direction.y, d_err.y);
        let dz = EFloat::new(ray_obj.direction.z, d_err.z);

        let a = dx * dx + dy * dy + dz * dz;
        let b = (dx * ox + dy * oy + dz * oz) * 2.0;
        let c = ox * ox + oy * oy + oz * oz - EFloat::from(self.radius) * EFloat::from(self.radius);

        // Solve quadratic equation for _t_ values
        let option = quadratic_ef(a, b, c);
        if option == None {
            return None;
        }
        let (t0, t1) = option.unwrap();

        // Check quadratic shape _t0_ and _t1_ for nearest intersection
        if t0.upper_bound() > ray_obj.max_length || t1.lower_bound() <= 0.0 {
            return None;
        }
        let mut t_shape_hit = t0;
        if t_shape_hit.lower_bound() < 0.0 {
            t_shape_hit = t1;
            if t_shape_hit.upper_bound() > ray_obj.max_length {
                return None;
            }
        }

        // Compute sphere hit position and phi
        let mut p_hit = ray_obj.at(t_shape_hit.value());

        // Refine sphere intersection point
        p_hit *= self.radius / p_hit.distance(Point3::zero());
        if p_hit.x == 0.0 && p_hit.y == 0.0 {
            p_hit.x = 1e-5 * self.radius;
        }
        let mut phi = Float::atan2(p_hit.y, p_hit.x);
        if phi < 0.0 {
            phi += 2.0 * Float::PI()
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
            if p_hit.x == 0.0 && p_hit.y == 0.0 {
                p_hit.x = 1e-5 * self.radius;
            }
            let mut phi = Float::atan2(p_hit.y, p_hit.x);
            if phi < 0.0 {
                phi += 2.0 * Float::PI()
            }
            if (self.z_min > -self.radius && p_hit.z < self.z_min)
                || (self.z_max < self.radius && p_hit.z > self.z_max)
                || phi > self.phi_max
            {
                return None;
            }
        }
        Some((t_shape_hit, p_hit, phi))
    }
}
