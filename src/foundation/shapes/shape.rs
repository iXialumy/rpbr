use crate::foundation::geometry::bounds::Bounds3;
use crate::foundation::geometry::ray::Ray;
use crate::foundation::pbr::Float;
use crate::foundation::shapes::sphere::Sphere;
use crate::foundation::shapes::surface_interaction::SurfaceInteraction;
use crate::foundation::transforms::transform::Transform;

pub struct Intersection {
    pub interaction: SurfaceInteraction,
    pub distance: Float,
}

pub enum Shape {
    Sphr(Sphere),
}

impl Shape {
    pub fn object_to_world(self) -> Transform {
        match self {
            Shape::Sphr(shape) => shape.object_to_world(),
        }
    }

    pub fn object_bounds(self) -> Bounds3 {
        match self {
            Shape::Sphr(shape) => shape.object_bounds(),
        }
    }

    pub fn world_bounds(self) -> Bounds3 {
        (self.object_to_world())(self.object_bounds())
    }

    pub fn area(&self) -> Float {
        match self {
            Shape::Sphr(shape) => shape.area(),
        }
    }

    pub fn intersect(&self, ray: Ray, test_alpha_texture: bool) -> Option<Intersection> {
        match self {
            Shape::Sphr(shape) => shape.intersect(ray, test_alpha_texture),
        }
    }

    pub fn intersect_p(&self, ray: Ray, test_alpha_texture: bool) -> bool {
        match self {
            Shape::Sphr(shape) => shape.intersect_p(ray, test_alpha_texture),
        }
    }
}
