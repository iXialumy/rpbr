use crate::foundation::geometry::bounds::Bounds3;
use crate::foundation::geometry::ray::Ray;
use crate::foundation::pbr::Float;
use crate::foundation::shapes::surface_interaction::SurfaceInteraction;
use crate::foundation::transform::transform::Transform;

pub struct Intersection {
    pub interaction: SurfaceInteraction,
    pub distance: Float,
}

pub trait Shape: Copy {
    fn object_to_world(self) -> Transform;
    fn object_bounds(self) -> Bounds3;
    fn world_bounds(self) -> Bounds3 {
        (self.object_to_world())(self.object_bounds())
    }
    fn area(&self) -> Float;
    fn intersect(&self, ray: Ray, test_alpha_texture: bool) -> Option<Intersection>;
    fn intersect_p(&self, ray: Ray, test_alpha_texture: bool) -> bool;
}
