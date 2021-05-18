use crate::foundation::geometry::normal::Normal3;
use crate::foundation::geometry::point::Point3;
use crate::foundation::geometry::vector::Vector3;
use crate::foundation::pbr::Float;

pub trait Interaction {
    fn point(&self) -> Point3;
    fn time(&self) -> Float;
    fn error(&self) -> Vector3;
    fn wo(&self) -> Vector3;
    fn normal(&self) -> Normal3;
    fn is_surface_interaction(&self) -> bool;
}
