use crate::foundation::geometry::normal::Normal3;
use crate::foundation::geometry::point::{Point2, Point3};
use crate::foundation::geometry::vector::Vector3;
use crate::foundation::shapes::interaction::Interaction;
use num_traits::{AsPrimitive, Float, FromPrimitive};
use std::fmt::Debug;

pub struct Shading<T: Float + AsPrimitive<f64> + FromPrimitive> {
    n: Normal3<T>,
    dpdu: Vector3<T>,
    dpdv: Vector3<T>,
    dndu: Normal3<T>,
    dndv: Normal3<T>,
}

pub struct CommonInteraction<T: Float + FromPrimitive + AsPrimitive<f64>> {
    pub(crate) point: Point3<T>,
    pub(crate) time: T,
    pub(crate) error: Vector3<T>,
    pub(crate) wo: Vector3<T>,
    pub(crate) normal: Normal3<T>,
}

pub struct SurfaceInteraction<T: Float + FromPrimitive + AsPrimitive<f64> + Debug> {
    pub(crate) interaction: CommonInteraction<T>,
    pub(crate) uv: Point2<T>,
    pub(crate) dpdu: Vector3<T>,
    pub(crate) dpdv: Vector3<T>,
    pub(crate) dndu: Normal3<T>,
    pub(crate) dndv: Normal3<T>,
    pub(crate) shading: Shading<T>,
}

impl<T: Float + FromPrimitive + AsPrimitive<f64> + Debug> SurfaceInteraction<T> {
    pub fn new(
        p: Point3<T>,
        p_error: Vector3<T>,
        uv: Point2<T>,
        wo: Vector3<T>,
        dpdu: Vector3<T>,
        dpdv: Vector3<T>,
        dndu: Normal3<T>,
        dndv: Normal3<T>,
        time: T,
        // shape: Option<& dyn Shape<T>>,
    ) -> Self {
        let cross: Normal3<T> = (dpdu.cross(dpdv)).into();
        let n = cross.normalize();

        SurfaceInteraction {
            interaction: CommonInteraction {
                point: p,
                time,
                error: p_error,
                wo,
                normal: n,
            },
            uv,
            dpdu,
            dpdv,
            dndu,
            dndv,
            shading: Shading {
                n,
                dpdu,
                dpdv,
                dndu,
                dndv,
            },
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64> + Debug> Interaction<T> for SurfaceInteraction<T> {
    fn point(&self) -> Point3<T> {
        self.interaction.point
    }

    fn time(&self) -> T {
        self.interaction.time
    }

    fn error(&self) -> Vector3<T> {
        self.interaction.error
    }

    fn wo(&self) -> Vector3<T> {
        self.interaction.wo
    }

    fn normal(&self) -> Normal3<T> {
        self.interaction.normal
    }

    fn is_surface_interaction(&self) -> bool {
        self.normal() != Normal3::zero()
    }
}
