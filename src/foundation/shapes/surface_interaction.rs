use crate::foundation::geometry::normal::Normal3;
use crate::foundation::geometry::point::{Point2, Point3};
use crate::foundation::geometry::vector::Vector3;
use crate::foundation::pbr::Float;
use crate::foundation::shapes::interaction::Interaction;

pub struct Shading {
    pub n: Normal3,
    pub dpdu: Vector3,
    pub dpdv: Vector3,
    pub dndu: Normal3,
    pub dndv: Normal3,
}

pub struct CommonInteraction {
    pub(crate) point: Point3,
    pub(crate) time: Float,
    pub(crate) error: Vector3,
    pub(crate) wo: Vector3,
    pub(crate) normal: Normal3,
}

pub struct SurfaceInteraction {
    pub(crate) interaction: CommonInteraction,
    pub(crate) uv: Point2,
    pub(crate) dpdu: Vector3,
    pub(crate) dpdv: Vector3,
    pub(crate) dndu: Normal3,
    pub(crate) dndv: Normal3,
    pub(crate) shading: Shading,
}

impl SurfaceInteraction {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        p: Point3,
        p_error: Vector3,
        uv: Point2,
        wo: Vector3,
        dpdu: Vector3,
        dpdv: Vector3,
        dndu: Normal3,
        dndv: Normal3,
        time: Float,
        // shape: Option<& dyn Shape>,
    ) -> Self {
        let cross: Normal3 = (dpdu.cross(dpdv)).into();
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

impl Interaction for SurfaceInteraction {
    fn point(&self) -> Point3 {
        self.interaction.point
    }

    fn time(&self) -> Float {
        self.interaction.time
    }

    fn error(&self) -> Vector3 {
        self.interaction.error
    }

    fn wo(&self) -> Vector3 {
        self.interaction.wo
    }

    fn normal(&self) -> Normal3 {
        self.interaction.normal
    }

    fn is_surface_interaction(&self) -> bool {
        self.normal() != Normal3::zero()
    }
}
