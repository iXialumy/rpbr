use crate::foundation::geometry::point::Point3;
use crate::foundation::geometry::vector::Vector3;
use crate::foundation::pbr::Float;

#[derive(Copy, Clone)]
pub struct Ray {
    pub origin: Point3,
    pub direction: Vector3,
    pub max_length: Float,
    pub time: Float,
    //medium: Medium,
}

impl Ray {
    pub fn at(self, length: Float) -> Point3 {
        self.origin + (self.direction * length)
    }

    pub fn default() -> Ray {
        Ray {
            origin: Point3::zero(),
            direction: Vector3::empty(),
            max_length: Float::INFINITY,
            time: 0.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::foundation::geometry::point::Point3;
    use crate::foundation::geometry::ray::Ray;
    use crate::foundation::geometry::vector::Vector3;

    #[test]
    fn stretch_ray() {
        let r = Ray {
            origin: Point3 {
                x: 1.0,
                y: 1.0,
                z: 1.0,
            },
            direction: Vector3 {
                x: 0.0,
                y: 1.0,
                z: 0.0,
            },
            max_length: f32::INFINITY,
            time: 0.0,
        };
        let expected = Point3 {
            x: 1.0,
            y: 3.0,
            z: 1.0,
        };
        let actual = r.at(2.0);
        assert_eq!(expected, actual);
    }
}
