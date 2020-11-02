use crate::foundation::geometry::point::Point3;
use crate::foundation::geometry::vector::Vector3;
use num_traits::Float;

pub struct Ray<T: Float> {
    pub origin: Point3<T>,
    pub direction: Vector3<T>,
    pub max_length: T,
    pub time: T,
    //medium: Medium,
}

impl<T: Float> Ray<T> {
    pub fn at(self, length: T) -> Point3<T> {
        self.origin + (self.direction * length)
    }
}

impl Ray<f32> {
    pub fn new() -> Ray<f32> {
        Ray {
            origin: Point3::<f32>::new(),
            direction: Vector3::<f32>::empty(),
            max_length: f32::INFINITY,
            time: 0.0,
        }
    }
}

impl Ray<f64> {
    pub fn new() -> Ray<f64> {
        Ray {
            origin: Point3::<f64>::new(),
            direction: Vector3::<f64>::empty(),
            max_length: f64::INFINITY,
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
