use std::cmp::Ordering;

use num_traits::{AsPrimitive, Float, FromPrimitive};

use crate::foundation::geometry::point::Point3;
use crate::foundation::geometry::vector::Vector3;
use crate::foundation::transform::matrix::Matrix4x4;

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Transform<T: Float + FromPrimitive> {
    pub matrix: Matrix4x4<T>,
    pub inverse: Matrix4x4<T>,
}

impl<T: Float + FromPrimitive + AsPrimitive<f64> + Copy> Transform<T> {
    /// Create a transform with the identity matrix
    pub fn identity() -> Self {
        Self {
            matrix: Matrix4x4::identity(),
            inverse: Matrix4x4::identity(),
        }
    }

    pub fn new(matrix: Matrix4x4<T>) -> Self {
        Self {
            matrix,
            inverse: matrix.inverse(),
        }
    }

    pub fn inverse(transform: Self) -> Self {
        Self {
            matrix: transform.inverse,
            inverse: transform.matrix,
        }
    }

    pub fn transpose(transform: Self) -> Self {
        Self {
            matrix: transform.matrix.transpose(),
            inverse: transform.inverse.transpose(),
        }
    }

    pub fn is_identity(self) -> bool {
        self.matrix == Matrix4x4::<T>::identity()
    }

    pub fn translate(delta: Vector3<T>) -> Self {
        let matrix = Matrix4x4::matrix_from_floats(
            T::from_f64(1.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            delta.x,
            T::from_f64(0.0).unwrap(),
            T::from_f64(1.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            delta.y,
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(1.0).unwrap(),
            delta.z,
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(1.0).unwrap(),
        );

        let inverse = Matrix4x4::matrix_from_floats(
            T::from_f64(1.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            -delta.x,
            T::from_f64(0.0).unwrap(),
            T::from_f64(1.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            -delta.y,
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(1.0).unwrap(),
            -delta.z,
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(1.0).unwrap(),
        );

        Transform { matrix, inverse }
    }

    pub fn scale(x: T, y: T, z: T) -> Self {
        let matrix = Matrix4x4::matrix_from_floats(
            x,
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            y,
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            z,
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(1.0).unwrap(),
        );
        let inverse = Matrix4x4::matrix_from_floats(
            T::from_f64(1.0).unwrap() / x,
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(1.0).unwrap() / y,
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(1.0).unwrap() / z,
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(1.0).unwrap(),
        );
        Self { matrix, inverse }
    }

    pub fn rotate_x(theta: T) -> Self {
        let sin_theta = theta.to_radians().sin();
        let cos_theta = theta.to_radians().cos();

        let matrix = Matrix4x4::matrix_from_floats(
            T::from_f64(1.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            cos_theta,
            -sin_theta,
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            sin_theta,
            cos_theta,
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(1.0).unwrap(),
        );
        Self {
            matrix,
            inverse: matrix.transpose(),
        }
    }

    pub fn rotate_y(theta: T) -> Self {
        let sin_theta = theta.to_radians().sin();
        let cos_theta = theta.to_radians().cos();

        let matrix = Matrix4x4::matrix_from_floats(
            cos_theta,
            T::from_f64(0.0).unwrap(),
            sin_theta,
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(1.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            -sin_theta,
            T::from_f64(0.0).unwrap(),
            cos_theta,
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(1.0).unwrap(),
        );
        Self {
            matrix,
            inverse: matrix.transpose(),
        }
    }

    pub fn rotate(theta: T, axis: Vector3<T>) -> Self {
        let a = axis.normalize();
        let sin_theta = theta.to_radians().sin();
        let cos_theta = theta.to_radians().cos();
        let mut m = [[T::from_f64(0.0).unwrap(); 4]; 4];
        // Compute rotation of first basis vector
        m[0][0] = a.x * a.x + (T::from_f64(1.0).unwrap() - a.x * a.x) * cos_theta;
        m[0][1] = a.x * a.y * (T::from_f64(1.0).unwrap() - cos_theta) - a.z * sin_theta;
        m[0][2] = a.x * a.z * (T::from_f64(1.0).unwrap() - cos_theta) + a.y * sin_theta;
        m[0][3] = T::from_f64(0.0).unwrap();

        // Compute rotations of second and third basis vectors
        m[1][0] = a.x * a.y * (T::from_f64(1.0).unwrap() - cos_theta) + a.z * sin_theta;
        m[1][1] = a.y * a.y + (T::from_f64(1.0).unwrap() - a.y * a.y) * cos_theta;
        m[1][2] = a.y * a.z * (T::from_f64(1.0).unwrap() - cos_theta) - a.x * sin_theta;
        m[1][3] = T::from_f64(0.0).unwrap();

        m[2][0] = a.x * a.z * (T::from_f64(1.0).unwrap() - cos_theta) - a.y * sin_theta;
        m[2][1] = a.y * a.z * (T::from_f64(1.0).unwrap() - cos_theta) + a.x * sin_theta;
        m[2][2] = a.z * a.z + (T::from_f64(1.0).unwrap() - a.z * a.z) * cos_theta;
        m[2][3] = T::from_f64(0.0).unwrap();

        let matrix = m.into();
        Transform {
            matrix,
            inverse: matrix.transpose(),
        }
    }

    pub fn look_at(pos: Point3<T>, look: Point3<T>, up: Vector3<T>) -> Transform<T> {
        let mut camera_to_world = Matrix4x4::identity();
        camera_to_world.m[0][3] = pos.x;
        camera_to_world.m[1][3] = pos.y;
        camera_to_world.m[2][3] = pos.z;

        let dir = (look - pos).normalize();
        let right = up.normalize().cross(dir).normalize();
        let new_up = dir.cross(right);

        camera_to_world.m[0][0] = right.x;
        camera_to_world.m[1][0] = right.y;
        camera_to_world.m[2][0] = right.z;
        camera_to_world.m[3][0] = T::from_f64(0.0).unwrap();
        camera_to_world.m[0][1] = new_up.x;
        camera_to_world.m[1][1] = new_up.y;
        camera_to_world.m[2][1] = new_up.z;
        camera_to_world.m[3][1] = T::from_f64(0.0).unwrap();
        camera_to_world.m[0][2] = dir.x;
        camera_to_world.m[1][2] = dir.y;
        camera_to_world.m[2][2] = dir.z;
        camera_to_world.m[3][2] = T::from_f64(0.0).unwrap();

        Transform {
            matrix: camera_to_world.inverse(),
            inverse: camera_to_world,
        }
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64>> From<[[T; 4]; 4]> for Transform<T> {
    fn from(array: [[T; 4]; 4]) -> Self {
        Transform::new(Matrix4x4::from(array))
    }
}

impl<T: Float + FromPrimitive> PartialOrd for Transform<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        for i in 0..4 {
            for j in 0..4 {
                if self.matrix.m[i][j] < other.matrix.m[i][j] {
                    return Some(Ordering::Less);
                }
                if self.matrix.m[i][j] > other.matrix.m[i][j] {
                    return Some(Ordering::Greater);
                }
            }
        }
        Some(Ordering::Equal)
    }
}

#[cfg(test)]
mod tests {
    use crate::foundation::transform::matrix::Matrix4x4;
    use crate::foundation::transform::transform::Transform;

    #[test]
    pub fn new() {
        let matrix = Matrix4x4 {
            m: [
                [2.0, -1.0, 0.0, 0.0],
                [1.0, 2.0, -2.0, 0.0],
                [0.0, -1.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        };

        let inverse = Matrix4x4 {
            m: [
                [0.0, 1.0, 2.0, 0.0],
                [-1.0, 2.0, 4.0, 0.0],
                [-1.0, 2.0, 5.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        };

        let actual = Transform::new(matrix);

        let expected = Transform { matrix, inverse };

        assert_eq!(expected, actual);
    }
}
