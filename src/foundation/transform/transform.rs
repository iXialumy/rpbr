use std::cmp::Ordering;
use std::fmt::Debug;
use std::ops::Mul;

use crate::foundation::geometry::bounds::Bounds3;
use crate::foundation::geometry::normal::Normal3;
use crate::foundation::geometry::point::Point3;
use crate::foundation::geometry::ray::Ray;
use crate::foundation::geometry::vector::Vector3;
use crate::foundation::pbr::Float;
use crate::foundation::shapes::interaction::Interaction;
use crate::foundation::shapes::surface_interaction::{
    CommonInteraction, Shading, SurfaceInteraction,
};
use crate::foundation::transform::matrix::Matrix4x4;
use crate::foundation::util::gamma;
use num::abs;

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Transform {
    pub matrix: Matrix4x4,
    pub inverse: Matrix4x4,
}

impl Transform {
    /// Create a transform with the identity matrix
    pub fn identity() -> Self {
        Self {
            matrix: Matrix4x4::identity(),
            inverse: Matrix4x4::identity(),
        }
    }

    pub fn new(matrix: Matrix4x4) -> Self {
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
        self.matrix == Matrix4x4::identity()
    }

    pub fn translate(delta: Vector3) -> Self {
        let matrix = Matrix4x4::matrix_from_floats(
            1.0, 0.0, 0.0, delta.x, 0.0, 1.0, 0.0, delta.y, 0.0, 0.0, 1.0, delta.z, 0.0, 0.0, 0.0,
            1.0,
        );

        let inverse = Matrix4x4::matrix_from_floats(
            1.0, 0.0, 0.0, -delta.x, 0.0, 1.0, 0.0, -delta.y, 0.0, 0.0, 1.0, -delta.z, 0.0, 0.0,
            0.0, 1.0,
        );

        Transform { matrix, inverse }
    }

    pub fn scale(x: Float, y: Float, z: Float) -> Self {
        let matrix = Matrix4x4::matrix_from_floats(
            x, 0.0, 0.0, 0.0, 0.0, y, 0.0, 0.0, 0.0, 0.0, z, 0.0, 0.0, 0.0, 0.0, 1.0,
        );
        let inverse = Matrix4x4::matrix_from_floats(
            1.0 / x,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0 / y,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0 / z,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
        );
        Self { matrix, inverse }
    }

    pub fn rotate_x(theta: Float) -> Self {
        let sin_theta = theta.to_radians().sin();
        let cos_theta = theta.to_radians().cos();

        let matrix = Matrix4x4::matrix_from_floats(
            1.0, 0.0, 0.0, 0.0, 0.0, cos_theta, -sin_theta, 0.0, 0.0, sin_theta, cos_theta, 0.0,
            0.0, 0.0, 0.0, 1.0,
        );
        Self {
            matrix,
            inverse: matrix.transpose(),
        }
    }

    pub fn rotate_y(theta: Float) -> Self {
        let sin_theta = theta.to_radians().sin();
        let cos_theta = theta.to_radians().cos();

        let matrix = Matrix4x4::matrix_from_floats(
            cos_theta, 0.0, sin_theta, 0.0, 0.0, 1.0, 0.0, 0.0, -sin_theta, 0.0, cos_theta, 0.0,
            0.0, 0.0, 0.0, 1.0,
        );
        Self {
            matrix,
            inverse: matrix.transpose(),
        }
    }

    pub fn rotate(theta: Float, axis: Vector3) -> Self {
        let a = axis.normalize();
        let sin_theta = theta.to_radians().sin();
        let cos_theta = theta.to_radians().cos();
        let mut m = [[0.0; 4]; 4];
        // Compute rotation of first basis vector
        m[0][0] = a.x * a.x + (1.0 - a.x * a.x) * cos_theta;
        m[0][1] = a.x * a.y * (1.0 - cos_theta) - a.z * sin_theta;
        m[0][2] = a.x * a.z * (1.0 - cos_theta) + a.y * sin_theta;
        m[0][3] = 0.0;

        // Compute rotations of second and third basis vectors
        m[1][0] = a.x * a.y * (1.0 - cos_theta) + a.z * sin_theta;
        m[1][1] = a.y * a.y + (1.0 - a.y * a.y) * cos_theta;
        m[1][2] = a.y * a.z * (1.0 - cos_theta) - a.x * sin_theta;
        m[1][3] = 0.0;

        m[2][0] = a.x * a.z * (1.0 - cos_theta) - a.y * sin_theta;
        m[2][1] = a.y * a.z * (1.0 - cos_theta) + a.x * sin_theta;
        m[2][2] = a.z * a.z + (1.0 - a.z * a.z) * cos_theta;
        m[2][3] = 0.0;

        let matrix = m.into();
        Transform {
            matrix,
            inverse: matrix.transpose(),
        }
    }

    pub fn look_at(pos: Point3, look: Point3, up: Vector3) -> Transform {
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
        camera_to_world.m[3][0] = 0.0;
        camera_to_world.m[0][1] = new_up.x;
        camera_to_world.m[1][1] = new_up.y;
        camera_to_world.m[2][1] = new_up.z;
        camera_to_world.m[3][1] = 0.0;
        camera_to_world.m[0][2] = dir.x;
        camera_to_world.m[1][2] = dir.y;
        camera_to_world.m[2][2] = dir.z;
        camera_to_world.m[3][2] = 0.0;

        Transform {
            matrix: camera_to_world.inverse(),
            inverse: camera_to_world,
        }
    }

    pub fn transform_point(self, point: Point3) -> Point3 {
        let x = point.x;
        let y = point.y;
        let z = point.z;
        let m = self.matrix.m;

        let xp = m[0][0] * x + m[0][1] * y + m[0][2] * z + m[0][3];
        let yp = m[1][0] * x + m[1][1] * y + m[1][2] * z + m[1][3];
        let zp = m[2][0] * x + m[2][1] * y + m[2][2] * z + m[2][3];
        let wp = m[3][0] * x + m[3][1] * y + m[3][2] * z + m[3][3];

        assert_ne!(wp, 0.0);
        if wp == 0.0 {
            Point3 {
                x: xp,
                y: yp,
                z: zp,
            }
        } else {
            Point3 {
                x: xp,
                y: yp,
                z: zp,
            } / wp
        }
    }

    pub fn transform_point_with_error(self, point: Point3) -> (Point3, Vector3) {
        let x = point.x;
        let y = point.y;
        let z = point.z;
        let m = self.matrix.m;

        let xp = m[0][0] * x + m[0][1] * y + m[0][2] * z + m[0][3];
        let yp = m[1][0] * x + m[1][1] * y + m[1][2] * z + m[1][3];
        let zp = m[2][0] * x + m[2][1] * y + m[2][2] * z + m[2][3];
        let wp = m[3][0] * x + m[3][1] * y + m[3][2] * z + m[3][3];

        let x_abs_sum = <Transform>::calc_x_abs_sum(&x, &y, &z, m);
        let y_abs_sum = <Transform>::calc_y_abs_sum(&x, &y, &z, m);
        let z_abs_sum = <Transform>::calc_z_abs_sum(&x, &y, &z, m);

        let error: Vector3 = Vector3::new(x_abs_sum, y_abs_sum, z_abs_sum) * gamma(3.0);

        debug_assert_eq!(wp, 0.0);

        let p = Point3 {
            x: xp,
            y: yp,
            z: zp,
        };
        if wp == 1.0 {
            (p, error)
        } else {
            (p / wp, error)
        }
    }

    pub fn transform_point_with_errors(self, point: Point3, error: Vector3) -> (Point3, Vector3) {
        let x = point.x;
        let y = point.y;
        let z = point.z;
        let m = self.matrix.m;

        let xp = m[0][0] * x + m[0][1] * y + m[0][2] * z + m[0][3];
        let yp = m[1][0] * x + m[1][1] * y + m[1][2] * z + m[1][3];
        let zp = m[2][0] * x + m[2][1] * y + m[2][2] * z + m[2][3];
        let wp = m[3][0] * x + m[3][1] * y + m[3][2] * z + m[3][3];

        let abs_err_x = (gamma(3.0) + 1.0)
            * (<Transform>::calc_x_abs_sum(&error.x, &error.y, &error.z, m))
            + gamma(3.0) * <Transform>::calc_x_abs_sum(&x, &y, &z, m);
        let abs_err_y = (gamma(3.0) + 1.0)
            * (<Transform>::calc_y_abs_sum(&error.x, &error.y, &error.z, m))
            + gamma(3.0) * <Transform>::calc_y_abs_sum(&x, &y, &z, m);
        let abs_err_z = (gamma(3.0) + 1.0)
            * (<Transform>::calc_z_abs_sum(&error.x, &error.y, &error.z, m))
            + gamma(3.0) * <Transform>::calc_z_abs_sum(&x, &y, &z, m);

        let error: Vector3 = Vector3::new(abs_err_x, abs_err_y, abs_err_z);

        debug_assert_eq!(wp, 0.0);

        let p = Point3 {
            x: xp,
            y: yp,
            z: zp,
        };
        if wp == 1.0 {
            (p, error)
        } else {
            (p / wp, error)
        }
    }

    fn calc_x_abs_sum(x: &Float, y: &Float, z: &Float, m: [[Float; 4]; 4]) -> Float {
        Float::abs(m[0][0] * *x)
            + Float::abs(m[0][1] * *y)
            + Float::abs(m[0][2] * *z)
            + Float::abs(m[0][3])
    }

    fn calc_y_abs_sum(x: &Float, y: &Float, z: &Float, m: [[Float; 4]; 4]) -> Float {
        Float::abs(m[1][0] * *x)
            + Float::abs(m[1][1] * *y)
            + Float::abs(m[1][2] * *z)
            + Float::abs(m[1][3])
    }

    fn calc_z_abs_sum(x: &Float, y: &Float, z: &Float, m: [[Float; 4]; 4]) -> Float {
        Float::abs(m[2][0] * *x)
            + Float::abs(m[2][1] * *y)
            + Float::abs(m[2][2] * *z)
            + Float::abs(m[2][3])
    }

    pub fn transform_vector(self, vector: Vector3) -> Vector3 {
        self(vector)
    }

    pub fn transform_vector_with_error(self, vector: Vector3) -> (Vector3, Vector3) {
        let err_x = gamma(3.0)
            * (abs(self.matrix.m[0][0] * vector.x)
                + abs(self.matrix.m[0][1] * vector.y)
                + abs(self.matrix.m[0][2] * vector.z));

        let err_y = gamma(3.0)
            * (abs(self.matrix.m[1][0] * vector.x)
                + abs(self.matrix.m[1][1] * vector.y)
                + abs(self.matrix.m[1][2] * vector.z));

        let err_z = gamma(3.0)
            * (abs(self.matrix.m[2][0] * vector.x)
                + abs(self.matrix.m[2][1] * vector.y)
                + abs(self.matrix.m[2][2] * vector.z));

        let v_out = Vector3 {
            x: self.matrix.m[0][0] * vector.x
                + self.matrix.m[0][1] * vector.y
                + self.matrix.m[0][2] * vector.z,
            y: self.matrix.m[1][0] * vector.x
                + self.matrix.m[1][1] * vector.y
                + self.matrix.m[1][2] * vector.z,
            z: self.matrix.m[2][0] * vector.x
                + self.matrix.m[2][1] * vector.y
                + self.matrix.m[2][2] * vector.z,
        };

        let err_out = Vector3 {
            x: err_x,
            y: err_y,
            z: err_z,
        };

        (v_out, err_out)
    }

    pub fn transform_ray_with_error(self, ray: Ray) -> (Ray, Vector3, Vector3) {
        let (mut origin, origin_error) = self.transform_point_with_error(ray.origin);
        let (direction, direction_error) = self.transform_vector_with_error(ray.direction);

        let max_length = ray.max_length;
        let length_squared = direction.length_squared();

        if length_squared > 0.0 {
            let dt = direction.abs().dot(origin_error) / length_squared;
            origin = origin + direction * dt;
        }

        let r = Ray {
            origin,
            direction,
            max_length,
            time: ray.time,
        };

        (r, origin_error, direction_error)
    }

    pub fn transform_surface_interaction(self, si: SurfaceInteraction) -> SurfaceInteraction {
        let (point, error) = self.transform_point_with_errors(si.point(), si.error());

        SurfaceInteraction {
            interaction: CommonInteraction {
                point,
                time: si.time(),
                error,
                wo: self(si.wo()).normalize(),
                normal: self(si.normal()).normalize(),
            },
            uv: si.uv,
            dpdu: self(si.dpdu),
            dpdv: self(si.dpdv),
            dndu: self(si.dndu),
            dndv: self(si.dndv),
            shading: Shading {
                n: self(si.shading.n),
                dpdu: self(si.shading.dpdu),
                dpdv: self(si.shading.dpdv),
                dndu: self(si.shading.dndu),
                dndv: self(si.shading.dndv),
            },
        }
    }

    /// Tells if a transformatin will swap the handedness of a coordinate system
    pub fn swaps_handedness(self) -> bool {
        self.matrix.det() < 0.0
    }
}

impl From<[[Float; 4]; 4]> for Transform {
    fn from(array: [[Float; 4]; 4]) -> Self {
        Transform::new(Matrix4x4::from(array))
    }
}

impl PartialOrd for Transform {
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

/// Transform a Point with a matrix
///
/// Example:
///
/// ```
/// #![feature(fn_traits)]
/// use rpbr::foundation::geometry::point::Point3;
/// use rpbr::foundation::transform::matrix::Matrix4x4;
/// use rpbr::foundation::transform::transform::Transform;
/// let p = Point3 {
///     x: 1.0,
///     y: 2.0,
///     z: 3.0
/// };
///
///
/// let m = Matrix4x4::from([[-5.0,   2.0,   3.0,  0.0]
///                         ,[ 1.0,   1.0,  -1.0,  2.0]
///                         ,[ 2.0,  -2.0,  -2.0,  0.0]
///                         ,[ 0.0,  -2.0,   1.0, -3.0]
///                         ]
/// );
///
/// let t = Transform::new(m);
///
/// let expected = Point3 {
///     x: -2.0,
///     y: -0.5,
///     z:  2.0
/// };
///
/// assert_eq!(expected, t(p));
/// ```
impl Fn<(Point3,)> for Transform {
    extern "rust-call" fn call(&self, args: (Point3,)) -> Point3 {
        let point = args.0;
        self.transform_point(point)
    }
}

impl FnMut<(Point3,)> for Transform {
    extern "rust-call" fn call_mut(&mut self, point: (Point3,)) -> Point3 {
        self.call(point)
    }
}

impl FnOnce<(Point3,)> for Transform {
    type Output = Point3;

    extern "rust-call" fn call_once(self, point: (Point3,)) -> Self::Output {
        self.call(point)
    }
}

impl Fn<(Vector3,)> for Transform {
    extern "rust-call" fn call(&self, args: (Vector3,)) -> Vector3 {
        let vector = args.0;
        let x = vector.x;
        let y = vector.y;
        let z = vector.z;
        let m = self.matrix.m;

        Vector3 {
            x: m[0][0] * x + m[0][1] * y + m[0][2] * z,
            y: m[1][0] * x + m[1][1] * y + m[1][2] * z,
            z: m[2][0] * x + m[2][1] * y + m[2][2] * z,
        }
    }
}

impl FnMut<(Vector3,)> for Transform {
    extern "rust-call" fn call_mut(&mut self, vector: (Vector3,)) -> Vector3 {
        self.call(vector)
    }
}

impl FnOnce<(Vector3,)> for Transform {
    type Output = Vector3;

    extern "rust-call" fn call_once(self, vector: (Vector3,)) -> Self::Output {
        self.call(vector)
    }
}

impl Fn<(Normal3,)> for Transform {
    extern "rust-call" fn call(&self, args: (Normal3,)) -> Normal3 {
        let normal = args.0;
        let x = normal.x;
        let y = normal.y;
        let z = normal.z;
        let m = self.matrix.m;

        Normal3 {
            x: m[0][0] * x + m[0][1] * y + m[0][2] * z,
            y: m[1][0] * x + m[1][1] * y + m[1][2] * z,
            z: m[2][0] * x + m[2][1] * y + m[2][2] * z,
        }
    }
}

impl FnMut<(Normal3,)> for Transform {
    extern "rust-call" fn call_mut(&mut self, normal: (Normal3,)) -> Normal3 {
        self.call(normal)
    }
}

impl FnOnce<(Normal3,)> for Transform {
    type Output = Normal3;

    extern "rust-call" fn call_once(self, normal: (Normal3,)) -> Self::Output {
        self.call(normal)
    }
}

impl Fn<(Ray,)> for Transform {
    extern "rust-call" fn call(&self, args: (Ray,)) -> Ray {
        let ray = args.0;
        let (mut origin, o_error) = self.transform_point_with_error(ray.origin);
        let direction = self(ray.direction);

        let length_squared = direction.length_squared();
        let mut max_length = ray.max_length;
        if length_squared > 0.0 {
            let dt = Vector3::dot(&(Vector3::abs(direction)), o_error) / length_squared;
            origin = origin + (direction * dt);
            max_length -= dt;
        }

        Ray {
            origin,
            direction,
            max_length,
            time: ray.time,
            // TODO medium: ray.medium,
        }
    }
}

impl FnMut<(Ray,)> for Transform {
    extern "rust-call" fn call_mut(&mut self, ray: (Ray,)) -> Ray {
        self.call(ray)
    }
}

impl FnOnce<(Ray,)> for Transform {
    type Output = Ray;

    extern "rust-call" fn call_once(self, ray: (Ray,)) -> Self::Output {
        self.call(ray)
    }
}

/// Transforms Bound
///
/// Example:
/// ```
/// use rpbr::foundation::geometry::bounds::Bounds3;
/// use rpbr::foundation::geometry::point::Point3;
/// use rpbr::foundation::transform::transform::Transform;
///
/// let bounds = Bounds3::new(Point3::zero(), Point3::new(1.0, 2.0, 3.0));
/// let transform = &Transform::scale(1.0, 1.0, 1.0);
///
/// let actual = (*transform)(bounds);
/// let expected = bounds;
///
/// assert_eq!(actual, expected);
/// ```
impl Fn<(Bounds3,)> for Transform {
    extern "rust-call" fn call(&self, bounds: (Bounds3,)) -> Bounds3 {
        let b = bounds.0;
        let mut ret = Bounds3::from(self(Point3 {
            x: b.p_min.y,
            y: b.p_min.y,
            z: b.p_min.z,
        }));
        ret = ret.union_point(self(Point3 {
            x: b.p_max.x,
            y: b.p_min.y,
            z: b.p_min.z,
        }));
        ret = ret.union_point(self(Point3 {
            x: b.p_min.y,
            y: b.p_max.y,
            z: b.p_min.z,
        }));
        ret = ret.union_point(self(Point3 {
            x: b.p_min.y,
            y: b.p_min.y,
            z: b.p_max.z,
        }));
        ret = ret.union_point(self(Point3 {
            x: b.p_min.y,
            y: b.p_max.y,
            z: b.p_max.z,
        }));
        ret = ret.union_point(self(Point3 {
            x: b.p_max.x,
            y: b.p_max.y,
            z: b.p_min.z,
        }));
        ret = ret.union_point(self(Point3 {
            x: b.p_max.x,
            y: b.p_min.y,
            z: b.p_max.z,
        }));
        ret = ret.union_point(self(Point3 {
            x: b.p_max.x,
            y: b.p_max.y,
            z: b.p_max.z,
        }));

        ret
    }
}

impl FnMut<(Bounds3,)> for Transform {
    extern "rust-call" fn call_mut(&mut self, bounds: (Bounds3,)) -> Bounds3 {
        self.call(bounds)
    }
}

impl FnOnce<(Bounds3,)> for Transform {
    type Output = Bounds3;

    extern "rust-call" fn call_once(self, bounds: (Bounds3,)) -> Self::Output {
        self.call(bounds)
    }
}

impl Mul for Transform {
    type Output = Transform;

    fn mul(self, other: Self) -> Self::Output {
        Self {
            matrix: self.matrix * other.matrix,
            inverse: other.inverse * self.inverse,
        }
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
