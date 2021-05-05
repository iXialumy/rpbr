use std::cmp::Ordering;
use std::fmt::Debug;
use std::ops::Mul;

use num_traits::{abs, AsPrimitive, Float, FromPrimitive, Signed};

use crate::foundation::geometry::bounds::Bounds3;
use crate::foundation::geometry::normal::Normal3;
use crate::foundation::geometry::point::Point3;
use crate::foundation::geometry::ray::Ray;
use crate::foundation::geometry::vector::Vector3;
use crate::foundation::transform::matrix::Matrix4x4;
use crate::foundation::util::gamma;
use num_traits::real::Real;

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Transform<T: Float + FromPrimitive> {
    pub matrix: Matrix4x4<T>,
    pub inverse: Matrix4x4<T>,
}

impl<T: Float + FromPrimitive + AsPrimitive<f64> + Copy + Debug + Signed> Transform<T> {
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

    pub fn transform_point_with_error(self, point: Point3<T>) -> (Point3<T>, Vector3<T>) {
        let x = point.x;
        let y = point.y;
        let z = point.z;
        let m = self.matrix.m;

        let xp = m[0][0] * x + m[0][1] * y + m[0][2] * z + m[0][3];
        let yp = m[1][0] * x + m[1][1] * y + m[1][2] * z + m[1][3];
        let zp = m[2][0] * x + m[2][1] * y + m[2][2] * z + m[2][3];
        let wp = m[3][0] * x + m[3][1] * y + m[3][2] * z + m[3][3];

        let x_abs_sum = Float::abs(m[0][0] * x)
            + Float::abs(m[0][1] * y)
            + Float::abs(m[0][2] * z)
            + Float::abs(m[0][3]);
        let y_abs_sum = Float::abs(m[1][0] * x)
            + Float::abs(m[1][1] * y)
            + Float::abs(m[1][2] * z)
            + Float::abs(m[1][3]);
        let z_abs_sum = Float::abs(m[2][0] * x)
            + Float::abs(m[2][1] * y)
            + Float::abs(m[2][2] * z)
            + Float::abs(m[2][3]);

        let error: Vector3<T> =
            Vector3::new(x_abs_sum, y_abs_sum, z_abs_sum) * gamma(T::from_f64(3.0).unwrap());

        debug_assert_eq!(wp, T::from_f64(0.0).unwrap());

        let p = Point3 {
            x: xp,
            y: yp,
            z: zp,
        };
        if wp == T::from_f64(1.0).unwrap() {
            (p, error)
        } else {
            (p / wp, error)
        }
    }

    pub fn transform_vector_with_error(self, vector: Vector3<T>) -> (Vector3<T>, Vector3<T>) {
        let err_x = gamma(T::from_f64(3.0).unwrap())
            * (abs(self.matrix.m[0][0] * vector.x)
                + abs(self.matrix.m[0][1] * vector.y)
                + abs(self.matrix.m[0][2] * vector.z));

        let err_y = gamma(T::from_f64(3.0).unwrap())
            * (abs(self.matrix.m[1][0] * vector.x)
                + abs(self.matrix.m[1][1] * vector.y)
                + abs(self.matrix.m[1][2] * vector.z));

        let err_z = gamma(T::from_f64(3.0).unwrap())
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

    pub fn transform_ray_with_error(self, ray: Ray<T>) -> (Ray<T>, Vector3<T>, Vector3<T>) {
        let (mut origin, origin_error) = self.transform_point_with_error(ray.origin);
        let (direction, direction_error) = self.transform_vector_with_error(ray.direction);

        let max_length = ray.max_length;
        let length_squared = direction.length_squared();

        if length_squared > T::zero() {
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

    /// Tells if a transformatin will swap the handedness of a coordinate system
    pub fn swaps_handedness(self) -> bool {
        self.matrix.det() < T::from_f64(0.0).unwrap()
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64> + Debug + Signed> From<[[T; 4]; 4]>
    for Transform<T>
{
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
impl<T: Float + FromPrimitive + AsPrimitive<f64> + Copy + Debug> Fn<(Point3<T>,)> for Transform<T> {
    extern "rust-call" fn call(&self, args: (Point3<T>,)) -> Point3<T> {
        let point = args.0;
        let x = point.x;
        let y = point.y;
        let z = point.z;
        let m = self.matrix.m;

        let xp = m[0][0] * x + m[0][1] * y + m[0][2] * z + m[0][3];
        let yp = m[1][0] * x + m[1][1] * y + m[1][2] * z + m[1][3];
        let zp = m[2][0] * x + m[2][1] * y + m[2][2] * z + m[2][3];
        let wp = m[3][0] * x + m[3][1] * y + m[3][2] * z + m[3][3];

        assert_ne!(wp, T::from_f64(0.0).unwrap());
        if wp == T::from_f64(0.0).unwrap() {
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
}

impl<T: Float + FromPrimitive + AsPrimitive<f64> + Copy + Debug> FnMut<(Point3<T>,)>
    for Transform<T>
{
    extern "rust-call" fn call_mut(&mut self, point: (Point3<T>,)) -> Point3<T> {
        self.call(point)
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64> + Copy + Debug> FnOnce<(Point3<T>,)>
    for Transform<T>
{
    type Output = Point3<T>;

    extern "rust-call" fn call_once(self, point: (Point3<T>,)) -> Self::Output {
        self.call(point)
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64> + Copy + Debug> Fn<(Vector3<T>,)>
    for Transform<T>
{
    extern "rust-call" fn call(&self, args: (Vector3<T>,)) -> Vector3<T> {
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

impl<T: Float + FromPrimitive + AsPrimitive<f64> + Copy + Debug> FnMut<(Vector3<T>,)>
    for Transform<T>
{
    extern "rust-call" fn call_mut(&mut self, vector: (Vector3<T>,)) -> Vector3<T> {
        self.call(vector)
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64> + Copy + Debug> FnOnce<(Vector3<T>,)>
    for Transform<T>
{
    type Output = Vector3<T>;

    extern "rust-call" fn call_once(self, vector: (Vector3<T>,)) -> Self::Output {
        self.call(vector)
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64> + Debug> Fn<(Normal3<T>,)> for Transform<T> {
    extern "rust-call" fn call(&self, args: (Normal3<T>,)) -> Normal3<T> {
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

impl<T: Float + FromPrimitive + AsPrimitive<f64> + Copy + Debug> FnMut<(Normal3<T>,)>
    for Transform<T>
{
    extern "rust-call" fn call_mut(&mut self, normal: (Normal3<T>,)) -> Normal3<T> {
        self.call(normal)
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64> + Copy + Debug> FnOnce<(Normal3<T>,)>
    for Transform<T>
{
    type Output = Normal3<T>;

    extern "rust-call" fn call_once(self, normal: (Normal3<T>,)) -> Self::Output {
        self.call(normal)
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64> + Debug + Signed> Fn<(Ray<T>,)> for Transform<T> {
    extern "rust-call" fn call(&self, args: (Ray<T>,)) -> Ray<T> {
        let ray = args.0;
        let (mut origin, o_error) = self.transform_point_with_error(ray.origin);
        let direction = self(ray.direction);

        let length_squared = direction.length_squared();
        let mut max_length = ray.max_length;
        if length_squared > T::from_f64(0.0).unwrap() {
            let dt = Vector3::dot(&direction.abs(), o_error) / length_squared;
            origin = origin + (direction * dt);
            max_length = max_length - dt;
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

impl<T: Float + FromPrimitive + AsPrimitive<f64> + Copy + Debug + Signed> FnMut<(Ray<T>,)>
    for Transform<T>
{
    extern "rust-call" fn call_mut(&mut self, ray: (Ray<T>,)) -> Ray<T> {
        self.call(ray)
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64> + Copy + Debug + Signed> FnOnce<(Ray<T>,)>
    for Transform<T>
{
    type Output = Ray<T>;

    extern "rust-call" fn call_once(self, ray: (Ray<T>,)) -> Self::Output {
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
/// let bounds = Bounds3::new(Point3::default(), Point3::new(1.0, 2.0, 3.0));
/// let transform = &Transform::scale(1.0, 1.0, 1.0);
///
/// let actual = (*transform)(bounds);
/// let expected = bounds;
///
/// assert_eq!(actual, expected);
/// ```
impl<T: Float + FromPrimitive + AsPrimitive<f64> + Copy + Debug> Fn<(Bounds3<T>,)>
    for Transform<T>
{
    extern "rust-call" fn call(&self, bounds: (Bounds3<T>,)) -> Bounds3<T> {
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

impl<T: Float + FromPrimitive + AsPrimitive<f64> + Copy + Debug> FnMut<(Bounds3<T>,)>
    for Transform<T>
{
    extern "rust-call" fn call_mut(&mut self, bounds: (Bounds3<T>,)) -> Bounds3<T> {
        self.call(bounds)
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64> + Copy + Debug> FnOnce<(Bounds3<T>,)>
    for Transform<T>
{
    type Output = Bounds3<T>;

    extern "rust-call" fn call_once(self, bounds: (Bounds3<T>,)) -> Self::Output {
        self.call(bounds)
    }
}

impl<T: Float + FromPrimitive + AsPrimitive<f64> + Copy + Debug> Mul for Transform<T> {
    type Output = Transform<T>;

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
