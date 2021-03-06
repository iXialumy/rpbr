use std::ops::Mul;
use std::ptr::swap;

use crate::foundation::pbr::Float;

#[derive(Debug, Copy, Clone, PartialOrd, PartialEq)]
pub struct Matrix4x4 {
    pub m: [[Float; 4]; 4],
}

impl Matrix4x4 {
    #[allow(clippy::too_many_arguments)]
    /// Create a matrix from 16 Floats
    pub fn matrix_from_floats(
        t00: Float,
        t01: Float,
        t02: Float,
        t03: Float,
        t10: Float,
        t11: Float,
        t12: Float,
        t13: Float,
        t20: Float,
        t21: Float,
        t22: Float,
        t23: Float,
        t30: Float,
        t31: Float,
        t32: Float,
        t33: Float,
    ) -> Self {
        Self {
            m: [
                [t00, t01, t02, t03],
                [t10, t11, t12, t13],
                [t20, t21, t22, t23],
                [t30, t31, t32, t33],
            ],
        }
    }

    /// Initializes a new identity matrix
    pub fn identity() -> Self {
        Self {
            m: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }

    /// Floatranspose the Matrix
    pub fn transpose(self) -> Self {
        Self {
            m: [
                [self.m[0][0], self.m[1][0], self.m[2][0], self.m[3][0]],
                [self.m[0][1], self.m[1][1], self.m[2][1], self.m[3][1]],
                [self.m[0][2], self.m[1][2], self.m[2][2], self.m[3][2]],
                [self.m[0][3], self.m[1][3], self.m[2][3], self.m[3][3]],
            ],
        }
    }

    pub fn inverse(self) -> Self {
        let mut indxc = [0i32; 4];
        let mut indxr = [0i32; 4];
        let mut ipiv = [0i32; 4];

        let mut minv = self.m;

        for i in 0..4 {
            let mut irow = 0;
            let mut icol = 0;
            let mut big = 0.0;

            // Choose pivot
            for j in 0..4 {
                if ipiv[j] != 1 {
                    for (k, &pivot) in ipiv.iter().enumerate() {
                        if pivot == 0 {
                            if minv[j][k].abs() >= big {
                                big = (minv[j][k]).abs();
                                irow = j;
                                icol = k;
                            }
                        } else if ipiv[k] > 1 {
                            panic!("Singular matrix in MatrixInvert");
                        }
                    }
                }
            }
            ipiv[icol] += 1;

            // Swap rows irow and icol for pivot
            if irow != icol {
                for k in 0..4 {
                    unsafe { swap(&mut minv[irow][k], &mut minv[icol][k]) }
                }
            }

            indxr[i] = irow as i32;
            indxc[i] = icol as i32;
            if minv[icol][icol] == 0.0 {
                panic!("Singular matrix in MatrixInvert")
            }

            // Set self.m[icol][icol] to one by scaling row icol appropriately
            let pivinv = 1.0 / minv[icol][icol];
            minv[icol][icol] = 1.0;
            for j in 0..4 {
                minv[icol][j] *= pivinv;
            }

            // Subtract this row from others to zero out their columns
            for j in 0..4 {
                if j != icol {
                    let save = minv[j][icol];
                    minv[j][icol] = 0.0;
                    for k in 0..4 {
                        minv[j][k] -= minv[icol][k] * save
                    }
                }
            }
        }
        // Swap columns to reflect permutation
        #[allow(clippy::needless_range_loop)]
        for j in (0..4).rev() {
            if indxr[j] != indxc[j] {
                for k in 0..4 {
                    unsafe {
                        swap(
                            &mut minv[k][indxr[j] as usize],
                            &mut minv[k][indxc[j] as usize],
                        )
                    }
                }
            }
        }

        Self { m: minv }
    }

    pub fn det(&self) -> Float {
        self.m[0][0] * (self.m[1][1] * self.m[2][2] - self.m[1][2] * self.m[2][1])
            - self.m[0][1] * (self.m[1][0] * self.m[2][2] - self.m[1][2] * self.m[2][0])
            + self.m[0][2] * (self.m[1][0] * self.m[2][1] - self.m[1][1] * self.m[2][0])
    }
}

impl Mul for Matrix4x4 {
    type Output = Self;

    #[allow(clippy::needless_range_loop)]
    fn mul(self, rhs: Self) -> Self::Output {
        let mut array: [[Float; 4]; 4] = Matrix4x4::identity().m;
        for i in 0..4 {
            for j in 0..4 {
                array[i][j] = self.m[i][0] * rhs.m[0][j]
                    + self.m[i][1] * rhs.m[1][j]
                    + self.m[i][2] * rhs.m[2][j]
                    + self.m[i][3] * rhs.m[3][j];
            }
        }
        Self { m: array }
    }
}

impl From<[[Float; 4]; 4]> for Matrix4x4 {
    fn from(array: [[Float; 4]; 4]) -> Self {
        Self { m: array }
    }
}

#[cfg(test)]
mod tests {
    use crate::foundation::transforms::matrix::Matrix4x4;

    #[test]
    pub fn invert_matrix() {
        let matrix = Matrix4x4 {
            m: [
                [2.0, -1.0, 0.0, 0.0],
                [1.0, 2.0, -2.0, 0.0],
                [0.0, -1.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        };

        assert_eq!(matrix * matrix.inverse(), Matrix4x4::identity());
    }

    #[test]
    pub fn transpose_matrix() {
        let matrix = Matrix4x4 {
            m: [
                [2.0, -1.0, 0.0, 0.0],
                [1.0, 2.0, -2.0, 0.0],
                [0.0, -1.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        };

        let expected = Matrix4x4 {
            m: [
                [2.0, 1.0, 0.0, 0.0],
                [-1.0, 2.0, -1.0, 0.0],
                [0.0, -2.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        };

        assert_eq!(matrix.transpose(), expected);
    }
}
