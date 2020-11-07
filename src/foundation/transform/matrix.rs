use std::ops::Mul;
use std::ptr::swap;

use num_traits::{Float, FromPrimitive};

struct Matrix4x4<T: Float + FromPrimitive> {
    m: [[T; 4]; 4],
}

impl<T: Float + FromPrimitive> Matrix4x4<T> {
    /// Initializes a new identity matrix
    pub fn identity() -> Self {
        Self {
            m: [
                [
                    T::from_f64(1.0).unwrap(),
                    T::from_f64(0.0).unwrap(),
                    T::from_f64(0.0).unwrap(),
                    T::from_f64(0.0).unwrap(),
                ],
                [
                    T::from_f64(0.0).unwrap(),
                    T::from_f64(1.0).unwrap(),
                    T::from_f64(0.0).unwrap(),
                    T::from_f64(0.0).unwrap(),
                ],
                [
                    T::from_f64(0.0).unwrap(),
                    T::from_f64(0.0).unwrap(),
                    T::from_f64(1.0).unwrap(),
                    T::from_f64(0.0).unwrap(),
                ],
                [
                    T::from_f64(0.0).unwrap(),
                    T::from_f64(0.0).unwrap(),
                    T::from_f64(0.0).unwrap(),
                    T::from_f64(1.0).unwrap(),
                ],
            ],
        }
    }

    /// Transpose the Matrix
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
            let mut big = T::from_f64(0.0).unwrap();

            // Choose pivot
            for j in 0..4 {
                if ipiv[j] != 1 {
                    for k in 0..4 {
                        if ipiv[k] == 0 {
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
            unsafe {
                if irow != icol {
                    for k in 0..4 {
                        swap(&mut minv[irow][k], &mut minv[icol][k])
                    }
                }
            }
            indxr[i] = irow as i32;
            indxc[i] = icol as i32;
            if minv[icol][icol] == T::from_f64(0.0).unwrap() {
                panic!("Singular matrix in MatrixInvert")
            }

            // Set self.m[icol][icol] to one by scaling row icol appropriately
            let pivinv = T::from_f64(1.0).unwrap() / minv[icol][icol];
            for j in 0..4 {
                minv[icol][j] = minv[icol][j] * pivinv;
            }

            // Subtract this row from others to zero out their columns
            for j in 0..4 {
                if j != icol {
                    let save = minv[j][icol];
                    minv[j][icol] = T::from_f64(0.0).unwrap();
                    for k in 0..4 {
                        minv[j][k] = minv[j][k] - minv[icol][k] * save
                    }
                }
            }
        }
        // Swap columns to reflect permutation
        unsafe {
            for j in 3..=0 {
                if indxr[j] != indxc[j] {
                    for k in 0..4 {
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
}

impl<T: Float + FromPrimitive> Mul for Matrix4x4<T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut array: [[T; 4]; 4] = Matrix4x4::identity().m;
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