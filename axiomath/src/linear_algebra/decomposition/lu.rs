//! LU decomposition with partial pivoting.
//!
//! LU decomposition factors a matrix A as A = PLU, where:
//! - P is a permutation matrix
//! - L is lower triangular with ones on the diagonal
//! - U is upper triangular
//!
//! # Algorithm
//!
//! We use the Doolittle algorithm with partial pivoting for numerical stability:
//!
//! 1. For each column k:
//!    - Find the pivot (largest absolute value in column k, rows k to n)
//!    - Swap rows if necessary (updates P)
//!    - Compute multipliers for rows below pivot
//!    - Update remaining submatrix
//!
//! # Applications
//!
//! - Solving linear systems Ax = b
//! - Computing determinants: det(A) = det(P) * det(L) * det(U)
//! - Matrix inversion
//!
//! # Complexity
//!
//! - Time: O(n³) for an n×n matrix
//! - Space: O(n²) for L, U, and P matrices
//!
//! # References
//!
//! - Golub & Van Loan, "Matrix Computations", Algorithm 3.4.1
//! - Trefethen & Bau, "Numerical Linear Algebra", Lecture 20

use crate::core::traits::Float;
use crate::linear_algebra::LinearAlgebraError;
use crate::matrix::Matrix;
use num_traits::Zero;

/// Performs LU decomposition with partial pivoting using the Doolittle algorithm.
///
/// Factors a square matrix A into A = PLU where:
/// - P is a permutation matrix (represented as a permutation vector)
/// - L is lower triangular with ones on the diagonal
/// - U is upper triangular
///
/// # Mathematical Definition
///
/// ```text
/// PA = LU
/// ```
///
/// where:
/// - P represents row permutations applied during pivoting
/// - L_ii = 1 for all i (unit lower triangular)
/// - L_ij = 0 for j > i
/// - U_ij = 0 for i > j
///
/// # Algorithm
///
/// Doolittle algorithm with partial pivoting:
///
/// ```text
/// for k = 1 to n:
///     Find pivot: i = argmax_{j≥k} |A_jk|
///     Swap rows k and i
///     for j = k+1 to n:
///         L_jk = A_jk / A_kk
///         for m = k to n:
///             A_jm = A_jm - L_jk * A_km
/// ```
///
/// # Numerical Stability
///
/// Partial pivoting ensures that all multipliers |L_ij| ≤ 1, which provides
/// good numerical stability in practice. The growth factor is typically small
/// for most matrices encountered in practice.
///
/// # Complexity
///
/// - Time: O(n³)
/// - Space: O(n²)
///
/// # Arguments
///
/// * `matrix` - Square matrix to decompose
///
/// # Returns
///
/// Returns `Ok((L, U, perm))` where:
/// - `L` is the lower triangular factor (with ones on diagonal)
/// - `U` is the upper triangular factor
/// - `perm` is the permutation vector representing P
///
/// Returns `Err(LinearAlgebraError::SingularMatrix)` if the matrix is singular
/// or numerically singular.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::decomposition::lu_decompose;
/// use axiomath::prelude::*;
///
/// let a = Matrix::from_array([
///     [2.0, 1.0, 1.0],
///     [4.0, 3.0, 3.0],
///     [8.0, 7.0, 9.0],
/// ]);
///
/// let (l, u, perm) = lu_decompose(&a).unwrap();
///
/// // Verify: PA = LU
/// // (multiply L * U and check against permuted A)
/// ```
///
/// # Panics
///
/// Panics if the matrix is not square.
pub fn lu_decompose<T: Float, const N: usize>(
    matrix: &Matrix<T, N, N>,
) -> Result<(Matrix<T, N, N>, Matrix<T, N, N>, Vec<usize>), LinearAlgebraError> {
    // Initialize L as identity, U as copy of A, and permutation as identity
    let mut l = Matrix::<T, N, N>::identity();
    let mut u = *matrix;
    let mut perm: Vec<usize> = (0..N).collect();

    // Tolerance for singularity check
    let epsilon = T::epsilon() * T::from(N).unwrap() * T::from(100.0).unwrap();

    // Doolittle algorithm with partial pivoting
    for k in 0..N {
        // Find pivot: row with largest absolute value in column k
        let mut pivot_row = k;
        let mut max_val = u[(k, k)].abs();

        for i in (k + 1)..N {
            let abs_val = u[(i, k)].abs();
            if abs_val > max_val {
                max_val = abs_val;
                pivot_row = i;
            }
        }

        // Check for singularity
        if max_val < epsilon {
            return Err(LinearAlgebraError::SingularMatrix);
        }

        // Swap rows if necessary
        if pivot_row != k {
            // Swap rows in U
            for j in 0..N {
                let temp = u[(k, j)];
                u[(k, j)] = u[(pivot_row, j)];
                u[(pivot_row, j)] = temp;
            }

            // Swap corresponding rows in L (only the already-computed part)
            for j in 0..k {
                let temp = l[(k, j)];
                l[(k, j)] = l[(pivot_row, j)];
                l[(pivot_row, j)] = temp;
            }

            // Update permutation vector
            perm.swap(k, pivot_row);
        }

        // Compute multipliers and eliminate
        for i in (k + 1)..N {
            // Compute multiplier
            let multiplier = u[(i, k)] / u[(k, k)];
            l[(i, k)] = multiplier;

            // Eliminate: row i -= multiplier * row k
            for j in k..N {
                u[(i, j)] = u[(i, j)] - multiplier * u[(k, j)];
            }
        }
    }

    Ok((l, u, perm))
}

/// Computes the determinant of a matrix using LU decomposition.
///
/// The determinant is computed as:
/// ```text
/// det(A) = sign(P) * det(L) * det(U) = sign(P) * product(U_ii)
/// ```
///
/// Since L has ones on the diagonal, det(L) = 1.
///
/// # Complexity
///
/// - Time: O(n³) (dominated by LU decomposition)
///
/// # Arguments
///
/// * `matrix` - Square matrix
///
/// # Returns
///
/// The determinant of the matrix, or `Err` if the matrix is singular.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::decomposition::lu::determinant_lu;
///
/// let a = Matrix::from_array([
///     [1.0, 2.0],
///     [3.0, 4.0],
/// ]);
///
/// let det = determinant_lu(&a).unwrap();
/// assert!((det - (-2.0)).abs() < 1e-10);
/// ```
pub fn determinant_lu<T: Float, const N: usize>(
    matrix: &Matrix<T, N, N>,
) -> Result<T, LinearAlgebraError> {
    let (_, u, perm) = lu_decompose(matrix)?;

    // Compute sign of permutation
    let mut sign = T::one();
    for i in 0..N {
        if perm[i] != i {
            // Count inversions
            for j in (i + 1)..N {
                if perm[j] < perm[i] {
                    sign = -sign;
                }
            }
        }
    }

    // Compute product of diagonal of U
    let mut det = sign;
    for i in 0..N {
        det = det * u[(i, i)];
    }

    Ok(det)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::multiply;
    use crate::utils::comparison::approximately_equal;

    #[test]
    fn test_lu_decompose_3x3() {
        let a = Matrix::from_rows([[2.0, 1.0, 1.0], [4.0, 3.0, 3.0], [8.0, 7.0, 9.0]]);

        let (l, u, perm) = lu_decompose(&a).unwrap();

        // L should be lower triangular with ones on diagonal
        for i in 0..3 {
            assert!((l[(i, i)] - 1.0_f64).abs() < 1e-10);
            for j in (i + 1)..3 {
                assert!(l[(i, j)].abs() < 1e-10_f64);
            }
        }

        // U should be upper triangular
        for i in 1..3 {
            for j in 0..i {
                assert!(u[(i, j)].abs() < 1e-10_f64);
            }
        }

        // Verify PA = LU by reconstructing
        let lu = multiply(&l, &u);

        // Apply permutation to original matrix
        let mut pa = Matrix::<f64, 3, 3>::zero();
        for i in 0..3 {
            for j in 0..3 {
                pa[(i, j)] = a[(perm[i], j)];
            }
        }

        // Check PA ≈ LU
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    approximately_equal(pa[(i, j)], lu[(i, j)], 1e-10, 1e-10),
                    "PA[{},{}] = {}, LU[{},{}] = {}",
                    i,
                    j,
                    pa[(i, j)],
                    i,
                    j,
                    lu[(i, j)]
                );
            }
        }
    }

    #[test]
    fn test_lu_decompose_identity() {
        let a = Matrix::<f64, 3, 3>::identity();
        let (l, u, perm) = lu_decompose(&a).unwrap();

        // For identity matrix: L = I, U = I, P = I
        for i in 0..3 {
            assert_eq!(perm[i], i);
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((l[(i, j)] - expected).abs() < 1e-10);
                assert!((u[(i, j)] - expected).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_lu_decompose_singular() {
        let a = Matrix::from_rows([[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [1.0, 1.0, 1.0]]);

        let result = lu_decompose(&a);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), LinearAlgebraError::SingularMatrix);
    }

    #[test]
    fn test_lu_determinant_2x2() {
        let a = Matrix::from_rows([[1.0, 2.0], [3.0, 4.0]]);

        let det = determinant_lu(&a).unwrap();
        // det = 1*4 - 2*3 = -2
        assert!((det - (-2.0_f64)).abs() < 1e-10);
    }

    #[test]
    fn test_lu_determinant_3x3() {
        let a = Matrix::from_rows([[2.0, 1.0, 1.0], [4.0, 3.0, 3.0], [8.0, 7.0, 9.0]]);

        let det = determinant_lu(&a).unwrap();
        // Compute manually: 2*(3*9-3*7) - 1*(4*9-3*8) + 1*(4*7-3*8) = 2*6 - 1*12 + 1*4 = 4
        assert!((det - 4.0_f64).abs() < 1e-10);
    }

    #[test]
    fn test_lu_decompose_with_pivoting() {
        // Matrix that requires pivoting
        let a = Matrix::from_rows([[0.0, 1.0], [1.0, 1.0]]);

        let (l, u, perm) = lu_decompose(&a).unwrap();

        // Verify PA = LU
        let lu = multiply(&l, &u);
        let mut pa = Matrix::<f64, 2, 2>::zero();
        for i in 0..2 {
            for j in 0..2 {
                pa[(i, j)] = a[(perm[i], j)];
            }
        }

        for i in 0..2 {
            for j in 0..2 {
                assert!(approximately_equal(pa[(i, j)], lu[(i, j)], 1e-10, 1e-10));
            }
        }
    }
}
