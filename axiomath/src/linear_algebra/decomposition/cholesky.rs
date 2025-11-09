//! Cholesky decomposition for symmetric positive definite matrices.
//!
//! Cholesky decomposition factors a symmetric positive definite matrix A as:
//! A = LL^T
//! where L is lower triangular with positive diagonal elements.
//!
//! # Advantages
//!
//! - About 2× faster than LU decomposition (only n³/3 operations vs 2n³/3)
//! - Numerically more stable than LU for positive definite matrices
//! - Useful for solving linear systems and computing determinants
//!
//! # Requirements
//!
//! The matrix must be:
//! 1. Symmetric: A = A^T
//! 2. Positive definite: x^T A x > 0 for all non-zero x
//!
//! # Applications
//!
//! - Solving linear systems with symmetric positive definite matrices
//! - Computing determinants: det(A) = (product of L_ii)²
//! - Generating correlated random variables
//! - Optimization algorithms
//!
//! # References
//!
//! - Golub & Van Loan, "Matrix Computations", Section 4.2
//! - Trefethen & Bau, "Numerical Linear Algebra", Lecture 23

use crate::core::traits::Float;
use crate::linear_algebra::LinearAlgebraError;
use crate::matrix::Matrix;
use num_traits::Zero;

/// Performs Cholesky decomposition on a symmetric positive definite matrix.
///
/// Factors A = LL^T where L is lower triangular with positive diagonal elements.
///
/// # Mathematical Definition
///
/// ```text
/// A = LL^T
/// ```
///
/// where:
/// - L is lower triangular: L_ij = 0 for j > i
/// - L_ii > 0 for all i (positive diagonal)
///
/// # Algorithm
///
/// Cholesky-Banachiewicz algorithm:
///
/// ```text
/// for j = 1 to n:
///     L_jj = sqrt(A_jj - sum_{k=1}^{j-1} L_jk²)
///     for i = j+1 to n:
///         L_ij = (A_ij - sum_{k=1}^{j-1} L_ik * L_jk) / L_jj
/// ```
///
/// # Numerical Stability
///
/// Cholesky decomposition is numerically stable for positive definite matrices.
/// If the algorithm encounters a negative value under the square root, the
/// matrix is not positive definite.
///
/// # Complexity
///
/// - Time: O(n³/3) ≈ half of LU decomposition
/// - Space: O(n²)
///
/// # Arguments
///
/// * `matrix` - Symmetric positive definite matrix (n×n)
///
/// # Returns
///
/// Returns `Ok(L)` where L is the lower triangular Cholesky factor.
///
/// Returns `Err` if:
/// - Matrix is not symmetric
/// - Matrix is not positive definite
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::decomposition::cholesky_decompose;
///
/// // Symmetric positive definite matrix
/// let a = Matrix::from_rows([
///     [4.0, 12.0, -16.0],
///     [12.0, 37.0, -43.0],
///     [-16.0, -43.0, 98.0],
/// ]);
///
/// let l = cholesky_decompose(&a).unwrap();
///
/// // Verify: A = LL^T
/// ```
pub fn cholesky_decompose<T: Float, const N: usize>(
    matrix: &Matrix<T, N, N>,
) -> Result<Matrix<T, N, N>, LinearAlgebraError> {
    let epsilon = T::epsilon() * T::from(N).unwrap() * T::from(100.0).unwrap();

    // Check if matrix is symmetric
    for i in 0..N {
        for j in (i + 1)..N {
            if (matrix[(i, j)] - matrix[(j, i)]).abs() > epsilon {
                return Err(LinearAlgebraError::NotSymmetric);
            }
        }
    }

    let mut l = Matrix::<T, N, N>::zero();

    // Cholesky-Banachiewicz algorithm
    for j in 0..N {
        // Compute L_jj
        let mut sum = matrix[(j, j)];
        for k in 0..j {
            sum = sum - l[(j, k)] * l[(j, k)];
        }

        // Check for positive definiteness
        if sum <= T::zero() {
            return Err(LinearAlgebraError::NotPositiveDefinite);
        }

        l[(j, j)] = sum.sqrt();

        // Compute L_ij for i > j
        for i in (j + 1)..N {
            let mut sum = matrix[(i, j)];
            for k in 0..j {
                sum = sum - l[(i, k)] * l[(j, k)];
            }
            l[(i, j)] = sum / l[(j, j)];
        }
    }

    Ok(l)
}

/// Solves a linear system Ax = b where A is symmetric positive definite using Cholesky decomposition.
///
/// This is more efficient than LU decomposition for SPD matrices.
///
/// # Algorithm
///
/// 1. Compute Cholesky decomposition: A = LL^T
/// 2. Solve Ly = b via forward substitution
/// 3. Solve L^T x = y via backward substitution
///
/// # Complexity
///
/// - Time: O(n³/3) for decomposition + O(n²) for solving = O(n³/3)
/// - About 2× faster than LU-based solving
///
/// # Arguments
///
/// * `matrix` - Symmetric positive definite matrix A (n×n)
/// * `b` - Right-hand side vector (n)
///
/// # Returns
///
/// Solution vector x such that Ax = b.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::vector::Vector;
/// use axiomath::linear_algebra::decomposition::cholesky::cholesky_solve;
///
/// let a = Matrix::from_rows([
///     [4.0, 2.0],
///     [2.0, 3.0],
/// ]);
///
/// let b = Vector::from([10.0, 9.0]);
/// let x = cholesky_solve(&a, &b).unwrap();
/// ```
pub fn cholesky_solve<T: Float, const N: usize>(
    matrix: &Matrix<T, N, N>,
    b: &crate::vector::Vector<T, N>,
) -> Result<crate::vector::Vector<T, N>, LinearAlgebraError> {
    use crate::linear_algebra::systems::triangular::{
        solve_lower_triangular, solve_upper_triangular,
    };
    use crate::matrix::transpose;

    // Step 1: Compute Cholesky decomposition
    let l = cholesky_decompose(matrix)?;

    // Step 2: Solve Ly = b
    let y = solve_lower_triangular(&l, b)?;

    // Step 3: Solve L^T x = y
    let lt = transpose(&l);
    let x = solve_upper_triangular(&lt, &y)?;

    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::{multiply, transpose};
    use crate::utils::comparison::approximately_equal;
    use crate::vector::Vector;

    #[test]
    fn test_cholesky_decompose_3x3() {
        // Symmetric positive definite matrix
        let a = Matrix::<f64, 3, 3>::from_rows([
            [4.0, 12.0, -16.0],
            [12.0, 37.0, -43.0],
            [-16.0, -43.0, 98.0],
        ]);

        let l = cholesky_decompose(&a).unwrap();

        // Verify L is lower triangular
        for i in 0..3 {
            for j in (i + 1)..3 {
                assert!(
                    l[(i, j)].abs() < 1e-10_f64,
                    "L[{},{}] = {}, should be zero",
                    i,
                    j,
                    l[(i, j)]
                );
            }
        }

        // Verify diagonal elements are positive
        for i in 0..3 {
            assert!(l[(i, i)] > 0.0, "L[{},{}] = {} should be positive", i, i, l[(i, i)]);
        }

        // Verify A = LL^T
        let lt = transpose(&l);
        let llt = multiply(&l, &lt);

        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    approximately_equal(a[(i, j)], llt[(i, j)], 1e-9, 1e-9),
                    "A[{},{}] = {}, LL^T[{},{}] = {}",
                    i,
                    j,
                    a[(i, j)],
                    i,
                    j,
                    llt[(i, j)]
                );
            }
        }
    }

    #[test]
    fn test_cholesky_decompose_2x2() {
        let a = Matrix::<f64, 2, 2>::from_rows([[4.0, 2.0], [2.0, 3.0]]);

        let l = cholesky_decompose(&a).unwrap();

        // Verify A = LL^T
        let lt = transpose(&l);
        let llt = multiply(&l, &lt);

        for i in 0..2 {
            for j in 0..2 {
                assert!(approximately_equal(a[(i, j)], llt[(i, j)], 1e-10, 1e-10));
            }
        }
    }

    #[test]
    fn test_cholesky_decompose_identity() {
        let a = Matrix::<f64, 3, 3>::identity();
        let l = cholesky_decompose(&a).unwrap();

        // For identity matrix, L should also be identity
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(approximately_equal(l[(i, j)], expected, 1e-10, 1e-10));
            }
        }
    }

    #[test]
    fn test_cholesky_decompose_not_symmetric() {
        // Not symmetric
        let a = Matrix::<f64, 2, 2>::from_rows([[4.0, 2.0], [3.0, 3.0]]);

        let result = cholesky_decompose(&a);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), LinearAlgebraError::NotSymmetric);
    }

    #[test]
    fn test_cholesky_decompose_not_positive_definite() {
        // Symmetric but not positive definite (has negative eigenvalue)
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [2.0, 1.0]]);

        let result = cholesky_decompose(&a);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), LinearAlgebraError::NotPositiveDefinite);
    }

    #[test]
    fn test_cholesky_solve_2x2() {
        let a = Matrix::<f64, 2, 2>::from_rows([[4.0, 2.0], [2.0, 3.0]]);

        let b = Vector::from([10.0, 9.0]);
        let x = cholesky_solve(&a, &b).unwrap();

        // Verify Ax = b
        let ax = crate::matrix::matrix_vector_multiply(&a, &x);
        for i in 0..2 {
            assert!(approximately_equal(ax[i], b[i], 1e-10, 1e-10));
        }
    }

    #[test]
    fn test_cholesky_solve_3x3() {
        let a = Matrix::<f64, 3, 3>::from_rows([
            [4.0, 12.0, -16.0],
            [12.0, 37.0, -43.0],
            [-16.0, -43.0, 98.0],
        ]);

        let b = Vector::from([1.0, 2.0, 3.0]);
        let x = cholesky_solve(&a, &b).unwrap();

        // Verify Ax = b
        let ax = crate::matrix::matrix_vector_multiply(&a, &x);
        for i in 0..3 {
            assert!(approximately_equal(ax[i], b[i], 1e-9, 1e-9));
        }
    }

    #[test]
    fn test_cholesky_diagonal_matrix() {
        // Diagonal positive definite matrix
        let a = Matrix::<f64, 3, 3>::from_rows([[4.0, 0.0, 0.0], [0.0, 9.0, 0.0], [0.0, 0.0, 16.0]]);

        let l = cholesky_decompose(&a).unwrap();

        // L should also be diagonal with sqrt of diagonal elements
        assert!(approximately_equal(l[(0, 0)], 2.0, 1e-10, 1e-10));
        assert!(approximately_equal(l[(1, 1)], 3.0, 1e-10, 1e-10));
        assert!(approximately_equal(l[(2, 2)], 4.0, 1e-10, 1e-10));

        // Off-diagonal should be zero
        for i in 0..3 {
            for j in 0..3 {
                if i != j {
                    assert!(l[(i, j)].abs() < 1e-10_f64);
                }
            }
        }
    }
}
