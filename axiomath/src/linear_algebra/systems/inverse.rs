//! Matrix inversion algorithms.
//!
//! This module implements matrix inversion using LU decomposition and
//! the Moore-Penrose pseudoinverse via SVD.
//!
//! # Methods
//!
//! - **Matrix Inversion**: For square non-singular matrices via LU decomposition
//! - **Pseudoinverse**: For rectangular or singular matrices via SVD (not yet implemented)
//!
//! # Mathematical Background
//!
//! For a square invertible matrix A, the inverse A⁻¹ satisfies:
//! ```text
//! A A⁻¹ = A⁻¹ A = I
//! ```
//!
//! # Algorithms
//!
//! Matrix inversion via LU decomposition:
//! 1. Compute PA = LU
//! 2. For each column i of the identity matrix:
//!    - Solve Ly = Pe_i (forward substitution)
//!    - Solve Ux = y (backward substitution)
//!    - Column i of A⁻¹ is x
//!
//! # Complexity
//!
//! - Time: O(n³) for LU decomposition + n × O(n²) for solving = O(n³)
//! - Space: O(n²)
//!
//! # Numerical Stability
//!
//! Matrix inversion is numerically sensitive. The condition number κ(A) = ||A|| ||A⁻¹||
//! measures sensitivity to errors. For ill-conditioned matrices (large κ), consider
//! using SVD-based pseudoinverse or solving Ax = b directly without computing A⁻¹.
//!
//! # References
//!
//! - Golub & Van Loan, "Matrix Computations", Section 3.2
//! - Trefethen & Bau, "Numerical Linear Algebra", Lecture 14

use crate::core::traits::Float;
use crate::linear_algebra::decomposition::lu_decompose;
use crate::linear_algebra::systems::triangular::{solve_lower_triangular, solve_upper_triangular};
use crate::linear_algebra::LinearAlgebraError;
use crate::matrix::Matrix;
use crate::vector::Vector;
use num_traits::Zero;

/// Computes the inverse of a square matrix using LU decomposition.
///
/// # Mathematical Definition
///
/// For an invertible matrix A, computes A⁻¹ such that:
/// ```text
/// A A⁻¹ = A⁻¹ A = I
/// ```
///
/// # Algorithm
///
/// 1. Compute LU decomposition with pivoting: PA = LU
/// 2. For each standard basis vector e_i:
///    - Apply permutation: b = Pe_i
///    - Solve Ly = b via forward substitution
///    - Solve Ux = y via backward substitution
///    - Column i of A⁻¹ is x
///
/// # Complexity
///
/// - Time: O(n³) for LU + n × O(n²) for solving = O(n³)
/// - Space: O(n²)
///
/// # Numerical Considerations
///
/// - For ill-conditioned matrices, the computed inverse may be inaccurate
/// - If you need to solve Ax = b, prefer direct solving over computing A⁻¹
/// - Use condition number estimation to assess numerical stability
///
/// # Arguments
///
/// * `matrix` - Square matrix to invert (n×n)
///
/// # Returns
///
/// Returns `Ok(A⁻¹)` where A⁻¹ is the matrix inverse.
///
/// Returns `Err(LinearAlgebraError::SingularMatrix)` if the matrix is
/// singular or numerically singular.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::systems::invert;
///
/// let a = Matrix::from_rows([
///     [4.0, 7.0],
///     [2.0, 6.0],
/// ]);
///
/// let a_inv = invert(&a).unwrap();
///
/// // Verify: A A⁻¹ ≈ I
/// ```
pub fn invert<T: Float, const N: usize>(
    matrix: &Matrix<T, N, N>,
) -> Result<Matrix<T, N, N>, LinearAlgebraError> {
    // Step 1: Compute LU decomposition with pivoting
    let (l, u, perm) = lu_decompose(matrix)?;

    // Step 2: Solve for each column of the inverse
    let mut inverse = Matrix::<T, N, N>::zero();

    for i in 0..N {
        // Create standard basis vector e_i
        let mut e_i = Vector::<T, N>::zero();
        e_i[i] = T::one();

        // Apply permutation: b = Pe_i
        let mut b = Vector::<T, N>::zero();
        for j in 0..N {
            b[j] = e_i[perm[j]];
        }

        // Solve Ly = b
        let y = solve_lower_triangular(&l, &b)?;

        // Solve Ux = y
        let x = solve_upper_triangular(&u, &y)?;

        // Store x as column i of the inverse
        for j in 0..N {
            inverse[(j, i)] = x[j];
        }
    }

    Ok(inverse)
}

/// Moore-Penrose pseudoinverse via SVD.
///
/// The pseudoinverse A⁺ generalizes matrix inversion to rectangular and
/// singular matrices. For a matrix A (m×n), A⁺ (n×m) satisfies:
///
/// ```text
/// A A⁺ A = A
/// A⁺ A A⁺ = A⁺
/// (A A⁺)ᵀ = A A⁺
/// (A⁺ A)ᵀ = A⁺ A
/// ```
///
/// # Algorithm
///
/// Via SVD: A = UΣVᵀ, then A⁺ = VΣ⁺Uᵀ where Σ⁺ is obtained by:
/// - Transposing Σ
/// - Taking reciprocal of non-zero singular values
/// - Setting very small singular values to zero (for numerical stability)
///
/// # Arguments
///
/// * `matrix` - Matrix to compute pseudoinverse for (m×n)
///
/// # Returns
///
/// Pseudoinverse A⁺ (n×m)
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::systems::pseudoinverse;
///
/// let a = Matrix::from_rows([
///     [1.0, 2.0],
///     [3.0, 4.0],
/// ]);
///
/// let a_pinv = pseudoinverse(&a).unwrap();
/// // For invertible matrices, pseudoinverse equals the inverse
/// ```
pub fn pseudoinverse<T: Float, const M: usize, const N: usize>(
    matrix: &Matrix<T, M, N>,
) -> Result<Matrix<T, N, M>, LinearAlgebraError> {
    use crate::linear_algebra::decomposition::svd::pseudoinverse as svd_pseudoinverse;

    // Use a tolerance based on machine epsilon and matrix size
    let tolerance = T::epsilon() * T::from(M.max(N)).unwrap() * T::from(10.0).unwrap();
    svd_pseudoinverse(matrix, tolerance)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::multiply;
    use crate::utils::comparison::approximately_equal;

    #[test]
    fn test_invert_2x2() {
        let a = Matrix::<f64, 2, 2>::from_rows([[4.0, 7.0], [2.0, 6.0]]);

        let a_inv = invert(&a).unwrap();

        // Verify A A⁻¹ = I
        let product = multiply(&a, &a_inv);
        for i in 0..2 {
            for j in 0..2 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    approximately_equal(product[(i, j)], expected, 1e-10, 1e-10),
                    "A A⁻¹[{},{}] = {}, expected {}",
                    i,
                    j,
                    product[(i, j)],
                    expected
                );
            }
        }

        // Verify A⁻¹ A = I
        let product2 = multiply(&a_inv, &a);
        for i in 0..2 {
            for j in 0..2 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    approximately_equal(product2[(i, j)], expected, 1e-10, 1e-10),
                    "A⁻¹ A[{},{}] = {}, expected {}",
                    i,
                    j,
                    product2[(i, j)],
                    expected
                );
            }
        }
    }

    #[test]
    fn test_invert_3x3() {
        let a = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [0.0, 1.0, 4.0],
            [5.0, 6.0, 0.0],
        ]);

        let a_inv = invert(&a).unwrap();

        // Verify A A⁻¹ = I
        let product = multiply(&a, &a_inv);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    approximately_equal(product[(i, j)], expected, 1e-10, 1e-10),
                    "A A⁻¹[{},{}] = {}, expected {}",
                    i,
                    j,
                    product[(i, j)],
                    expected
                );
            }
        }
    }

    #[test]
    fn test_invert_identity() {
        let a = Matrix::<f64, 3, 3>::identity();
        let a_inv = invert(&a).unwrap();

        // Inverse of identity is identity
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(approximately_equal(a_inv[(i, j)], expected, 1e-10, 1e-10));
            }
        }
    }

    #[test]
    fn test_invert_diagonal() {
        let a = Matrix::<f64, 3, 3>::from_diagonal([2.0, 3.0, 4.0]);
        let a_inv = invert(&a).unwrap();

        // Inverse of diagonal is diagonal with reciprocals
        assert!(approximately_equal(a_inv[(0, 0)], 0.5, 1e-10, 1e-10));
        assert!(approximately_equal(a_inv[(1, 1)], 1.0 / 3.0, 1e-10, 1e-10));
        assert!(approximately_equal(a_inv[(2, 2)], 0.25, 1e-10, 1e-10));

        // Off-diagonal should be zero
        for i in 0..3 {
            for j in 0..3 {
                if i != j {
                    assert!(a_inv[(i, j)].abs() < 1e-10_f64);
                }
            }
        }
    }

    #[test]
    fn test_invert_singular() {
        // Singular matrix (determinant = 0)
        let a = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [2.0, 4.0, 6.0],
            [1.0, 1.0, 1.0],
        ]);

        let result = invert(&a);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), LinearAlgebraError::SingularMatrix);
    }

    #[test]
    fn test_invert_with_pivoting() {
        // Matrix that requires pivoting during LU decomposition
        let a = Matrix::<f64, 2, 2>::from_rows([[0.0, 1.0], [1.0, 1.0]]);

        let a_inv = invert(&a).unwrap();

        // Verify A A⁻¹ = I
        let product = multiply(&a, &a_inv);
        for i in 0..2 {
            for j in 0..2 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(approximately_equal(product[(i, j)], expected, 1e-10, 1e-10));
            }
        }
    }

    #[test]
    fn test_invert_consistency() {
        // Test that (A⁻¹)⁻¹ = A
        let a = Matrix::<f64, 2, 2>::from_rows([[3.0, 2.0], [1.0, 4.0]]);

        let a_inv = invert(&a).unwrap();
        let a_inv_inv = invert(&a_inv).unwrap();

        // (A⁻¹)⁻¹ should equal A
        for i in 0..2 {
            for j in 0..2 {
                assert!(
                    approximately_equal(a_inv_inv[(i, j)], a[(i, j)], 1e-9, 1e-9),
                    "(A⁻¹)⁻¹[{},{}] = {}, A[{},{}] = {}",
                    i,
                    j,
                    a_inv_inv[(i, j)],
                    i,
                    j,
                    a[(i, j)]
                );
            }
        }
    }

    #[test]
    fn test_invert_determinant_property() {
        // Test that det(A⁻¹) = 1/det(A)
        use crate::matrix::determinant;

        let a = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 2.0, 1.0],
            [0.0, 1.0, 3.0],
            [2.0, 1.0, 1.0],
        ]);

        let a_inv = invert(&a).unwrap();

        let det_a = determinant(&a);
        let det_a_inv = determinant(&a_inv);

        // det(A⁻¹) * det(A) should equal 1
        assert!(
            approximately_equal(det_a * det_a_inv, 1.0, 1e-9, 1e-9),
            "det(A) * det(A⁻¹) = {} * {} = {}, expected 1.0",
            det_a,
            det_a_inv,
            det_a * det_a_inv
        );
    }
}
