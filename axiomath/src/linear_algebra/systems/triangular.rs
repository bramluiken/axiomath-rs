//! Triangular system solvers.
//!
//! This module implements forward and backward substitution algorithms for solving
//! triangular systems of linear equations. These are fundamental building blocks
//! for more complex linear system solvers.
//!
//! # Algorithms
//!
//! - **Forward Substitution**: Solves Lx = b where L is lower triangular
//! - **Backward Substitution**: Solves Ux = b where U is upper triangular
//!
//! # Complexity
//!
//! Both algorithms have O(n²) time complexity and O(n) space complexity.
//!
//! # References
//!
//! - Golub & Van Loan, "Matrix Computations", Section 3.1
//! - Trefethen & Bau, "Numerical Linear Algebra", Lecture 17

use crate::core::traits::Float;
use crate::linear_algebra::LinearAlgebraError;
use crate::matrix::Matrix;
use crate::vector::Vector;
use num_traits::Zero;

/// Solves a lower triangular system Lx = b using forward substitution.
///
/// # Mathematical Definition
///
/// For a lower triangular matrix L where L_ij = 0 for j > i, we solve:
///
/// ```text
/// L₁₁ x₁                                  = b₁
/// L₂₁ x₁ + L₂₂ x₂                        = b₂
/// L₃₁ x₁ + L₃₂ x₂ + L₃₃ x₃              = b₃
/// ...
/// ```
///
/// # Algorithm
///
/// ```text
/// for i = 1 to n:
///     sum = b_i
///     for j = 1 to i-1:
///         sum -= L_ij * x_j
///     x_i = sum / L_ii
/// ```
///
/// # Numerical Stability
///
/// Forward substitution is numerically stable. The only potential issue is
/// division by zero if any diagonal element is zero (matrix is singular).
///
/// # Complexity
///
/// - Time: O(n²)
/// - Space: O(n)
///
/// # Arguments
///
/// * `l` - Lower triangular matrix (n×n)
/// * `b` - Right-hand side vector (n)
///
/// # Returns
///
/// Solution vector x such that Lx = b.
///
/// Returns `Err(LinearAlgebraError::SingularMatrix)` if any diagonal element
/// of L is zero or numerically zero.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::vector::Vector;
/// use axiomath::linear_algebra::systems::solve_lower_triangular;
///
/// let l = Matrix::from_rows([
///     [2.0, 0.0, 0.0],
///     [1.0, 3.0, 0.0],
///     [4.0, 2.0, 5.0],
/// ]);
///
/// let b = Vector::from([4.0, 10.0, 31.0]);
/// let x = solve_lower_triangular(&l, &b).unwrap();
///
/// // x should be [2.0, 8/3, 3.0]
/// ```
///
/// # Panics
///
/// Panics if the matrix is not square or if dimensions don't match.
pub fn solve_lower_triangular<T: Float, const N: usize>(
    l: &Matrix<T, N, N>,
    b: &Vector<T, N>,
) -> Result<Vector<T, N>, LinearAlgebraError> {
    let epsilon = T::epsilon() * T::from(N).unwrap() * T::from(100.0).unwrap();
    let mut x = Vector::<T, N>::zero();

    for i in 0..N {
        // Check for zero diagonal (singular matrix)
        if l[(i, i)].abs() < epsilon {
            return Err(LinearAlgebraError::SingularMatrix);
        }

        // Compute: x_i = (b_i - sum(L_ij * x_j for j < i)) / L_ii
        let mut sum = b[i];
        for j in 0..i {
            sum = sum - l[(i, j)] * x[j];
        }
        x[i] = sum / l[(i, i)];
    }

    Ok(x)
}

/// Solves an upper triangular system Ux = b using backward substitution.
///
/// # Mathematical Definition
///
/// For an upper triangular matrix U where U_ij = 0 for i > j, we solve:
///
/// ```text
/// U₁₁ x₁ + U₁₂ x₂ + U₁₃ x₃ + ... + U₁ₙ xₙ = b₁
///          U₂₂ x₂ + U₂₃ x₃ + ... + U₂ₙ xₙ = b₂
///                   U₃₃ x₃ + ... + U₃ₙ xₙ = b₃
///                                      ...
///                                 Uₙₙ xₙ = bₙ
/// ```
///
/// # Algorithm
///
/// ```text
/// for i = n down to 1:
///     sum = b_i
///     for j = i+1 to n:
///         sum -= U_ij * x_j
///     x_i = sum / U_ii
/// ```
///
/// # Numerical Stability
///
/// Backward substitution is numerically stable. The only potential issue is
/// division by zero if any diagonal element is zero (matrix is singular).
///
/// # Complexity
///
/// - Time: O(n²)
/// - Space: O(n)
///
/// # Arguments
///
/// * `u` - Upper triangular matrix (n×n)
/// * `b` - Right-hand side vector (n)
///
/// # Returns
///
/// Solution vector x such that Ux = b.
///
/// Returns `Err(LinearAlgebraError::SingularMatrix)` if any diagonal element
/// of U is zero or numerically zero.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::vector::Vector;
/// use axiomath::linear_algebra::systems::solve_upper_triangular;
///
/// let u = Matrix::from_rows([
///     [2.0, 1.0, 3.0],
///     [0.0, 4.0, 2.0],
///     [0.0, 0.0, 5.0],
/// ]);
///
/// let b = Vector::from([17.0, 18.0, 10.0]);
/// let x = solve_upper_triangular(&u, &b).unwrap();
///
/// // x should be [2.0, 3.0, 2.0]
/// ```
///
/// # Panics
///
/// Panics if the matrix is not square or if dimensions don't match.
pub fn solve_upper_triangular<T: Float, const N: usize>(
    u: &Matrix<T, N, N>,
    b: &Vector<T, N>,
) -> Result<Vector<T, N>, LinearAlgebraError> {
    let epsilon = T::epsilon() * T::from(N).unwrap() * T::from(100.0).unwrap();
    let mut x = Vector::<T, N>::zero();

    // Iterate backwards from last row to first
    for i in (0..N).rev() {
        // Check for zero diagonal (singular matrix)
        if u[(i, i)].abs() < epsilon {
            return Err(LinearAlgebraError::SingularMatrix);
        }

        // Compute: x_i = (b_i - sum(U_ij * x_j for j > i)) / U_ii
        let mut sum = b[i];
        for j in (i + 1)..N {
            sum = sum - u[(i, j)] * x[j];
        }
        x[i] = sum / u[(i, i)];
    }

    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::matrix_vector_multiply;
    use crate::utils::comparison::approximately_equal;

    #[test]
    fn test_solve_lower_triangular_3x3() {
        let l = Matrix::from_rows([[2.0, 0.0, 0.0], [1.0, 3.0, 0.0], [4.0, 2.0, 5.0]]);

        let b = Vector::from([4.0, 10.0, 31.0]);
        let x = solve_lower_triangular(&l, &b).unwrap();

        // Verify Lx = b
        let lx = matrix_vector_multiply(&l, &x);
        for i in 0..3 {
            assert!(
                approximately_equal(lx[i], b[i], 1e-10, 1e-10),
                "Lx[{}] = {}, b[{}] = {}",
                i,
                lx[i],
                i,
                b[i]
            );
        }
    }

    #[test]
    fn test_solve_upper_triangular_3x3() {
        let u = Matrix::from_rows([[2.0, 1.0, 3.0], [0.0, 4.0, 2.0], [0.0, 0.0, 5.0]]);

        let b = Vector::from([17.0, 18.0, 10.0]);
        let x = solve_upper_triangular(&u, &b).unwrap();

        // Verify Ux = b
        let ux = matrix_vector_multiply(&u, &x);
        for i in 0..3 {
            assert!(
                approximately_equal(ux[i], b[i], 1e-10, 1e-10),
                "Ux[{}] = {}, b[{}] = {}",
                i,
                ux[i],
                i,
                b[i]
            );
        }
    }

    #[test]
    fn test_solve_lower_triangular_identity() {
        let l = Matrix::<f64, 3, 3>::identity();
        let b = Vector::from([1.0, 2.0, 3.0]);
        let x = solve_lower_triangular(&l, &b).unwrap();

        // For identity matrix, x should equal b
        for i in 0..3 {
            assert!((x[i] - b[i]).abs() < 1e-10);
        }
    }

    #[test]
    fn test_solve_upper_triangular_identity() {
        let u = Matrix::<f64, 3, 3>::identity();
        let b = Vector::from([1.0, 2.0, 3.0]);
        let x = solve_upper_triangular(&u, &b).unwrap();

        // For identity matrix, x should equal b
        for i in 0..3 {
            assert!((x[i] - b[i]).abs() < 1e-10);
        }
    }

    #[test]
    fn test_solve_lower_triangular_singular() {
        let l = Matrix::from_rows([[2.0, 0.0, 0.0], [1.0, 0.0, 0.0], [4.0, 2.0, 5.0]]);

        let b = Vector::from([4.0, 10.0, 31.0]);
        let result = solve_lower_triangular(&l, &b);

        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), LinearAlgebraError::SingularMatrix);
    }

    #[test]
    fn test_solve_upper_triangular_singular() {
        let u = Matrix::from_rows([[2.0, 1.0, 3.0], [0.0, 4.0, 2.0], [0.0, 0.0, 0.0]]);

        let b = Vector::from([17.0, 18.0, 10.0]);
        let result = solve_upper_triangular(&u, &b);

        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), LinearAlgebraError::SingularMatrix);
    }

    #[test]
    fn test_solve_lower_triangular_unit_diagonal() {
        // Lower triangular with ones on diagonal (common in LU decomposition)
        let l = Matrix::from_rows([[1.0, 0.0, 0.0], [0.5, 1.0, 0.0], [2.0, 1.0, 1.0]]);

        let b = Vector::from([3.0, 4.0, 10.0]);
        let x = solve_lower_triangular(&l, &b).unwrap();

        // Verify Lx = b
        let lx = matrix_vector_multiply(&l, &x);
        for i in 0..3 {
            assert!(approximately_equal(lx[i], b[i], 1e-10, 1e-10));
        }
    }

    #[test]
    fn test_solve_upper_triangular_diagonal() {
        // Diagonal matrix (special case of upper triangular)
        let u = Matrix::from_rows([[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]]);

        let b = Vector::from([6.0, 9.0, 12.0]);
        let x = solve_upper_triangular(&u, &b).unwrap();

        // x should be [3.0, 3.0, 3.0]
        assert!((x[0] - 3.0_f64).abs() < 1e-10);
        assert!((x[1] - 3.0_f64).abs() < 1e-10);
        assert!((x[2] - 3.0_f64).abs() < 1e-10);
    }
}
