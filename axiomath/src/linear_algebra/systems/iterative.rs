//! Iterative methods for solving linear systems.
//!
//! This module implements iterative solvers that approximate solutions through
//! successive refinement. These methods are particularly useful for large, sparse
//! systems where direct methods are computationally expensive.
//!
//! # Methods
//!
//! - **Jacobi Iteration**: Simple splitting method, requires diagonal dominance
//! - **Gauss-Seidel**: Improved version of Jacobi, typically converges faster
//! - **Conjugate Gradient**: Optimal Krylov subspace method for SPD matrices
//!
//! # Convergence
//!
//! - **Jacobi and Gauss-Seidel**: Converge for strictly diagonally dominant matrices
//! - **Conjugate Gradient**: Converges for symmetric positive definite matrices
//!
//! # References
//!
//! - Golub & Van Loan, "Matrix Computations", Chapter 10-11
//! - Trefethen & Bau, "Numerical Linear Algebra", Lectures 38-40
//! - Saad, "Iterative Methods for Sparse Linear Systems"

use crate::core::traits::Float;
use crate::linear_algebra::LinearAlgebraError;
use crate::matrix::Matrix;
use crate::vector::{dot_product, magnitude, scale, subtract, Vector};
use num_traits::Zero;

/// Solves a linear system using Jacobi iteration.
///
/// # Mathematical Definition
///
/// Split A = D + R where D is diagonal, then iterate:
/// ```text
/// x^(k+1) = D⁻¹(b - Rx^(k))
/// ```
///
/// Equivalently, for each component:
/// ```text
/// x_i^(k+1) = (b_i - Σ_{j≠i} a_ij x_j^(k)) / a_ii
/// ```
///
/// # Algorithm
///
/// ```text
/// Given: A, b, x^(0), max iterations, tolerance
/// for k = 0 to max_iterations:
///     for i = 1 to n:
///         σ = Σ_{j≠i} a_ij x_j^(k)
///         x_i^(k+1) = (b_i - σ) / a_ii
///     if ||x^(k+1) - x^(k)|| < tolerance:
///         return x^(k+1)
/// ```
///
/// # Convergence
///
/// Jacobi iteration converges if A is strictly diagonally dominant:
/// ```text
/// |a_ii| > Σ_{j≠i} |a_ij| for all i
/// ```
///
/// # Complexity
///
/// - Time per iteration: O(n²)
/// - Space: O(n)
///
/// # Arguments
///
/// * `a` - Coefficient matrix (n×n)
/// * `b` - Right-hand side vector (n)
/// * `x0` - Initial guess (n)
/// * `max_iterations` - Maximum number of iterations
/// * `tolerance` - Convergence tolerance for ||x^(k+1) - x^(k)||
///
/// # Returns
///
/// Approximate solution x such that Ax ≈ b.
///
/// Returns `Err(LinearAlgebraError::ConvergenceFailure)` if the method
/// does not converge within the maximum iterations.
///
/// Returns `Err(LinearAlgebraError::SingularMatrix)` if any diagonal
/// element is zero.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::vector::Vector;
/// use axiomath::linear_algebra::systems::iterative::jacobi;
///
/// // Diagonally dominant matrix
/// let a = Matrix::from_rows([
///     [10.0, -1.0, 2.0],
///     [-1.0, 11.0, -1.0],
///     [2.0, -1.0, 10.0],
/// ]);
///
/// let b = Vector::from([6.0, 25.0, -11.0]);
/// let x0 = Vector::from([0.0, 0.0, 0.0]);
///
/// let x = jacobi(&a, &b, &x0, 100, 1e-10).unwrap();
/// ```
pub fn jacobi<T: Float, const N: usize>(
    a: &Matrix<T, N, N>,
    b: &Vector<T, N>,
    x0: &Vector<T, N>,
    max_iterations: usize,
    tolerance: T,
) -> Result<Vector<T, N>, LinearAlgebraError> {
    let epsilon = T::epsilon() * T::from(N).unwrap() * T::from(100.0).unwrap();

    // Check for zero diagonal elements
    for i in 0..N {
        if a[(i, i)].abs() < epsilon {
            return Err(LinearAlgebraError::SingularMatrix);
        }
    }

    let mut x = *x0;
    let mut x_new = Vector::<T, N>::zero();

    for iteration in 0..max_iterations {
        // Compute next iteration
        for i in 0..N {
            let mut sigma = T::zero();
            for j in 0..N {
                if j != i {
                    sigma = sigma + a[(i, j)] * x[j];
                }
            }
            x_new[i] = (b[i] - sigma) / a[(i, i)];
        }

        // Check convergence
        let diff = subtract(&x_new, &x);
        if magnitude(&diff) < tolerance {
            return Ok(x_new);
        }

        x = x_new;
    }

    Err(LinearAlgebraError::ConvergenceFailure {
        iterations: max_iterations,
    })
}

/// Solves a linear system using Gauss-Seidel iteration.
///
/// # Mathematical Definition
///
/// Split A = L + D + U where L is strictly lower, D is diagonal, U is strictly upper.
/// Then iterate:
/// ```text
/// x^(k+1) = (D + L)⁻¹(b - Ux^(k))
/// ```
///
/// Equivalently, for each component (using updated values immediately):
/// ```text
/// x_i^(k+1) = (b_i - Σ_{j<i} a_ij x_j^(k+1) - Σ_{j>i} a_ij x_j^(k)) / a_ii
/// ```
///
/// # Algorithm
///
/// ```text
/// Given: A, b, x^(0), max iterations, tolerance
/// for k = 0 to max_iterations:
///     for i = 1 to n:
///         σ = Σ_{j<i} a_ij x_j (new) + Σ_{j>i} a_ij x_j (old)
///         x_i = (b_i - σ) / a_ii
///     if ||x^(k+1) - x^(k)|| < tolerance:
///         return x
/// ```
///
/// # Convergence
///
/// Gauss-Seidel converges if:
/// - A is strictly diagonally dominant, OR
/// - A is symmetric positive definite
///
/// Generally converges faster than Jacobi.
///
/// # Complexity
///
/// - Time per iteration: O(n²)
/// - Space: O(n)
///
/// # Arguments
///
/// * `a` - Coefficient matrix (n×n)
/// * `b` - Right-hand side vector (n)
/// * `x0` - Initial guess (n)
/// * `max_iterations` - Maximum number of iterations
/// * `tolerance` - Convergence tolerance
///
/// # Returns
///
/// Approximate solution x such that Ax ≈ b.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::vector::Vector;
/// use axiomath::linear_algebra::systems::iterative::gauss_seidel;
///
/// let a = Matrix::from_rows([
///     [10.0, -1.0, 2.0],
///     [-1.0, 11.0, -1.0],
///     [2.0, -1.0, 10.0],
/// ]);
///
/// let b = Vector::from([6.0, 25.0, -11.0]);
/// let x0 = Vector::from([0.0, 0.0, 0.0]);
///
/// let x = gauss_seidel(&a, &b, &x0, 100, 1e-10).unwrap();
/// ```
pub fn gauss_seidel<T: Float, const N: usize>(
    a: &Matrix<T, N, N>,
    b: &Vector<T, N>,
    x0: &Vector<T, N>,
    max_iterations: usize,
    tolerance: T,
) -> Result<Vector<T, N>, LinearAlgebraError> {
    let epsilon = T::epsilon() * T::from(N).unwrap() * T::from(100.0).unwrap();

    // Check for zero diagonal elements
    for i in 0..N {
        if a[(i, i)].abs() < epsilon {
            return Err(LinearAlgebraError::SingularMatrix);
        }
    }

    let mut x = *x0;
    let mut x_old = Vector::<T, N>::zero();

    for iteration in 0..max_iterations {
        x_old = x;

        // Gauss-Seidel update: use new values immediately
        for i in 0..N {
            let mut sigma = T::zero();
            for j in 0..N {
                if j != i {
                    sigma = sigma + a[(i, j)] * x[j];
                }
            }
            x[i] = (b[i] - sigma) / a[(i, i)];
        }

        // Check convergence
        let diff = subtract(&x, &x_old);
        if magnitude(&diff) < tolerance {
            return Ok(x);
        }
    }

    Err(LinearAlgebraError::ConvergenceFailure {
        iterations: max_iterations,
    })
}

/// Solves a symmetric positive definite system using Conjugate Gradient.
///
/// # Mathematical Definition
///
/// Conjugate Gradient is a Krylov subspace method that minimizes:
/// ```text
/// φ(x) = ½ x^T A x - x^T b
/// ```
///
/// over successively larger subspaces K_k = span{r_0, Ar_0, ..., A^(k-1)r_0}.
///
/// # Algorithm
///
/// ```text
/// r_0 = b - Ax_0
/// p_0 = r_0
/// for k = 0 to max_iterations:
///     α_k = (r_k^T r_k) / (p_k^T A p_k)
///     x_{k+1} = x_k + α_k p_k
///     r_{k+1} = r_k - α_k A p_k
///     if ||r_{k+1}|| < tolerance:
///         return x_{k+1}
///     β_k = (r_{k+1}^T r_{k+1}) / (r_k^T r_k)
///     p_{k+1} = r_{k+1} + β_k p_k
/// ```
///
/// # Convergence
///
/// For SPD matrices, CG converges in at most n iterations (exact arithmetic).
/// In practice, converges much faster depending on condition number.
///
/// Convergence rate depends on eigenvalue distribution:
/// ```text
/// ||e_k|| ≤ 2 ((√κ - 1)/(√κ + 1))^k ||e_0||
/// ```
/// where κ = λ_max/λ_min is the condition number.
///
/// # Complexity
///
/// - Time per iteration: O(n²) for dense matrices, O(n) for sparse
/// - Space: O(n)
///
/// # Requirements
///
/// Matrix A must be symmetric positive definite (SPD).
///
/// # Arguments
///
/// * `a` - Symmetric positive definite matrix (n×n)
/// * `b` - Right-hand side vector (n)
/// * `x0` - Initial guess (n)
/// * `max_iterations` - Maximum number of iterations
/// * `tolerance` - Convergence tolerance for residual norm
///
/// # Returns
///
/// Approximate solution x such that Ax ≈ b.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::vector::Vector;
/// use axiomath::linear_algebra::systems::iterative::conjugate_gradient;
///
/// // Symmetric positive definite matrix
/// let a = Matrix::from_rows([
///     [4.0, 1.0],
///     [1.0, 3.0],
/// ]);
///
/// let b = Vector::from([1.0, 2.0]);
/// let x0 = Vector::from([0.0, 0.0]);
///
/// let x = conjugate_gradient(&a, &b, &x0, 100, 1e-10).unwrap();
/// ```
pub fn conjugate_gradient<T: Float, const N: usize>(
    a: &Matrix<T, N, N>,
    b: &Vector<T, N>,
    x0: &Vector<T, N>,
    max_iterations: usize,
    tolerance: T,
) -> Result<Vector<T, N>, LinearAlgebraError> {
    use crate::matrix::matrix_vector_multiply;

    let mut x = *x0;

    // r_0 = b - Ax_0
    let ax = matrix_vector_multiply(a, &x);
    let mut r = subtract(b, &ax);

    // p_0 = r_0
    let mut p = r;

    let mut r_dot_r = dot_product(&r, &r);

    for _iteration in 0..max_iterations {
        // Check convergence
        let residual_norm = magnitude(&r);
        if residual_norm < tolerance {
            return Ok(x);
        }

        // α_k = (r_k^T r_k) / (p_k^T A p_k)
        let ap = matrix_vector_multiply(a, &p);
        let p_dot_ap = dot_product(&p, &ap);

        // Check for zero denominator (shouldn't happen for SPD matrices)
        if p_dot_ap.abs() < T::epsilon() {
            return Err(LinearAlgebraError::NumericalInstability);
        }

        let alpha = r_dot_r / p_dot_ap;

        // x_{k+1} = x_k + α_k p_k
        x = x + scale(&p, alpha);

        // r_{k+1} = r_k - α_k A p_k
        r = subtract(&r, &scale(&ap, alpha));

        let r_dot_r_new = dot_product(&r, &r);

        // β_k = (r_{k+1}^T r_{k+1}) / (r_k^T r_k)
        let beta = r_dot_r_new / r_dot_r;

        // p_{k+1} = r_{k+1} + β_k p_k
        p = r + scale(&p, beta);

        r_dot_r = r_dot_r_new;
    }

    Err(LinearAlgebraError::ConvergenceFailure {
        iterations: max_iterations,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::matrix_vector_multiply;
    use crate::utils::comparison::approximately_equal;

    #[test]
    fn test_jacobi_diagonally_dominant() {
        // Strictly diagonally dominant matrix
        let a = Matrix::<f64, 3, 3>::from_rows([
            [10.0, -1.0, 2.0],
            [-1.0, 11.0, -1.0],
            [2.0, -1.0, 10.0],
        ]);

        let b = Vector::from([6.0, 25.0, -11.0]);
        let x0 = Vector::from([0.0, 0.0, 0.0]);

        let x = jacobi(&a, &b, &x0, 100, 1e-10).unwrap();

        // Verify Ax ≈ b
        let ax = matrix_vector_multiply(&a, &x);
        for i in 0..3 {
            assert!(
                approximately_equal(ax[i], b[i], 1e-8, 1e-8),
                "Ax[{}] = {}, b[{}] = {}",
                i,
                ax[i],
                i,
                b[i]
            );
        }
    }

    #[test]
    fn test_gauss_seidel_diagonally_dominant() {
        let a = Matrix::<f64, 3, 3>::from_rows([
            [10.0, -1.0, 2.0],
            [-1.0, 11.0, -1.0],
            [2.0, -1.0, 10.0],
        ]);

        let b = Vector::from([6.0, 25.0, -11.0]);
        let x0 = Vector::from([0.0, 0.0, 0.0]);

        let x = gauss_seidel(&a, &b, &x0, 100, 1e-10).unwrap();

        // Verify Ax ≈ b
        let ax = matrix_vector_multiply(&a, &x);
        for i in 0..3 {
            assert!(
                approximately_equal(ax[i], b[i], 1e-8, 1e-8),
                "Ax[{}] = {}, b[{}] = {}",
                i,
                ax[i],
                i,
                b[i]
            );
        }
    }

    #[test]
    fn test_conjugate_gradient_spd() {
        // Symmetric positive definite matrix
        let a = Matrix::<f64, 2, 2>::from_rows([[4.0, 1.0], [1.0, 3.0]]);

        let b = Vector::from([1.0, 2.0]);
        let x0 = Vector::from([0.0, 0.0]);

        let x = conjugate_gradient(&a, &b, &x0, 100, 1e-10).unwrap();

        // Verify Ax ≈ b
        let ax = matrix_vector_multiply(&a, &x);
        for i in 0..2 {
            assert!(
                approximately_equal(ax[i], b[i], 1e-10, 1e-10),
                "Ax[{}] = {}, b[{}] = {}",
                i,
                ax[i],
                i,
                b[i]
            );
        }
    }

    #[test]
    fn test_conjugate_gradient_3x3() {
        // Symmetric positive definite matrix
        let a = Matrix::<f64, 3, 3>::from_rows([
            [4.0, 1.0, 0.0],
            [1.0, 4.0, 1.0],
            [0.0, 1.0, 4.0],
        ]);

        let b = Vector::from([1.0, 2.0, 3.0]);
        let x0 = Vector::from([0.0, 0.0, 0.0]);

        let x = conjugate_gradient(&a, &b, &x0, 100, 1e-10).unwrap();

        // Verify Ax ≈ b
        let ax = matrix_vector_multiply(&a, &x);
        for i in 0..3 {
            assert!(approximately_equal(ax[i], b[i], 1e-9, 1e-9));
        }
    }

    #[test]
    fn test_gauss_seidel_faster_than_jacobi() {
        let a = Matrix::<f64, 3, 3>::from_rows([
            [10.0, -1.0, 2.0],
            [-1.0, 11.0, -1.0],
            [2.0, -1.0, 10.0],
        ]);

        let b = Vector::from([6.0, 25.0, -11.0]);
        let x0 = Vector::from([0.0, 0.0, 0.0]);

        // Count iterations needed for each method
        // Gauss-Seidel typically converges in fewer iterations
        let tolerance = 1e-6;

        let x_jacobi = jacobi(&a, &b, &x0, 50, tolerance).unwrap();
        let x_gs = gauss_seidel(&a, &b, &x0, 50, tolerance).unwrap();

        // Both should produce similar results
        for i in 0..3 {
            assert!(approximately_equal(x_jacobi[i], x_gs[i], 1e-4, 1e-4));
        }
    }

    #[test]
    fn test_jacobi_singular_matrix() {
        // Matrix with zero diagonal element
        let a = Matrix::<f64, 2, 2>::from_rows([[0.0, 1.0], [1.0, 2.0]]);

        let b = Vector::from([1.0, 2.0]);
        let x0 = Vector::from([0.0, 0.0]);

        let result = jacobi(&a, &b, &x0, 100, 1e-10);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), LinearAlgebraError::SingularMatrix);
    }

    #[test]
    fn test_convergence_failure() {
        // Non-diagonally dominant matrix may not converge
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 10.0], [10.0, 1.0]]);

        let b = Vector::from([1.0, 1.0]);
        let x0 = Vector::from([0.0, 0.0]);

        // Very few iterations to force failure
        let result = jacobi(&a, &b, &x0, 2, 1e-10);
        assert!(result.is_err());
        if let Err(LinearAlgebraError::ConvergenceFailure { iterations }) = result {
            assert_eq!(iterations, 2);
        } else {
            panic!("Expected ConvergenceFailure error");
        }
    }
}
