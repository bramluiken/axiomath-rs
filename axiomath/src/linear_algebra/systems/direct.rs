//! Direct methods for solving linear systems.
//!
//! This module implements direct methods that compute exact solutions (up to
//! round-off error) in a finite number of steps.
//!
//! # Methods
//!
//! - **LU Decomposition**: Factor A = PLU, then solve via forward/backward substitution
//! - **Gaussian Elimination**: Classical elimination with partial pivoting
//!
//! # References
//!
//! - Golub & Van Loan, "Matrix Computations", Chapter 3
//! - Trefethen & Bau, "Numerical Linear Algebra", Lectures 20-21

use crate::core::traits::Float;
use crate::linear_algebra::decomposition::lu_decompose;
use crate::linear_algebra::systems::triangular::{solve_lower_triangular, solve_upper_triangular};
use crate::linear_algebra::LinearAlgebraError;
use crate::matrix::Matrix;
use crate::vector::Vector;

/// Solves a linear system Ax = b using LU decomposition.
///
/// # Mathematical Definition
///
/// Given A (n×n) and b (n×1), find x such that:
/// ```text
/// Ax = b
/// ```
///
/// # Algorithm
///
/// 1. Compute LU decomposition with pivoting: PA = LU
/// 2. Rewrite system as: PAx = Pb, or LUx = Pb
/// 3. Let y = Ux, solve Ly = Pb using forward substitution
/// 4. Solve Ux = y using backward substitution
///
/// # Complexity
///
/// - Time: O(n³) for LU decomposition + O(n²) for solving = O(n³)
/// - Space: O(n²) for L and U matrices
///
/// # Numerical Stability
///
/// This method is numerically stable due to partial pivoting in the LU
/// decomposition. The condition number of the matrix affects accuracy.
///
/// # Arguments
///
/// * `a` - Coefficient matrix (n×n)
/// * `b` - Right-hand side vector (n)
///
/// # Returns
///
/// Solution vector x such that Ax = b.
///
/// Returns `Err` if:
/// - Matrix A is singular or numerically singular
/// - Dimensions don't match
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::vector::Vector;
/// use axiomath::linear_algebra::systems::solve_linear_system;
///
/// let a = Matrix::from_rows([
///     [2.0, 1.0, 1.0],
///     [4.0, 3.0, 3.0],
///     [8.0, 7.0, 9.0],
/// ]);
///
/// let b = Vector::from([4.0, 10.0, 24.0]);
/// let x = solve_linear_system(&a, &b).unwrap();
///
/// // Verify: Ax ≈ b
/// ```
///
/// # Panics
///
/// Panics if the matrix is not square.
pub fn solve_linear_system<T: Float, const N: usize>(
    a: &Matrix<T, N, N>,
    b: &Vector<T, N>,
) -> Result<Vector<T, N>, LinearAlgebraError> {
    // Step 1: Compute LU decomposition with pivoting
    let (l, u, perm) = lu_decompose(a)?;

    // Step 2: Apply permutation to b to get Pb
    let mut pb = Vector::<T, N>::zero();
    for i in 0..N {
        pb[i] = b[perm[i]];
    }

    // Step 3: Solve Ly = Pb using forward substitution
    let y = solve_lower_triangular(&l, &pb)?;

    // Step 4: Solve Ux = y using backward substitution
    let x = solve_upper_triangular(&u, &y)?;

    Ok(x)
}

/// Solves a linear system Ax = b using Gaussian elimination with partial pivoting.
///
/// # Mathematical Definition
///
/// Given A (n×n) and b (n×1), find x such that:
/// ```text
/// Ax = b
/// ```
///
/// # Algorithm
///
/// Classic Gaussian elimination:
///
/// ```text
/// Form augmented matrix [A|b]
/// for k = 1 to n-1:
///     Find pivot: i = argmax_{j≥k} |A_jk|
///     Swap rows k and i
///     for i = k+1 to n:
///         m = A_ik / A_kk
///         Row i -= m * Row k
/// Solve via backward substitution
/// ```
///
/// # Complexity
///
/// - Time: O(n³)
/// - Space: O(n²) for augmented matrix
///
/// # Arguments
///
/// * `a` - Coefficient matrix (n×n)
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
/// use axiomath::linear_algebra::systems::gaussian_elimination;
///
/// let a = Matrix::from_rows([
///     [2.0, 1.0],
///     [1.0, 3.0],
/// ]);
///
/// let b = Vector::from([5.0, 7.0]);
/// let x = gaussian_elimination(&a, &b).unwrap();
/// ```
pub fn gaussian_elimination<T: Float, const N: usize>(
    a: &Matrix<T, N, N>,
    b: &Vector<T, N>,
) -> Result<Vector<T, N>, LinearAlgebraError> {
    // Create augmented matrix [A|b]
    let mut aug = *a;
    let mut b_mut = *b;

    let epsilon = T::epsilon() * T::from(N).unwrap() * T::from(100.0).unwrap();

    // Forward elimination with partial pivoting
    for k in 0..N - 1 {
        // Find pivot
        let mut pivot_row = k;
        let mut max_val = aug[(k, k)].abs();

        for i in (k + 1)..N {
            let abs_val = aug[(i, k)].abs();
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
            for j in 0..N {
                let temp = aug[(k, j)];
                aug[(k, j)] = aug[(pivot_row, j)];
                aug[(pivot_row, j)] = temp;
            }
            let temp = b_mut[k];
            b_mut[k] = b_mut[pivot_row];
            b_mut[pivot_row] = temp;
        }

        // Eliminate column k below diagonal
        for i in (k + 1)..N {
            let multiplier = aug[(i, k)] / aug[(k, k)];
            for j in k..N {
                aug[(i, j)] = aug[(i, j)] - multiplier * aug[(k, j)];
            }
            b_mut[i] = b_mut[i] - multiplier * b_mut[k];
        }
    }

    // Check last diagonal element
    if aug[(N - 1, N - 1)].abs() < epsilon {
        return Err(LinearAlgebraError::SingularMatrix);
    }

    // Backward substitution
    solve_upper_triangular(&aug, &b_mut)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::matrix_vector_multiply;
    use crate::utils::comparison::approximately_equal;

    #[test]
    fn test_solve_linear_system_3x3() {
        let a = Matrix::from_rows([[2.0, 1.0, 1.0], [4.0, 3.0, 3.0], [8.0, 7.0, 9.0]]);

        let b = Vector::from([4.0, 10.0, 24.0]);
        let x = solve_linear_system(&a, &b).unwrap();

        // Verify Ax = b
        let ax = matrix_vector_multiply(&a, &x);
        for i in 0..3 {
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
    fn test_solve_linear_system_2x2() {
        let a = Matrix::from_rows([[3.0, 2.0], [1.0, 4.0]]);

        let b = Vector::from([13.0, 9.0]);
        let x = solve_linear_system(&a, &b).unwrap();

        // Verify Ax = b instead of checking explicit values
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
    fn test_solve_linear_system_identity() {
        let a = Matrix::<f64, 3, 3>::identity();
        let b = Vector::from([1.0, 2.0, 3.0]);
        let x = solve_linear_system(&a, &b).unwrap();

        // For identity matrix, x = b
        for i in 0..3 {
            assert!((x[i] - b[i]).abs() < 1e-10);
        }
    }

    #[test]
    fn test_solve_linear_system_singular() {
        let a = Matrix::from_rows([[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [1.0, 1.0, 1.0]]);

        let b = Vector::from([1.0, 2.0, 3.0]);
        let result = solve_linear_system(&a, &b);

        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), LinearAlgebraError::SingularMatrix);
    }

    #[test]
    fn test_gaussian_elimination_2x2() {
        let a = Matrix::from_rows([[2.0, 1.0], [1.0, 3.0]]);

        let b = Vector::from([5.0, 7.0]);
        let x = gaussian_elimination(&a, &b).unwrap();

        // Verify solution
        let ax = matrix_vector_multiply(&a, &x);
        for i in 0..2 {
            assert!(approximately_equal(ax[i], b[i], 1e-10, 1e-10));
        }
    }

    #[test]
    fn test_gaussian_elimination_3x3() {
        let a = Matrix::from_rows([[1.0, 2.0, 3.0], [2.0, 5.0, 3.0], [1.0, 0.0, 8.0]]);

        let b = Vector::from([3.0, 4.0, 1.0]);
        let x = gaussian_elimination(&a, &b).unwrap();

        // Verify solution
        let ax = matrix_vector_multiply(&a, &x);
        for i in 0..3 {
            assert!(approximately_equal(ax[i], b[i], 1e-10, 1e-10));
        }
    }

    #[test]
    fn test_gaussian_elimination_with_pivoting() {
        // Matrix that requires pivoting
        let a = Matrix::from_rows([[0.0, 1.0], [1.0, 1.0]]);

        let b = Vector::from([1.0, 2.0]);
        let x = gaussian_elimination(&a, &b).unwrap();

        // Verify solution
        let ax = matrix_vector_multiply(&a, &x);
        for i in 0..2 {
            assert!(approximately_equal(ax[i], b[i], 1e-10, 1e-10));
        }
    }

    #[test]
    fn test_consistency_lu_vs_gaussian() {
        let a = Matrix::from_rows([[2.0, 1.0, 3.0], [4.0, 1.0, 2.0], [1.0, 3.0, 4.0]]);

        let b = Vector::from([10.0, 12.0, 14.0]);

        let x1 = solve_linear_system(&a, &b).unwrap();
        let x2 = gaussian_elimination(&a, &b).unwrap();

        // Both methods should give the same solution
        for i in 0..3 {
            assert!(approximately_equal(x1[i], x2[i], 1e-10, 1e-10));
        }
    }
}
