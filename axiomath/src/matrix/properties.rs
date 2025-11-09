//! Matrix properties and characteristic functions.
//!
//! This module provides functions for computing matrix properties such as:
//! - Trace (sum of diagonal elements)
//! - Determinant
//! - Boolean property checks (symmetric, diagonal, triangular, etc.)
//!
//! # Mathematical Background
//!
//! ## Trace
//!
//! The trace of a square matrix is the sum of its diagonal elements:
//! ```text
//! tr(A) = Σ(i=1 to n) A_ii
//! ```
//!
//! ## Determinant
//!
//! The determinant is a scalar value that encodes important properties:
//! - det(A) = 0 if and only if A is singular (not invertible)
//! - det(AB) = det(A)·det(B)
//! - det(A^T) = det(A)
//!
//! # Examples
//!
//! ```
//! use axiomath::matrix::{Matrix, trace, determinant};
//!
//! let m = Matrix::<f64, 2, 2>::from_rows([
//!     [1.0, 2.0],
//!     [3.0, 4.0],
//! ]);
//!
//! assert_eq!(trace(&m), 5.0);  // 1 + 4
//! assert_eq!(determinant(&m), -2.0);  // 1*4 - 2*3
//! ```

use crate::core::{Real, Scalar};
use crate::matrix::Matrix;
use num_traits::Zero;

/// Computes the trace of a square matrix.
///
/// # Mathematical Definition
///
/// ```text
/// tr(A) = Σ(i=1 to N) A_ii
/// ```
///
/// # Properties
///
/// - tr(A + B) = tr(A) + tr(B)
/// - tr(λA) = λ·tr(A)
/// - tr(A^T) = tr(A)
/// - tr(AB) = tr(BA) (cyclic property)
///
/// # Complexity
///
/// - Time: O(N)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, trace};
///
/// let m = Matrix::<f64, 3, 3>::from_rows([
///     [1.0, 2.0, 3.0],
///     [4.0, 5.0, 6.0],
///     [7.0, 8.0, 9.0],
/// ]);
///
/// assert_eq!(trace(&m), 15.0);  // 1 + 5 + 9
/// ```
pub fn trace<T: Scalar, const N: usize>(matrix: &Matrix<T, N, N>) -> T {
    let mut sum = T::zero();
    for i in 0..N {
        sum += matrix[(i, i)];
    }
    sum
}

/// Computes the determinant of a 2×2 matrix.
///
/// # Mathematical Definition
///
/// ```text
/// det([a  b]) = ad - bc
///     ([c  d])
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, determinant_2x2};
///
/// let m = Matrix::<f64, 2, 2>::from_rows([
///     [1.0, 2.0],
///     [3.0, 4.0],
/// ]);
///
/// assert_eq!(determinant_2x2(&m), -2.0);  // 1*4 - 2*3
/// ```
pub fn determinant_2x2<T: Scalar>(matrix: &Matrix<T, 2, 2>) -> T {
    matrix[(0, 0)] * matrix[(1, 1)] - matrix[(0, 1)] * matrix[(1, 0)]
}

/// Computes the determinant of a 3×3 matrix using the rule of Sarrus.
///
/// # Mathematical Definition
///
/// ```text
/// det(A) = a₁₁a₂₂a₃₃ + a₁₂a₂₃a₃₁ + a₁₃a₂₁a₃₂
///        - a₁₃a₂₂a₃₁ - a₁₂a₂₁a₃₃ - a₁₁a₂₃a₃₂
/// ```
///
/// # Algorithm
///
/// Uses the rule of Sarrus for direct computation of 3×3 determinants.
///
/// # Complexity
///
/// - Time: O(1) (constant 12 multiplications and 5 additions/subtractions)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, determinant_3x3};
///
/// let m = Matrix::<f64, 3, 3>::from_rows([
///     [1.0, 2.0, 3.0],
///     [4.0, 5.0, 6.0],
///     [7.0, 8.0, 9.0],
/// ]);
///
/// assert_eq!(determinant_3x3(&m), 0.0);  // Singular matrix
/// ```
pub fn determinant_3x3<T: Scalar>(matrix: &Matrix<T, 3, 3>) -> T {
    let a = matrix[(0, 0)];
    let b = matrix[(0, 1)];
    let c = matrix[(0, 2)];
    let d = matrix[(1, 0)];
    let e = matrix[(1, 1)];
    let f = matrix[(1, 2)];
    let g = matrix[(2, 0)];
    let h = matrix[(2, 1)];
    let i = matrix[(2, 2)];

    // Rule of Sarrus
    a * e * i + b * f * g + c * d * h - c * e * g - b * d * i - a * f * h
}

/// Computes the determinant of a 4×4 matrix using Laplace expansion.
///
/// # Algorithm
///
/// Expands along the first row using cofactors of 3×3 minors.
///
/// # Complexity
///
/// - Time: O(N!) but optimized for 4×4
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, determinant_4x4};
///
/// let m = Matrix::<f64, 4, 4>::from_rows([
///     [1.0, 0.0, 0.0, 0.0],
///     [0.0, 2.0, 0.0, 0.0],
///     [0.0, 0.0, 3.0, 0.0],
///     [0.0, 0.0, 0.0, 4.0],
/// ]);
///
/// assert_eq!(determinant_4x4(&m), 24.0);  // Product of diagonal
/// ```
pub fn determinant_4x4<T: Scalar>(matrix: &Matrix<T, 4, 4>) -> T {
    let mut det = T::zero();

    for j in 0..4 {
        // Extract 3×3 minor by removing row 0 and column j
        let minor = minor_3x3_from_4x4(matrix, 0, j);
        let cofactor = determinant_3x3(&minor);
        let sign = if j % 2 == 0 { T::one() } else { T::zero() - T::one() };
        det += sign * matrix[(0, j)] * cofactor;
    }

    det
}

// Helper to extract a 3×3 minor from a 4×4 matrix
fn minor_3x3_from_4x4<T: Scalar>(
    matrix: &Matrix<T, 4, 4>,
    row_remove: usize,
    col_remove: usize,
) -> Matrix<T, 3, 3> {
    Matrix::from_fn(|i, j| {
        let source_i = if i < row_remove { i } else { i + 1 };
        let source_j = if j < col_remove { j } else { j + 1 };
        matrix[(source_i, source_j)]
    })
}

/// Computes the determinant of an N×N matrix.
///
/// # Algorithm
///
/// - For 1×1: Returns the single element
/// - For 2×2: Direct formula
/// - For 3×3: Rule of Sarrus
/// - For 4×4: Laplace expansion
/// - For larger: Returns zero (not implemented - use LU decomposition in Phase 5)
///
/// # Warning
///
/// For matrices larger than 4×4, this returns `T::zero()` as a placeholder.
/// For production use with large matrices, use LU decomposition (Phase 5).
///
/// # Complexity
///
/// - Time: O(1) for N ≤ 4
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, determinant};
///
/// let m = Matrix::<f64, 2, 2>::from_rows([
///     [1.0, 2.0],
///     [3.0, 4.0],
/// ]);
///
/// assert_eq!(determinant(&m), -2.0);
/// ```
pub fn determinant<T: Scalar, const N: usize>(matrix: &Matrix<T, N, N>) -> T {
    match N {
        1 => matrix[(0, 0)],
        2 => determinant_2x2_generic(matrix),
        3 => determinant_3x3_generic(matrix),
        4 => determinant_4x4_generic(matrix),
        _ => T::zero(), // Placeholder - implement LU decomposition in Phase 5
    }
}

// Generic wrappers that work with any size (needed to satisfy the type checker)
fn determinant_2x2_generic<T: Scalar, const N: usize>(matrix: &Matrix<T, N, N>) -> T {
    matrix[(0, 0)] * matrix[(1, 1)] - matrix[(0, 1)] * matrix[(1, 0)]
}

fn determinant_3x3_generic<T: Scalar, const N: usize>(matrix: &Matrix<T, N, N>) -> T {
    let a = matrix[(0, 0)];
    let b = matrix[(0, 1)];
    let c = matrix[(0, 2)];
    let d = matrix[(1, 0)];
    let e = matrix[(1, 1)];
    let f = matrix[(1, 2)];
    let g = matrix[(2, 0)];
    let h = matrix[(2, 1)];
    let i = matrix[(2, 2)];

    a * e * i + b * f * g + c * d * h - c * e * g - b * d * i - a * f * h
}

fn determinant_4x4_generic<T: Scalar, const N: usize>(matrix: &Matrix<T, N, N>) -> T {
    let mut det = T::zero();

    for j in 0..4 {
        let cofactor = minor_3x3_from_nxn_generic(matrix, 0, j);
        let cofactor_det = determinant_3x3(&cofactor);
        let sign = if j % 2 == 0 { T::one() } else { T::zero() - T::one() };
        det += sign * matrix[(0, j)] * cofactor_det;
    }

    det
}

fn minor_3x3_from_nxn_generic<T: Scalar, const N: usize>(
    matrix: &Matrix<T, N, N>,
    row_remove: usize,
    col_remove: usize,
) -> Matrix<T, 3, 3> {
    Matrix::from_fn(|i, j| {
        let source_i = if i < row_remove { i } else { i + 1 };
        let source_j = if j < col_remove { j } else { j + 1 };
        matrix[(source_i, source_j)]
    })
}

/// Checks if a matrix is symmetric.
///
/// # Mathematical Definition
///
/// A matrix A is symmetric if A = A^T, i.e., A_ij = A_ji for all i, j.
///
/// # Complexity
///
/// - Time: O(N²)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, is_symmetric};
///
/// let sym = Matrix::<f64, 3, 3>::from_rows([
///     [1.0, 2.0, 3.0],
///     [2.0, 4.0, 5.0],
///     [3.0, 5.0, 6.0],
/// ]);
/// assert!(is_symmetric(&sym));
///
/// let not_sym = Matrix::<f64, 2, 2>::from_rows([
///     [1.0, 2.0],
///     [3.0, 4.0],
/// ]);
/// assert!(!is_symmetric(&not_sym));
/// ```
pub fn is_symmetric<T: Scalar, const N: usize>(matrix: &Matrix<T, N, N>) -> bool {
    for i in 0..N {
        for j in (i + 1)..N {
            if matrix[(i, j)] != matrix[(j, i)] {
                return false;
            }
        }
    }
    true
}

/// Checks if a matrix is diagonal.
///
/// # Mathematical Definition
///
/// A matrix A is diagonal if A_ij = 0 for all i ≠ j.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, is_diagonal};
///
/// let diag = Matrix::<f64, 3, 3>::from_diagonal([1.0, 2.0, 3.0]);
/// assert!(is_diagonal(&diag));
///
/// let identity = Matrix::<f64, 3, 3>::identity();
/// assert!(is_diagonal(&identity));
///
/// let not_diag = Matrix::<f64, 2, 2>::from_rows([
///     [1.0, 2.0],
///     [0.0, 3.0],
/// ]);
/// assert!(!is_diagonal(&not_diag));
/// ```
pub fn is_diagonal<T: Scalar, const N: usize>(matrix: &Matrix<T, N, N>) -> bool {
    for i in 0..N {
        for j in 0..N {
            if i != j && matrix[(i, j)] != T::zero() {
                return false;
            }
        }
    }
    true
}

/// Checks if a matrix is upper triangular.
///
/// # Mathematical Definition
///
/// A matrix A is upper triangular if A_ij = 0 for all i > j.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, is_upper_triangular};
///
/// let upper = Matrix::<f64, 3, 3>::from_rows([
///     [1.0, 2.0, 3.0],
///     [0.0, 4.0, 5.0],
///     [0.0, 0.0, 6.0],
/// ]);
/// assert!(is_upper_triangular(&upper));
/// ```
pub fn is_upper_triangular<T: Scalar, const N: usize>(matrix: &Matrix<T, N, N>) -> bool {
    for i in 0..N {
        for j in 0..i {
            if matrix[(i, j)] != T::zero() {
                return false;
            }
        }
    }
    true
}

/// Checks if a matrix is lower triangular.
///
/// # Mathematical Definition
///
/// A matrix A is lower triangular if A_ij = 0 for all i < j.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, is_lower_triangular};
///
/// let lower = Matrix::<f64, 3, 3>::from_rows([
///     [1.0, 0.0, 0.0],
///     [2.0, 3.0, 0.0],
///     [4.0, 5.0, 6.0],
/// ]);
/// assert!(is_lower_triangular(&lower));
/// ```
pub fn is_lower_triangular<T: Scalar, const N: usize>(matrix: &Matrix<T, N, N>) -> bool {
    for i in 0..N {
        for j in (i + 1)..N {
            if matrix[(i, j)] != T::zero() {
                return false;
            }
        }
    }
    true
}

/// Checks if a matrix is singular (determinant is zero).
///
/// # Mathematical Definition
///
/// A matrix A is singular if det(A) = 0.
///
/// # Note
///
/// For floating-point matrices, uses an epsilon tolerance to handle
/// numerical errors.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, is_singular};
///
/// let singular = Matrix::<f64, 2, 2>::from_rows([
///     [1.0, 2.0],
///     [2.0, 4.0],  // Row 2 = 2 * Row 1
/// ]);
/// assert!(is_singular(&singular, 1e-10));
///
/// let non_singular = Matrix::<f64, 2, 2>::identity();
/// assert!(!is_singular(&non_singular, 1e-10));
/// ```
pub fn is_singular<T: Real, const N: usize>(matrix: &Matrix<T, N, N>, epsilon: T) -> bool {
    determinant(matrix).abs() < epsilon
}

/// Checks if a matrix is invertible (non-singular).
///
/// # Mathematical Definition
///
/// A matrix A is invertible if det(A) ≠ 0.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, is_invertible};
///
/// let invertible = Matrix::<f64, 2, 2>::identity();
/// assert!(is_invertible(&invertible, 1e-10));
/// ```
pub fn is_invertible<T: Real, const N: usize>(matrix: &Matrix<T, N, N>, epsilon: T) -> bool {
    !is_singular(matrix, epsilon)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trace_2x2() {
        let m = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        assert_eq!(trace(&m), 5.0);
    }

    #[test]
    fn test_trace_3x3() {
        let m = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ]);
        assert_eq!(trace(&m), 15.0);
    }

    #[test]
    fn test_trace_identity() {
        let identity = Matrix::<f64, 4, 4>::identity();
        assert_eq!(trace(&identity), 4.0);
    }

    #[test]
    fn test_trace_additive() {
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let b = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);

        assert_eq!(trace(&(a + b)), trace(&a) + trace(&b));
    }

    #[test]
    fn test_determinant_2x2() {
        let m = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        assert_eq!(determinant_2x2(&m), -2.0);
    }

    #[test]
    fn test_determinant_2x2_identity() {
        let identity = Matrix::<f64, 2, 2>::identity();
        assert_eq!(determinant_2x2(&identity), 1.0);
    }

    #[test]
    fn test_determinant_3x3() {
        let m = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ]);
        assert_eq!(determinant_3x3(&m), 0.0); // Singular
    }

    #[test]
    fn test_determinant_3x3_identity() {
        let identity = Matrix::<f64, 3, 3>::identity();
        assert_eq!(determinant_3x3(&identity), 1.0);
    }

    #[test]
    fn test_determinant_3x3_example() {
        let m = Matrix::<f64, 3, 3>::from_rows([
            [2.0, -1.0, 0.0],
            [-1.0, 2.0, -1.0],
            [0.0, -1.0, 2.0],
        ]);
        assert_eq!(determinant_3x3(&m), 4.0);
    }

    #[test]
    fn test_determinant_general() {
        let m2 = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        assert_eq!(determinant(&m2), -2.0);

        let m3 = Matrix::<f64, 3, 3>::from_rows([
            [2.0, -1.0, 0.0],
            [-1.0, 2.0, -1.0],
            [0.0, -1.0, 2.0],
        ]);
        assert_eq!(determinant(&m3), 4.0);
    }

    #[test]
    fn test_determinant_4x4() {
        let m = Matrix::<f64, 4, 4>::from_rows([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 2.0, 0.0, 0.0],
            [0.0, 0.0, 3.0, 0.0],
            [0.0, 0.0, 0.0, 4.0],
        ]);
        assert_eq!(determinant(&m), 24.0); // Product of diagonal elements
    }

    #[test]
    fn test_is_symmetric() {
        let sym = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [2.0, 4.0, 5.0],
            [3.0, 5.0, 6.0],
        ]);
        assert!(is_symmetric(&sym));

        let not_sym = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        assert!(!is_symmetric(&not_sym));
    }

    #[test]
    fn test_is_diagonal() {
        let diag = Matrix::<f64, 3, 3>::from_diagonal([1.0, 2.0, 3.0]);
        assert!(is_diagonal(&diag));

        let identity = Matrix::<f64, 3, 3>::identity();
        assert!(is_diagonal(&identity));

        let zero = Matrix::<f64, 3, 3>::zero();
        assert!(is_diagonal(&zero));

        let not_diag = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [0.0, 3.0]]);
        assert!(!is_diagonal(&not_diag));
    }

    #[test]
    fn test_is_upper_triangular() {
        let upper = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [0.0, 4.0, 5.0],
            [0.0, 0.0, 6.0],
        ]);
        assert!(is_upper_triangular(&upper));

        let identity = Matrix::<f64, 3, 3>::identity();
        assert!(is_upper_triangular(&identity));

        let not_upper = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        assert!(!is_upper_triangular(&not_upper));
    }

    #[test]
    fn test_is_lower_triangular() {
        let lower = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 0.0, 0.0],
            [2.0, 3.0, 0.0],
            [4.0, 5.0, 6.0],
        ]);
        assert!(is_lower_triangular(&lower));

        let identity = Matrix::<f64, 3, 3>::identity();
        assert!(is_lower_triangular(&identity));

        let not_lower = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [0.0, 3.0]]);
        assert!(!is_lower_triangular(&not_lower));
    }

    #[test]
    fn test_is_singular() {
        let singular = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [2.0, 4.0]]);
        assert!(is_singular(&singular, 1e-10));

        let non_singular = Matrix::<f64, 2, 2>::identity();
        assert!(!is_singular(&non_singular, 1e-10));
    }

    #[test]
    fn test_is_invertible() {
        let invertible = Matrix::<f64, 2, 2>::identity();
        assert!(is_invertible(&invertible, 1e-10));

        let not_invertible = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [2.0, 4.0]]);
        assert!(!is_invertible(&not_invertible, 1e-10));
    }

    #[test]
    fn test_diagonal_is_symmetric() {
        let diag = Matrix::<f64, 3, 3>::from_diagonal([1.0, 2.0, 3.0]);
        assert!(is_diagonal(&diag));
        assert!(is_symmetric(&diag));
    }

    #[test]
    fn test_diagonal_is_both_triangular() {
        let diag = Matrix::<f64, 3, 3>::from_diagonal([1.0, 2.0, 3.0]);
        assert!(is_upper_triangular(&diag));
        assert!(is_lower_triangular(&diag));
    }
}
