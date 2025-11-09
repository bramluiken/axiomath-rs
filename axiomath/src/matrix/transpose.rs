//! Matrix transpose operations.
//!
//! This module provides functions for transposing matrices, which swaps
//! rows and columns.
//!
//! # Mathematical Background
//!
//! The transpose of a matrix A, denoted A^T or A', is obtained by reflecting
//! the matrix over its main diagonal:
//!
//! ```text
//! If A = [a  b  c]     then A^T = [a  d]
//!        [d  e  f]               [b  e]
//!                                 [c  f]
//! ```
//!
//! Formally: (A^T)_ij = A_ji
//!
//! # Properties
//!
//! - (A^T)^T = A (involutive)
//! - (A + B)^T = A^T + B^T (distributive over addition)
//! - (λA)^T = λA^T (compatible with scalar multiplication)
//! - (AB)^T = B^T A^T (reverses order in products)
//! - det(A^T) = det(A) (determinant unchanged)
//!
//! # Examples
//!
//! ```
//! use axiomath::matrix::{Matrix, transpose};
//!
//! let a = Matrix::<f64, 2, 3>::from_rows([
//!     [1.0, 2.0, 3.0],
//!     [4.0, 5.0, 6.0],
//! ]);
//!
//! let at = transpose(&a);
//! assert_eq!(at.dimensions(), (3, 2));
//! assert_eq!(at[(0, 0)], 1.0);
//! assert_eq!(at[(2, 1)], 6.0);
//! ```

use crate::core::Scalar;
use crate::matrix::Matrix;

/// Transposes a matrix by swapping rows and columns.
///
/// # Mathematical Definition
///
/// For an M×N matrix A, the transpose A^T is an N×M matrix where:
/// ```text
/// (A^T)_ij = A_ji
/// ```
///
/// # Type Transformation
///
/// ```text
/// Matrix<T, M, N> → Matrix<T, N, M>
/// ```
///
/// Note that the dimensions are swapped.
///
/// # Complexity
///
/// - Time: O(M×N)
/// - Space: O(M×N) (creates a new matrix)
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, transpose};
///
/// let a = Matrix::<f64, 2, 3>::from_rows([
///     [1.0, 2.0, 3.0],
///     [4.0, 5.0, 6.0],
/// ]);
///
/// let at = transpose(&a);
///
/// // Dimensions are swapped
/// assert_eq!(at.rows(), 3);
/// assert_eq!(at.cols(), 2);
///
/// // Elements are reflected across diagonal
/// assert_eq!(at[(0, 0)], 1.0);
/// assert_eq!(at[(0, 1)], 4.0);
/// assert_eq!(at[(1, 0)], 2.0);
/// assert_eq!(at[(2, 1)], 6.0);
/// ```
pub fn transpose<T: Scalar, const ROWS: usize, const COLS: usize>(
    matrix: &Matrix<T, ROWS, COLS>,
) -> Matrix<T, COLS, ROWS> {
    Matrix::<T, COLS, ROWS>::from_fn(|row, col| matrix[(col, row)])
}

/// Transposes a square matrix in-place.
///
/// # Mathematical Definition
///
/// Swaps elements A_ij ↔ A_ji for all i ≠ j.
///
/// # Complexity
///
/// - Time: O(N²)
/// - Space: O(1) (in-place modification)
///
/// # Advantages
///
/// - More memory efficient than creating a new matrix
/// - Faster for large square matrices
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
///
/// let mut a = Matrix::<f64, 3, 3>::from_rows([
///     [1.0, 2.0, 3.0],
///     [4.0, 5.0, 6.0],
///     [7.0, 8.0, 9.0],
/// ]);
///
/// a.transpose_in_place();
///
/// assert_eq!(a[(0, 1)], 4.0);
/// assert_eq!(a[(1, 0)], 2.0);
/// assert_eq!(a[(0, 2)], 7.0);
/// assert_eq!(a[(2, 0)], 3.0);
/// ```
pub fn transpose_in_place<T: Scalar, const N: usize>(matrix: &mut Matrix<T, N, N>) {
    for i in 0..N {
        for j in (i + 1)..N {
            // Swap A_ij and A_ji
            let temp = matrix[(i, j)];
            matrix[(i, j)] = matrix[(j, i)];
            matrix[(j, i)] = temp;
        }
    }
}

impl<T: Scalar, const N: usize> Matrix<T, N, N> {
    /// Transposes the matrix in-place (square matrices only).
    ///
    /// This is a convenience method that calls `transpose_in_place`.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::matrix::Matrix;
    ///
    /// let mut m = Matrix::<f64, 2, 2>::from_rows([
    ///     [1.0, 2.0],
    ///     [3.0, 4.0],
    /// ]);
    ///
    /// m.transpose_in_place();
    ///
    /// assert_eq!(m[(0, 1)], 3.0);
    /// assert_eq!(m[(1, 0)], 2.0);
    /// ```
    pub fn transpose_in_place(&mut self) {
        transpose_in_place(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_transpose_2x3() {
        let a = Matrix::<f64, 2, 3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);

        let at = transpose(&a);

        assert_eq!(at.rows(), 3);
        assert_eq!(at.cols(), 2);
        assert_eq!(at[(0, 0)], 1.0);
        assert_eq!(at[(0, 1)], 4.0);
        assert_eq!(at[(1, 0)], 2.0);
        assert_eq!(at[(1, 1)], 5.0);
        assert_eq!(at[(2, 0)], 3.0);
        assert_eq!(at[(2, 1)], 6.0);
    }

    #[test]
    fn test_transpose_square() {
        let a = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ]);

        let at = transpose(&a);

        assert_eq!(at.rows(), 3);
        assert_eq!(at.cols(), 3);
        assert_eq!(at[(0, 0)], 1.0);
        assert_eq!(at[(0, 1)], 4.0);
        assert_eq!(at[(0, 2)], 7.0);
        assert_eq!(at[(1, 0)], 2.0);
        assert_eq!(at[(2, 0)], 3.0);
        assert_eq!(at[(2, 2)], 9.0);
    }

    #[test]
    fn test_transpose_involutive() {
        // (A^T)^T = A
        let a = Matrix::<f64, 2, 3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);

        let at = transpose(&a);
        let att = transpose(&at);

        assert_eq!(att, a);
    }

    #[test]
    fn test_transpose_distributive() {
        // (A + B)^T = A^T + B^T
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let b = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);

        let left = transpose(&(a + b));
        let right = transpose(&a) + transpose(&b);

        assert_eq!(left, right);
    }

    #[test]
    fn test_transpose_scalar_multiplication() {
        // (λA)^T = λA^T
        let a = Matrix::<f64, 2, 3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        let scalar = 2.5;

        let left = transpose(&(a * scalar));
        let right = transpose(&a) * scalar;

        assert_eq!(left, right);
    }

    #[test]
    fn test_transpose_identity() {
        let identity = Matrix::<f64, 4, 4>::identity();
        let transposed = transpose(&identity);

        assert_eq!(transposed, identity);
    }

    #[test]
    fn test_transpose_diagonal() {
        let diag = Matrix::<f64, 3, 3>::from_diagonal([1.0, 2.0, 3.0]);
        let transposed = transpose(&diag);

        assert_eq!(transposed, diag);
    }

    #[test]
    fn test_transpose_in_place_2x2() {
        let mut m = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);

        transpose_in_place(&mut m);

        assert_eq!(m[(0, 0)], 1.0);
        assert_eq!(m[(0, 1)], 3.0);
        assert_eq!(m[(1, 0)], 2.0);
        assert_eq!(m[(1, 1)], 4.0);
    }

    #[test]
    fn test_transpose_in_place_3x3() {
        let mut m = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ]);

        transpose_in_place(&mut m);

        assert_eq!(m[(0, 0)], 1.0);
        assert_eq!(m[(0, 1)], 4.0);
        assert_eq!(m[(0, 2)], 7.0);
        assert_eq!(m[(1, 0)], 2.0);
        assert_eq!(m[(1, 1)], 5.0);
        assert_eq!(m[(1, 2)], 8.0);
        assert_eq!(m[(2, 0)], 3.0);
        assert_eq!(m[(2, 1)], 6.0);
        assert_eq!(m[(2, 2)], 9.0);
    }

    #[test]
    fn test_transpose_in_place_method() {
        let mut m = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);

        m.transpose_in_place();

        assert_eq!(m[(0, 1)], 3.0);
        assert_eq!(m[(1, 0)], 2.0);
    }

    #[test]
    fn test_transpose_in_place_involutive() {
        let original = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ]);

        let mut m = original;
        transpose_in_place(&mut m);
        transpose_in_place(&mut m);

        assert_eq!(m, original);
    }

    #[test]
    fn test_transpose_in_place_symmetric() {
        // Symmetric matrix should be unchanged by transpose
        let mut m = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [2.0, 4.0, 5.0],
            [3.0, 5.0, 6.0],
        ]);

        let original = m;
        transpose_in_place(&mut m);

        assert_eq!(m, original);
    }

    #[test]
    fn test_transpose_vs_transpose_in_place() {
        let original = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ]);

        let via_transpose = transpose(&original);

        let mut via_in_place = original;
        transpose_in_place(&mut via_in_place);

        assert_eq!(via_transpose, via_in_place);
    }

    #[test]
    fn test_transpose_zero_matrix() {
        let zero = Matrix::<f64, 3, 4>::zero();
        let transposed = transpose(&zero);

        assert_eq!(transposed.rows(), 4);
        assert_eq!(transposed.cols(), 3);
        assert!(transposed.is_zero());
    }
}
