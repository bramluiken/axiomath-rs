//! Matrix constructors for creating special matrices.
//!
//! This module provides functions for creating common matrix types:
//! - Zero matrix (all elements are 0)
//! - Identity matrix (I_ij = δ_ij, where δ is the Kronecker delta)
//! - Diagonal matrix (non-zero only on main diagonal)
//! - Matrices from rows or columns
//! - Scalar matrices (scalar multiple of identity)
//!
//! # Examples
//!
//! ```
//! use axiomath::matrix::Matrix;
//!
//! // Zero matrix
//! let zero = Matrix::<f64, 3, 3>::zero();
//!
//! // Identity matrix
//! let identity = Matrix::<f64, 3, 3>::identity();
//!
//! // Diagonal matrix
//! let diag = Matrix::<f64, 3, 3>::from_diagonal([1.0, 2.0, 3.0]);
//! ```

use crate::core::Scalar;
use crate::matrix::Matrix;
use num_traits::{One, Zero};

impl<T: Scalar, const ROWS: usize, const COLS: usize> Matrix<T, ROWS, COLS> {
    /// Creates a zero matrix (all elements are zero).
    ///
    /// # Mathematical Definition
    ///
    /// ```text
    /// O = [0  0  ...  0]
    ///     [0  0  ...  0]
    ///     [... ... ... ...]
    ///     [0  0  ...  0]
    /// ```
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::matrix::Matrix;
    ///
    /// let zero = Matrix::<f64, 2, 3>::zero();
    /// assert!(zero.is_zero());
    /// ```
    pub fn zero() -> Self {
        Self {
            data: [[T::zero(); COLS]; ROWS],
        }
    }

    /// Creates a matrix from row vectors.
    ///
    /// # Arguments
    ///
    /// * `rows` - Array of row arrays
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::matrix::Matrix;
    ///
    /// let m = Matrix::<f64, 2, 3>::from_rows([
    ///     [1.0, 2.0, 3.0],
    ///     [4.0, 5.0, 6.0],
    /// ]);
    /// assert_eq!(m[(0, 1)], 2.0);
    /// assert_eq!(m[(1, 2)], 6.0);
    /// ```
    pub fn from_rows(rows: [[T; COLS]; ROWS]) -> Self {
        Self { data: rows }
    }

    /// Creates a matrix from column vectors.
    ///
    /// # Arguments
    ///
    /// * `cols` - Array of column arrays
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::matrix::Matrix;
    ///
    /// let m = Matrix::<f64, 2, 3>::from_columns([
    ///     [1.0, 4.0],
    ///     [2.0, 5.0],
    ///     [3.0, 6.0],
    /// ]);
    /// assert_eq!(m[(0, 0)], 1.0);
    /// assert_eq!(m[(0, 1)], 2.0);
    /// assert_eq!(m[(1, 2)], 6.0);
    /// ```
    pub fn from_columns(cols: [[T; ROWS]; COLS]) -> Self {
        let mut data = [[T::zero(); COLS]; ROWS];
        for col_idx in 0..COLS {
            for row_idx in 0..ROWS {
                data[row_idx][col_idx] = cols[col_idx][row_idx];
            }
        }
        Self { data }
    }

    /// Creates a matrix with all elements set to the same value.
    ///
    /// # Arguments
    ///
    /// * `value` - The value to fill the matrix with
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::matrix::Matrix;
    ///
    /// let m = Matrix::<f64, 2, 3>::from_scalar(5.0);
    /// assert_eq!(m[(0, 0)], 5.0);
    /// assert_eq!(m[(1, 2)], 5.0);
    /// ```
    pub fn from_scalar(value: T) -> Self {
        Self {
            data: [[value; COLS]; ROWS],
        }
    }

    /// Creates a matrix using a function that maps indices to values.
    ///
    /// # Arguments
    ///
    /// * `f` - Function taking (row, col) and returning the element value
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::matrix::Matrix;
    ///
    /// // Create a matrix where each element equals row + col
    /// let m = Matrix::<f64, 3, 3>::from_fn(|row, col| (row + col) as f64);
    /// assert_eq!(m[(0, 0)], 0.0);
    /// assert_eq!(m[(1, 2)], 3.0);
    /// assert_eq!(m[(2, 2)], 4.0);
    /// ```
    pub fn from_fn<F>(mut f: F) -> Self
    where
        F: FnMut(usize, usize) -> T,
    {
        let mut data = [[T::zero(); COLS]; ROWS];
        for row in 0..ROWS {
            for col in 0..COLS {
                data[row][col] = f(row, col);
            }
        }
        Self { data }
    }
}

// Methods specific to square matrices
impl<T: Scalar, const N: usize> Matrix<T, N, N> {
    /// Creates an identity matrix.
    ///
    /// # Mathematical Definition
    ///
    /// The identity matrix I satisfies:
    /// ```text
    /// I_ij = δ_ij = { 1  if i = j
    ///               { 0  if i ≠ j
    /// ```
    ///
    /// where δ_ij is the Kronecker delta.
    ///
    /// # Properties
    ///
    /// - For any matrix A: AI = IA = A
    /// - det(I) = 1
    /// - I⁻¹ = I
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::matrix::Matrix;
    ///
    /// let identity = Matrix::<f64, 3, 3>::identity();
    /// assert_eq!(identity[(0, 0)], 1.0);
    /// assert_eq!(identity[(1, 1)], 1.0);
    /// assert_eq!(identity[(0, 1)], 0.0);
    /// ```
    pub fn identity() -> Self {
        let mut data = [[T::zero(); N]; N];
        for i in 0..N {
            data[i][i] = T::one();
        }
        Self { data }
    }

    /// Creates a diagonal matrix from a vector of diagonal elements.
    ///
    /// # Mathematical Definition
    ///
    /// For diagonal vector d = [d₀, d₁, ..., d_{n-1}]:
    /// ```text
    /// D_ij = { d_i  if i = j
    ///        { 0    if i ≠ j
    /// ```
    ///
    /// # Arguments
    ///
    /// * `diagonal` - Array of diagonal elements
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::matrix::Matrix;
    ///
    /// let diag = Matrix::<f64, 3, 3>::from_diagonal([1.0, 2.0, 3.0]);
    /// assert_eq!(diag[(0, 0)], 1.0);
    /// assert_eq!(diag[(1, 1)], 2.0);
    /// assert_eq!(diag[(2, 2)], 3.0);
    /// assert_eq!(diag[(0, 1)], 0.0);
    /// ```
    pub fn from_diagonal(diagonal: [T; N]) -> Self {
        let mut data = [[T::zero(); N]; N];
        for i in 0..N {
            data[i][i] = diagonal[i];
        }
        Self { data }
    }

    /// Creates a scalar matrix (scalar multiple of the identity).
    ///
    /// # Mathematical Definition
    ///
    /// ```text
    /// λI = [λ  0  ...  0]
    ///      [0  λ  ...  0]
    ///      [... ... ... ...]
    ///      [0  0  ...  λ]
    /// ```
    ///
    /// # Arguments
    ///
    /// * `scalar` - The scalar value to place on the diagonal
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::matrix::Matrix;
    ///
    /// let scaled = Matrix::<f64, 3, 3>::scalar(5.0);
    /// assert_eq!(scaled[(0, 0)], 5.0);
    /// assert_eq!(scaled[(2, 2)], 5.0);
    /// assert_eq!(scaled[(0, 1)], 0.0);
    /// ```
    pub fn scalar(scalar: T) -> Self {
        let mut data = [[T::zero(); N]; N];
        for i in 0..N {
            data[i][i] = scalar;
        }
        Self { data }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero_matrix() {
        let zero = Matrix::<f64, 3, 4>::zero();
        for row in 0..3 {
            for col in 0..4 {
                assert_eq!(zero[(row, col)], 0.0);
            }
        }
        assert!(zero.is_zero());
    }

    #[test]
    fn test_from_rows() {
        let m = Matrix::<f64, 2, 3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        assert_eq!(m[(0, 0)], 1.0);
        assert_eq!(m[(0, 1)], 2.0);
        assert_eq!(m[(0, 2)], 3.0);
        assert_eq!(m[(1, 0)], 4.0);
        assert_eq!(m[(1, 1)], 5.0);
        assert_eq!(m[(1, 2)], 6.0);
    }

    #[test]
    fn test_from_columns() {
        let m = Matrix::<f64, 2, 3>::from_columns([[1.0, 4.0], [2.0, 5.0], [3.0, 6.0]]);
        assert_eq!(m[(0, 0)], 1.0);
        assert_eq!(m[(0, 1)], 2.0);
        assert_eq!(m[(0, 2)], 3.0);
        assert_eq!(m[(1, 0)], 4.0);
        assert_eq!(m[(1, 1)], 5.0);
        assert_eq!(m[(1, 2)], 6.0);
    }

    #[test]
    fn test_from_scalar() {
        let m = Matrix::<f64, 2, 3>::from_scalar(7.0);
        for row in 0..2 {
            for col in 0..3 {
                assert_eq!(m[(row, col)], 7.0);
            }
        }
    }

    #[test]
    fn test_from_fn() {
        let m = Matrix::<f64, 3, 3>::from_fn(|row, col| (row * 3 + col) as f64);
        assert_eq!(m[(0, 0)], 0.0);
        assert_eq!(m[(0, 1)], 1.0);
        assert_eq!(m[(0, 2)], 2.0);
        assert_eq!(m[(1, 0)], 3.0);
        assert_eq!(m[(2, 2)], 8.0);
    }

    #[test]
    fn test_identity_2x2() {
        let identity = Matrix::<f64, 2, 2>::identity();
        assert_eq!(identity[(0, 0)], 1.0);
        assert_eq!(identity[(1, 1)], 1.0);
        assert_eq!(identity[(0, 1)], 0.0);
        assert_eq!(identity[(1, 0)], 0.0);
    }

    #[test]
    fn test_identity_3x3() {
        let identity = Matrix::<f64, 3, 3>::identity();
        for i in 0..3 {
            for j in 0..3 {
                if i == j {
                    assert_eq!(identity[(i, j)], 1.0);
                } else {
                    assert_eq!(identity[(i, j)], 0.0);
                }
            }
        }
    }

    #[test]
    fn test_from_diagonal() {
        let diag = Matrix::<f64, 3, 3>::from_diagonal([1.0, 2.0, 3.0]);
        assert_eq!(diag[(0, 0)], 1.0);
        assert_eq!(diag[(1, 1)], 2.0);
        assert_eq!(diag[(2, 2)], 3.0);
        assert_eq!(diag[(0, 1)], 0.0);
        assert_eq!(diag[(0, 2)], 0.0);
        assert_eq!(diag[(1, 0)], 0.0);
        assert_eq!(diag[(1, 2)], 0.0);
        assert_eq!(diag[(2, 0)], 0.0);
        assert_eq!(diag[(2, 1)], 0.0);
    }

    #[test]
    fn test_scalar_matrix() {
        let scaled = Matrix::<f64, 4, 4>::scalar(5.0);
        for i in 0..4 {
            for j in 0..4 {
                if i == j {
                    assert_eq!(scaled[(i, j)], 5.0);
                } else {
                    assert_eq!(scaled[(i, j)], 0.0);
                }
            }
        }
    }

    #[test]
    fn test_from_rows_and_columns_equivalence() {
        let from_rows = Matrix::<f64, 2, 3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        let from_cols = Matrix::<f64, 2, 3>::from_columns([[1.0, 4.0], [2.0, 5.0], [3.0, 6.0]]);
        assert_eq!(from_rows, from_cols);
    }
}
