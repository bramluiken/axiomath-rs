//! Core matrix type definition using const generics.
//!
//! This module defines the fundamental `Matrix<T, ROWS, COLS>` type, which provides
//! compile-time dimension checking for matrix operations.
//!
//! # Mathematical Background
//!
//! A matrix is a rectangular array of numbers arranged in rows and columns.
//! For a matrix A with M rows and N columns (an M×N matrix):
//!
//! ```text
//! A = [a₁₁  a₁₂  ...  a₁ₙ]
//!     [a₂₁  a₂₂  ...  a₂ₙ]
//!     [... ...  ...  ...]
//!     [aₘ₁  aₘ₂  ...  aₘₙ]
//! ```
//!
//! Each element a_ij is identified by its row index i and column index j.
//!
//! # Design Decisions
//!
//! - **Row-major storage**: Data is stored as `[[T; COLS]; ROWS]` for cache efficiency
//! - **Const generics**: Dimensions are compile-time constants for type safety
//! - **Zero-copy indexing**: Direct access to underlying array where possible
//! - **Generic over T**: Works with any type implementing `Scalar` trait
//!
//! # Examples
//!
//! ```
//! use axiomath::matrix::Matrix;
//! use axiomath::core::Scalar;
//!
//! // Create a 2×3 matrix
//! let m = Matrix::<f64, 2, 3>::from_rows([
//!     [1.0, 2.0, 3.0],
//!     [4.0, 5.0, 6.0],
//! ]);
//!
//! // Access elements
//! assert_eq!(m[(0, 0)], 1.0);
//! assert_eq!(m[(1, 2)], 6.0);
//!
//! // Get dimensions
//! assert_eq!(m.rows(), 2);
//! assert_eq!(m.cols(), 3);
//! ```

use crate::core::Scalar;
use num_traits::Zero;
use std::fmt;
use std::ops::{Index, IndexMut};

/// A generic matrix with compile-time dimensions.
///
/// # Type Parameters
///
/// - `T`: The element type (must implement `Scalar`)
/// - `ROWS`: Number of rows (compile-time constant)
/// - `COLS`: Number of columns (compile-time constant)
///
/// # Storage
///
/// Elements are stored in row-major order as a 2D array `[[T; COLS]; ROWS]`.
/// This provides:
/// - Contiguous memory layout
/// - Cache-friendly access patterns
/// - Zero-cost abstractions
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
///
/// // 2×2 identity matrix
/// let identity = Matrix::<f64, 2, 2>::identity();
/// assert_eq!(identity[(0, 0)], 1.0);
/// assert_eq!(identity[(0, 1)], 0.0);
/// ```
#[derive(Copy, Clone, PartialEq, Eq, Hash)]
pub struct Matrix<T, const ROWS: usize, const COLS: usize> {
    /// Internal data storage in row-major format
    pub(crate) data: [[T; COLS]; ROWS],
}

impl<T: Scalar, const ROWS: usize, const COLS: usize> Matrix<T, ROWS, COLS> {
    /// Returns the number of rows in the matrix.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::matrix::Matrix;
    ///
    /// let m = Matrix::<f64, 3, 4>::zero();
    /// assert_eq!(m.rows(), 3);
    /// ```
    #[inline]
    pub const fn rows(&self) -> usize {
        ROWS
    }

    /// Returns the number of columns in the matrix.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::matrix::Matrix;
    ///
    /// let m = Matrix::<f64, 3, 4>::zero();
    /// assert_eq!(m.cols(), 4);
    /// ```
    #[inline]
    pub const fn cols(&self) -> usize {
        COLS
    }

    /// Returns the dimensions as a tuple (rows, cols).
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::matrix::Matrix;
    ///
    /// let m = Matrix::<f64, 2, 3>::zero();
    /// assert_eq!(m.dimensions(), (2, 3));
    /// ```
    #[inline]
    pub const fn dimensions(&self) -> (usize, usize) {
        (ROWS, COLS)
    }

    /// Returns a reference to the underlying data array.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::matrix::Matrix;
    ///
    /// let m = Matrix::<f64, 2, 2>::identity();
    /// let data = m.as_array();
    /// assert_eq!(data[0][0], 1.0);
    /// ```
    #[inline]
    pub const fn as_array(&self) -> &[[T; COLS]; ROWS] {
        &self.data
    }

    /// Returns a mutable reference to the underlying data array.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::matrix::Matrix;
    ///
    /// let mut m = Matrix::<f64, 2, 2>::zero();
    /// m.as_array_mut()[0][0] = 42.0;
    /// assert_eq!(m[(0, 0)], 42.0);
    /// ```
    #[inline]
    pub fn as_array_mut(&mut self) -> &mut [[T; COLS]; ROWS] {
        &mut self.data
    }

    /// Returns a specific row as a fixed-size array.
    ///
    /// # Panics
    ///
    /// Panics if `row >= ROWS`.
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
    /// assert_eq!(m.row(1), &[4.0, 5.0, 6.0]);
    /// ```
    #[inline]
    pub fn row(&self, row: usize) -> &[T; COLS] {
        &self.data[row]
    }

    /// Returns a mutable reference to a specific row.
    ///
    /// # Panics
    ///
    /// Panics if `row >= ROWS`.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::matrix::Matrix;
    ///
    /// let mut m = Matrix::<f64, 2, 3>::zero();
    /// m.row_mut(0)[1] = 5.0;
    /// assert_eq!(m[(0, 1)], 5.0);
    /// ```
    #[inline]
    pub fn row_mut(&mut self, row: usize) -> &mut [T; COLS] {
        &mut self.data[row]
    }

    /// Extracts a specific column as a vector.
    ///
    /// # Panics
    ///
    /// Panics if `col >= COLS`.
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
    /// assert_eq!(m.column(1), [2.0, 5.0]);
    /// ```
    pub fn column(&self, col: usize) -> [T; ROWS] {
        let mut result = [T::zero(); ROWS];
        for row in 0..ROWS {
            result[row] = self.data[row][col];
        }
        result
    }

    /// Applies a function to each element of the matrix.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::matrix::Matrix;
    ///
    /// let m = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
    /// let doubled = m.map(|x| x * 2.0);
    /// assert_eq!(doubled[(0, 0)], 2.0);
    /// assert_eq!(doubled[(1, 1)], 8.0);
    /// ```
    pub fn map<F, U: Scalar>(self, mut f: F) -> Matrix<U, ROWS, COLS>
    where
        F: FnMut(T) -> U,
    {
        let mut result = Matrix::<U, ROWS, COLS>::zero();
        for row in 0..ROWS {
            for col in 0..COLS {
                result.data[row][col] = f(self.data[row][col]);
            }
        }
        result
    }

    /// Applies a function element-wise to two matrices.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::matrix::Matrix;
    ///
    /// let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
    /// let b = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);
    /// let sum = a.zip_with(b, |x, y| x + y);
    /// assert_eq!(sum[(0, 0)], 6.0);
    /// assert_eq!(sum[(1, 1)], 12.0);
    /// ```
    pub fn zip_with<F, U: Scalar>(self, other: Matrix<T, ROWS, COLS>, mut f: F) -> Matrix<U, ROWS, COLS>
    where
        F: FnMut(T, T) -> U,
    {
        let mut result = Matrix::<U, ROWS, COLS>::zero();
        for row in 0..ROWS {
            for col in 0..COLS {
                result.data[row][col] = f(self.data[row][col], other.data[row][col]);
            }
        }
        result
    }

    /// Checks if all elements are zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::matrix::Matrix;
    ///
    /// let zero = Matrix::<f64, 2, 2>::zero();
    /// assert!(zero.is_zero());
    ///
    /// let non_zero = Matrix::<f64, 2, 2>::identity();
    /// assert!(!non_zero.is_zero());
    /// ```
    pub fn is_zero(&self) -> bool {
        for row in 0..ROWS {
            for col in 0..COLS {
                if self.data[row][col] != T::zero() {
                    return false;
                }
            }
        }
        true
    }
}

// Indexing with (row, col) tuple
impl<T: Scalar, const ROWS: usize, const COLS: usize> Index<(usize, usize)> for Matrix<T, ROWS, COLS> {
    type Output = T;

    /// Access element at (row, col).
    ///
    /// # Panics
    ///
    /// Panics if row >= ROWS or col >= COLS.
    #[inline]
    fn index(&self, (row, col): (usize, usize)) -> &Self::Output {
        &self.data[row][col]
    }
}

impl<T: Scalar, const ROWS: usize, const COLS: usize> IndexMut<(usize, usize)> for Matrix<T, ROWS, COLS> {
    /// Mutably access element at (row, col).
    ///
    /// # Panics
    ///
    /// Panics if row >= ROWS or col >= COLS.
    #[inline]
    fn index_mut(&mut self, (row, col): (usize, usize)) -> &mut Self::Output {
        &mut self.data[row][col]
    }
}

// Default implementation
impl<T: Scalar, const ROWS: usize, const COLS: usize> Default for Matrix<T, ROWS, COLS> {
    fn default() -> Self {
        Self::zero()
    }
}

// Debug formatting
impl<T: Scalar + fmt::Debug, const ROWS: usize, const COLS: usize> fmt::Debug for Matrix<T, ROWS, COLS> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Matrix<{}, {}>", ROWS, COLS)?;
        f.debug_list().entries(self.data.iter()).finish()
    }
}

// Display formatting
impl<T: Scalar + fmt::Display, const ROWS: usize, const COLS: usize> fmt::Display for Matrix<T, ROWS, COLS> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "[")?;
        for row in 0..ROWS {
            write!(f, "  [")?;
            for col in 0..COLS {
                if col > 0 {
                    write!(f, ", ")?;
                }
                write!(f, "{}", self.data[row][col])?;
            }
            writeln!(f, "]")?;
        }
        write!(f, "]")
    }
}

// From array conversion
impl<T: Scalar, const ROWS: usize, const COLS: usize> From<[[T; COLS]; ROWS]> for Matrix<T, ROWS, COLS> {
    fn from(data: [[T; COLS]; ROWS]) -> Self {
        Self { data }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dimensions() {
        let m = Matrix::<f64, 3, 4>::zero();
        assert_eq!(m.rows(), 3);
        assert_eq!(m.cols(), 4);
        assert_eq!(m.dimensions(), (3, 4));
    }

    #[test]
    fn test_indexing() {
        let m = Matrix::<f64, 2, 3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        assert_eq!(m[(0, 0)], 1.0);
        assert_eq!(m[(0, 1)], 2.0);
        assert_eq!(m[(1, 2)], 6.0);
    }

    #[test]
    fn test_indexing_mut() {
        let mut m = Matrix::<f64, 2, 2>::zero();
        m[(0, 1)] = 42.0;
        m[(1, 0)] = 7.0;
        assert_eq!(m[(0, 1)], 42.0);
        assert_eq!(m[(1, 0)], 7.0);
    }

    #[test]
    fn test_row_access() {
        let m = Matrix::<f64, 2, 3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        assert_eq!(m.row(0), &[1.0, 2.0, 3.0]);
        assert_eq!(m.row(1), &[4.0, 5.0, 6.0]);
    }

    #[test]
    fn test_row_mut() {
        let mut m = Matrix::<f64, 2, 3>::zero();
        m.row_mut(1)[2] = 99.0;
        assert_eq!(m[(1, 2)], 99.0);
    }

    #[test]
    fn test_column_access() {
        let m = Matrix::<f64, 2, 3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        assert_eq!(m.column(0), [1.0, 4.0]);
        assert_eq!(m.column(1), [2.0, 5.0]);
        assert_eq!(m.column(2), [3.0, 6.0]);
    }

    #[test]
    fn test_map() {
        let m = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let doubled = m.map(|x| x * 2.0);
        assert_eq!(doubled[(0, 0)], 2.0);
        assert_eq!(doubled[(0, 1)], 4.0);
        assert_eq!(doubled[(1, 0)], 6.0);
        assert_eq!(doubled[(1, 1)], 8.0);
    }

    #[test]
    fn test_zip_with() {
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let b = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);
        let sum = a.zip_with(b, |x, y| x + y);
        assert_eq!(sum[(0, 0)], 6.0);
        assert_eq!(sum[(0, 1)], 8.0);
        assert_eq!(sum[(1, 0)], 10.0);
        assert_eq!(sum[(1, 1)], 12.0);
    }

    #[test]
    fn test_is_zero() {
        let zero = Matrix::<f64, 3, 3>::zero();
        assert!(zero.is_zero());

        let mut non_zero = Matrix::<f64, 3, 3>::zero();
        non_zero[(1, 1)] = 1.0;
        assert!(!non_zero.is_zero());
    }

    #[test]
    fn test_from_array() {
        let data = [[1.0, 2.0], [3.0, 4.0]];
        let m = Matrix::from(data);
        assert_eq!(m[(0, 0)], 1.0);
        assert_eq!(m[(1, 1)], 4.0);
    }

    #[test]
    fn test_as_array() {
        let m = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let arr = m.as_array();
        assert_eq!(arr[0][1], 2.0);
        assert_eq!(arr[1][0], 3.0);
    }

    #[test]
    fn test_default() {
        let m: Matrix<f64, 2, 2> = Default::default();
        assert!(m.is_zero());
    }

    #[test]
    fn test_clone_and_copy() {
        let m1 = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let m2 = m1; // Copy
        let m3 = m1.clone();
        assert_eq!(m1, m2);
        assert_eq!(m1, m3);
    }

    #[test]
    fn test_equality() {
        let m1 = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let m2 = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let m3 = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 5.0]]);

        assert_eq!(m1, m2);
        assert_ne!(m1, m3);
    }
}
