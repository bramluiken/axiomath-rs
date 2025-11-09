//! Elementary row operations for matrix manipulation.
//!
//! This module provides the three elementary row operations used in Gaussian elimination,
//! Gauss-Jordan elimination, and other matrix algorithms:
//!
//! 1. **Row Swap** - Exchange two rows
//! 2. **Row Scale** - Multiply a row by a nonzero scalar
//! 3. **Row Addition** - Add a multiple of one row to another
//!
//! These operations are fundamental building blocks for:
//! - Computing Reduced Row Echelon Form (RREF)
//! - Gaussian elimination
//! - Matrix inversion via augmented matrices
//! - Finding null spaces
//!
//! # Mathematical Properties
//!
//! Elementary row operations preserve:
//! - The row space of the matrix
//! - The solution set of the linear system Ax = b
//! - The rank of the matrix
//!
//! # References
//!
//! - Strang, "Linear Algebra and Its Applications", Section 2.2
//! - Lay, "Linear Algebra and Its Applications", Section 1.1

use crate::core::Scalar;
use crate::matrix::Matrix;

/// Swaps two rows in a matrix.
///
/// Exchanges row `row1` with row `row2`. If `row1 == row2`, the operation
/// is a no-op (no change).
///
/// # Mathematical Notation
///
/// ```text
/// R_i ↔ R_j
/// ```
///
/// # Arguments
///
/// * `matrix` - Mutable reference to the matrix
/// * `row1` - Index of first row (0-indexed)
/// * `row2` - Index of second row (0-indexed)
///
/// # Panics
///
/// Panics if `row1 >= M` or `row2 >= M`.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::row_operations::swap_rows;
///
/// let mut m = Matrix::<f64, 3, 2>::from_rows([
///     [1.0, 2.0],
///     [3.0, 4.0],
///     [5.0, 6.0],
/// ]);
///
/// swap_rows(&mut m, 0, 2);
///
/// assert_eq!(m[(0, 0)], 5.0);
/// assert_eq!(m[(2, 0)], 1.0);
/// ```
pub fn swap_rows<T: Scalar, const M: usize, const N: usize>(
    matrix: &mut Matrix<T, M, N>,
    row1: usize,
    row2: usize,
) {
    if row1 == row2 {
        return;
    }

    assert!(row1 < M, "row1 index {} out of bounds for matrix with {} rows", row1, M);
    assert!(row2 < M, "row2 index {} out of bounds for matrix with {} rows", row2, M);

    for j in 0..N {
        let temp = matrix[(row1, j)];
        matrix[(row1, j)] = matrix[(row2, j)];
        matrix[(row2, j)] = temp;
    }
}

/// Scales a row by a scalar.
///
/// Multiplies every element in row `row` by `scalar`.
///
/// # Mathematical Notation
///
/// ```text
/// R_i → c · R_i
/// ```
///
/// # Arguments
///
/// * `matrix` - Mutable reference to the matrix
/// * `row` - Index of the row to scale (0-indexed)
/// * `scalar` - Scalar multiplier
///
/// # Panics
///
/// Panics if `row >= M`.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::row_operations::scale_row;
///
/// let mut m = Matrix::<f64, 2, 2>::from_rows([
///     [2.0, 4.0],
///     [1.0, 3.0],
/// ]);
///
/// scale_row(&mut m, 0, 0.5);
///
/// assert_eq!(m[(0, 0)], 1.0);
/// assert_eq!(m[(0, 1)], 2.0);
/// ```
pub fn scale_row<T: Scalar, const M: usize, const N: usize>(
    matrix: &mut Matrix<T, M, N>,
    row: usize,
    scalar: T,
) {
    assert!(row < M, "row index {} out of bounds for matrix with {} rows", row, M);

    for j in 0..N {
        matrix[(row, j)] = matrix[(row, j)] * scalar;
    }
}

/// Adds a scaled row to another row.
///
/// Replaces `target_row` with `target_row + scalar * source_row`.
/// The source row remains unchanged.
///
/// # Mathematical Notation
///
/// ```text
/// R_i → R_i + c · R_j
/// ```
///
/// # Arguments
///
/// * `matrix` - Mutable reference to the matrix
/// * `target_row` - Index of the row to be modified (0-indexed)
/// * `source_row` - Index of the row to add from (0-indexed)
/// * `scalar` - Scalar multiplier for the source row
///
/// # Panics
///
/// Panics if `target_row >= M` or `source_row >= M`.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::row_operations::add_scaled_row;
///
/// let mut m = Matrix::<f64, 2, 2>::from_rows([
///     [1.0, 2.0],
///     [3.0, 4.0],
/// ]);
///
/// // Eliminate first column of row 1: R₁ → R₁ - 3R₀
/// add_scaled_row(&mut m, 1, 0, -3.0);
///
/// assert_eq!(m[(1, 0)], 0.0);  // 3 - 3*1 = 0
/// assert_eq!(m[(1, 1)], -2.0); // 4 - 3*2 = -2
/// ```
pub fn add_scaled_row<T: Scalar, const M: usize, const N: usize>(
    matrix: &mut Matrix<T, M, N>,
    target_row: usize,
    source_row: usize,
    scalar: T,
) {
    assert!(target_row < M, "target_row index {} out of bounds for matrix with {} rows", target_row, M);
    assert!(source_row < M, "source_row index {} out of bounds for matrix with {} rows", source_row, M);

    for j in 0..N {
        matrix[(target_row, j)] = matrix[(target_row, j)] + scalar * matrix[(source_row, j)];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_swap_rows() {
        let mut m = Matrix::<f64, 3, 2>::from_rows([
            [1.0, 2.0],
            [3.0, 4.0],
            [5.0, 6.0],
        ]);

        swap_rows(&mut m, 0, 2);

        // Row 0 should now have values from original row 2
        assert_eq!(m[(0, 0)], 5.0);
        assert_eq!(m[(0, 1)], 6.0);

        // Row 2 should now have values from original row 0
        assert_eq!(m[(2, 0)], 1.0);
        assert_eq!(m[(2, 1)], 2.0);

        // Row 1 should be unchanged
        assert_eq!(m[(1, 0)], 3.0);
        assert_eq!(m[(1, 1)], 4.0);
    }

    #[test]
    fn test_swap_rows_same() {
        let original = Matrix::<f64, 2, 2>::from_rows([
            [1.0, 2.0],
            [3.0, 4.0],
        ]);

        let mut m = original;
        swap_rows(&mut m, 1, 1); // Swap row with itself

        // Matrix should be unchanged
        for i in 0..2 {
            for j in 0..2 {
                assert_eq!(m[(i, j)], original[(i, j)]);
            }
        }
    }

    #[test]
    fn test_scale_row() {
        let mut m = Matrix::<f64, 2, 2>::from_rows([
            [2.0, 4.0],
            [1.0, 3.0],
        ]);

        scale_row(&mut m, 0, 0.5);

        // Row 0 should be scaled by 0.5
        assert_eq!(m[(0, 0)], 1.0);
        assert_eq!(m[(0, 1)], 2.0);

        // Row 1 should be unchanged
        assert_eq!(m[(1, 0)], 1.0);
        assert_eq!(m[(1, 1)], 3.0);
    }

    #[test]
    fn test_scale_row_by_negative() {
        let mut m = Matrix::<f64, 2, 2>::from_rows([
            [1.0, 2.0],
            [3.0, 4.0],
        ]);

        scale_row(&mut m, 1, -2.0);

        // Row 1 should be scaled by -2
        assert_eq!(m[(1, 0)], -6.0);
        assert_eq!(m[(1, 1)], -8.0);

        // Row 0 unchanged
        assert_eq!(m[(0, 0)], 1.0);
        assert_eq!(m[(0, 1)], 2.0);
    }

    #[test]
    fn test_add_scaled_row() {
        let mut m = Matrix::<f64, 2, 2>::from_rows([
            [1.0, 2.0],
            [3.0, 4.0],
        ]);

        // Eliminate first column of row 1: R₁ → R₁ - 3R₀
        add_scaled_row(&mut m, 1, 0, -3.0);

        // Row 1 should be modified
        assert_eq!(m[(1, 0)], 0.0);  // 3 - 3*1 = 0
        assert_eq!(m[(1, 1)], -2.0); // 4 - 3*2 = -2

        // Row 0 should be unchanged (source row is not modified)
        assert_eq!(m[(0, 0)], 1.0);
        assert_eq!(m[(0, 1)], 2.0);
    }

    #[test]
    fn test_add_scaled_row_zero_scalar() {
        let original = Matrix::<f64, 2, 2>::from_rows([
            [1.0, 2.0],
            [3.0, 4.0],
        ]);

        let mut m = original;
        add_scaled_row(&mut m, 1, 0, 0.0);

        // Matrix should be unchanged when scalar is 0
        for i in 0..2 {
            for j in 0..2 {
                assert_eq!(m[(i, j)], original[(i, j)]);
            }
        }
    }

    #[test]
    fn test_row_operations_sequence() {
        // Test a sequence of operations like in Gaussian elimination
        let mut m = Matrix::<f64, 3, 3>::from_rows([
            [2.0, 4.0, 6.0],
            [1.0, 3.0, 5.0],
            [3.0, 7.0, 9.0],
        ]);

        // Step 1: Scale row 0 to make pivot = 1
        scale_row(&mut m, 0, 0.5);
        assert_eq!(m[(0, 0)], 1.0);

        // Step 2: Eliminate column 0 below pivot
        add_scaled_row(&mut m, 1, 0, -1.0);
        add_scaled_row(&mut m, 2, 0, -3.0);

        assert_eq!(m[(1, 0)], 0.0);
        assert_eq!(m[(2, 0)], 0.0);

        // Step 3: Swap rows 1 and 2 if needed (test swap in sequence)
        if m[(1, 1)].abs() < m[(2, 1)].abs() {
            swap_rows(&mut m, 1, 2);
        }

        // Verify row 0 is still correct after all operations
        assert_eq!(m[(0, 0)], 1.0);
        assert_eq!(m[(0, 1)], 2.0);
        assert_eq!(m[(0, 2)], 3.0);
    }
}
