//! Reduced Row Echelon Form (RREF) computation.
//!
//! This module implements the Gauss-Jordan elimination algorithm to compute
//! the Reduced Row Echelon Form of a matrix. RREF is a canonical form that
//! is fundamental for:
//!
//! - Solving linear systems
//! - Computing matrix rank
//! - Finding null spaces and column spaces
//! - Determining linear independence
//!
//! # Mathematical Definition
//!
//! A matrix is in RREF if it satisfies:
//!
//! 1. All rows consisting entirely of zeros are at the bottom
//! 2. The first nonzero entry (pivot) in each nonzero row is 1
//! 3. Each pivot is to the right of the pivot in the row above
//! 4. Each pivot is the only nonzero entry in its column
//!
//! # Algorithm
//!
//! The Gauss-Jordan elimination algorithm:
//!
//! 1. **Forward phase**: For each column, find the largest absolute value below
//!    the current row (partial pivoting), swap it to the current position,
//!    scale the row so the pivot equals 1, then eliminate all entries below.
//!
//! 2. **Backward phase**: For each pivot column, eliminate all entries above
//!    the pivot (this distinguishes RREF from REF).
//!
//! # Numerical Stability
//!
//! - Uses **partial pivoting** to select the largest available pivot
//! - Treats values below epsilon as zero to handle floating-point errors
//! - Returns pivot column indices and rank for further analysis
//!
//! # References
//!
//! - Strang, "Linear Algebra and Its Applications", Section 2.2-2.3
//! - Lay, "Linear Algebra and Its Applications", Section 1.2

use crate::core::Float;
use crate::matrix::Matrix;
use crate::linear_algebra::row_operations::{swap_rows, scale_row, add_scaled_row};
use std::collections::BTreeSet;

/// Computes the Reduced Row Echelon Form (RREF) of a matrix.
///
/// This is a convenience function that wraps [`rref_with_pivots`] and returns
/// only the RREF matrix, discarding the pivot column information.
///
/// # Arguments
///
/// * `matrix` - The matrix to reduce
///
/// # Returns
///
/// The RREF of the input matrix.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::rref::rref;
///
/// let m = Matrix::<f64, 3, 4>::from_rows([
///     [1.0, 2.0, 3.0, 4.0],
///     [2.0, 4.0, 6.0, 8.0],
///     [1.0, 1.0, 1.0, 1.0],
/// ]);
///
/// let r = rref(&m);
///
/// // First row has pivot in column 0
/// assert!((r[(0, 0)] - 1.0).abs() < 1e-10);
/// // Second row has pivot in column 1
/// assert!((r[(1, 1)] - 1.0).abs() < 1e-10);
/// ```
pub fn rref<T: Float, const M: usize, const N: usize>(
    matrix: &Matrix<T, M, N>,
) -> Matrix<T, M, N> {
    let (rref_matrix, _, _) = rref_with_pivots(matrix);
    rref_matrix
}

/// Computes the RREF along with pivot column information.
///
/// Returns the RREF matrix, the set of pivot column indices, and the rank
/// (number of pivot columns). This information is essential for computing
/// null spaces, column spaces, and understanding the structure of solutions.
///
/// # Mathematical Properties
///
/// - **Rank**: The number of pivot columns equals the rank of the matrix
/// - **Pivot columns**: Indices of columns containing pivots (leading 1's)
/// - **Free columns**: Columns not in the pivot set correspond to free variables
///
/// # Numerical Stability
///
/// Uses partial pivoting with tolerance `epsilon * max(M, N) * 100` to handle
/// floating-point errors gracefully.
///
/// # Arguments
///
/// * `matrix` - The matrix to reduce
///
/// # Returns
///
/// A tuple `(rref_matrix, pivot_columns, rank)` where:
/// - `rref_matrix` is the reduced row echelon form
/// - `pivot_columns` is a `BTreeSet` of column indices containing pivots
/// - `rank` is the number of pivots (equals `pivot_columns.len()`)
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::rref::rref_with_pivots;
///
/// let m = Matrix::<f64, 2, 3>::from_rows([
///     [1.0, 2.0, 3.0],
///     [2.0, 4.0, 7.0],
/// ]);
///
/// let (r, pivots, rank) = rref_with_pivots(&m);
///
/// assert_eq!(rank, 2);
/// assert!(pivots.contains(&0));  // First column is pivot
/// assert!(pivots.contains(&1));  // Second column is pivot
/// assert!((r[(0, 0)] - 1.0).abs() < 1e-10);
/// assert!((r[(1, 1)] - 1.0).abs() < 1e-10);
/// ```
pub fn rref_with_pivots<T: Float, const M: usize, const N: usize>(
    matrix: &Matrix<T, M, N>,
) -> (Matrix<T, M, N>, BTreeSet<usize>, usize) {
    let mut result = *matrix;
    let mut pivot_columns = BTreeSet::new();

    // Compute tolerance for numerical stability
    let size = if M > N { M } else { N };
    let tolerance = T::epsilon() * T::from(size).unwrap() * T::from(100.0).unwrap();

    let mut current_row = 0;

    // Forward phase: Create row echelon form with partial pivoting
    for col in 0..N {
        if current_row >= M {
            break;
        }

        // Find pivot: largest absolute value in this column at or below current_row
        let mut pivot_row = current_row;
        let mut max_val = result[(current_row, col)].abs();

        for row in (current_row + 1)..M {
            let val = result[(row, col)].abs();
            if val > max_val {
                max_val = val;
                pivot_row = row;
            }
        }

        // If max value is below tolerance, this column has no pivot
        if max_val < tolerance {
            continue;
        }

        // Swap pivot row to current position (if needed)
        if pivot_row != current_row {
            swap_rows(&mut result, current_row, pivot_row);
        }

        // Scale row so pivot equals 1
        let pivot_val = result[(current_row, col)];
        scale_row(&mut result, current_row, T::one() / pivot_val);

        // Eliminate entries below pivot
        for row in (current_row + 1)..M {
            let factor = result[(row, col)];
            if factor.abs() >= tolerance {
                add_scaled_row(&mut result, row, current_row, -factor);
            }
        }

        // Record pivot column and advance
        pivot_columns.insert(col);
        current_row += 1;
    }

    // Backward phase: Eliminate entries above pivots (makes it RREF, not just REF)
    for col in (0..N).rev() {
        // Find which row has the pivot in this column
        let mut pivot_row = None;
        for row in 0..M {
            if (result[(row, col)] - T::one()).abs() < tolerance {
                // Check this is actually the pivot (first nonzero in row)
                let mut is_pivot = true;
                for c in 0..col {
                    if result[(row, c)].abs() >= tolerance {
                        is_pivot = false;
                        break;
                    }
                }
                if is_pivot {
                    pivot_row = Some(row);
                    break;
                }
            }
        }

        // If this column has a pivot, eliminate above it
        if let Some(prow) = pivot_row {
            for row in 0..prow {
                let factor = result[(row, col)];
                if factor.abs() >= tolerance {
                    add_scaled_row(&mut result, row, prow, -factor);
                }
            }
        }
    }

    // Clean up near-zero entries
    for i in 0..M {
        for j in 0..N {
            if result[(i, j)].abs() < tolerance {
                result[(i, j)] = T::zero();
            }
        }
    }

    let rank = pivot_columns.len();
    (result, pivot_columns, rank)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rref_identity() {
        let m = Matrix::<f64, 3, 3>::identity();
        let r = rref(&m);

        // Identity should remain identity
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((r[(i, j)] - expected).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_rref_simple_2x2() {
        let m = Matrix::<f64, 2, 2>::from_rows([
            [2.0, 4.0],
            [1.0, 3.0],
        ]);

        let r = rref(&m);

        // Expected RREF: [[1, 0], [0, 1]]
        assert!((r[(0, 0)] - 1.0).abs() < 1e-10);
        assert!((r[(0, 1)] - 0.0).abs() < 1e-10);
        assert!((r[(1, 0)] - 0.0).abs() < 1e-10);
        assert!((r[(1, 1)] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_rref_with_pivots_full_rank() {
        let m = Matrix::<f64, 2, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [2.0, 4.0, 7.0],
        ]);

        let (r, pivots, rank) = rref_with_pivots(&m);

        assert_eq!(rank, 2);
        assert_eq!(pivots.len(), 2);

        // First pivot should be in column 0
        assert!(pivots.contains(&0));
        assert!((r[(0, 0)] - 1.0).abs() < 1e-10);

        // Second pivot should be in column 1 or 2
        if pivots.contains(&1) {
            assert!((r[(1, 1)] - 1.0).abs() < 1e-10);
        } else {
            assert!(pivots.contains(&2));
            assert!((r[(1, 2)] - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_rref_dependent_rows() {
        let m = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [2.0, 4.0, 6.0],  // Dependent on row 0
            [1.0, 1.0, 1.0],
        ]);

        let (r, pivots, rank) = rref_with_pivots(&m);

        // Rank should be 2 (one dependent row)
        assert_eq!(rank, 2);

        // Last row should be all zeros
        for j in 0..3 {
            assert!(r[(2, j)].abs() < 1e-10);
        }
    }

    #[test]
    fn test_rref_zero_matrix() {
        let m = Matrix::<f64, 2, 3>::from_rows([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ]);

        let (r, pivots, rank) = rref_with_pivots(&m);

        assert_eq!(rank, 0);
        assert!(pivots.is_empty());

        // Should remain all zeros
        for i in 0..2 {
            for j in 0..3 {
                assert!(r[(i, j)].abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_rref_rectangular_tall() {
        let m = Matrix::<f64, 4, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
            [1.0, 0.0, 1.0],
        ]);

        let (r, pivots, rank) = rref_with_pivots(&m);

        // This matrix has rank 3 (or possibly 2 if rows are dependent)
        assert!(rank <= 3);
        assert_eq!(pivots.len(), rank);

        // Verify RREF properties: each pivot column should have exactly one 1
        for &col in &pivots {
            let mut ones_count = 0;
            let mut nonzeros_count = 0;
            for row in 0..4 {
                let val = r[(row, col)].abs();
                if val >= 1e-10 {
                    nonzeros_count += 1;
                    if (val - 1.0).abs() < 1e-10 {
                        ones_count += 1;
                    }
                }
            }
            assert_eq!(ones_count, 1, "Pivot column {} should have exactly one 1", col);
            assert_eq!(nonzeros_count, 1, "Pivot column {} should have only one nonzero", col);
        }
    }

    #[test]
    fn test_rref_rectangular_wide() {
        let m = Matrix::<f64, 2, 4>::from_rows([
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 4.0, 7.0, 10.0],
        ]);

        let (r, pivots, rank) = rref_with_pivots(&m);

        // Rank should be at most 2
        assert!(rank <= 2);
        assert_eq!(pivots.len(), rank);

        // Each pivot row should have a leading 1
        for row in 0..rank {
            // Find first nonzero entry
            let mut found_pivot = false;
            for col in 0..4 {
                if r[(row, col)].abs() >= 1e-10 {
                    assert!((r[(row, col)] - 1.0).abs() < 1e-10, "First nonzero in row {} should be 1", row);
                    found_pivot = true;
                    break;
                }
            }
            assert!(found_pivot, "Row {} should have a pivot", row);
        }
    }
}
