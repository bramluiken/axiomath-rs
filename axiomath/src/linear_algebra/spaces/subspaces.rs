//! Subspace operations for matrices.
//!
//! This module provides functions for computing fundamental subspaces
//! associated with matrices, including:
//! - Column space (range)
//! - Null space (kernel)
//! - Row space
//! - Left null space
//!
//! # Mathematical Background
//!
//! For a matrix A (m×n), the fundamental subspaces are:
//!
//! ## Column Space (Range)
//! ```text
//! Col(A) = {Ax : x ∈ ℝⁿ}
//! ```
//! The span of the columns of A. Dimension = rank(A).
//!
//! ## Null Space (Kernel)
//! ```text
//! Null(A) = {x ∈ ℝⁿ : Ax = 0}
//! ```
//! All vectors that map to zero. Dimension = n - rank(A).
//!
//! ## Row Space
//! ```text
//! Row(A) = Col(Aᵀ)
//! ```
//! The span of the rows of A. Dimension = rank(A).
//!
//! ## Left Null Space
//! ```text
//! Null(Aᵀ) = {y ∈ ℝᵐ : Aᵀy = 0}
//! ```
//! Dimension = m - rank(A).
//!
//! # Fundamental Theorem of Linear Algebra
//!
//! These four subspaces form two orthogonal pairs:
//! - Row(A) ⊥ Null(A) in ℝⁿ
//! - Col(A) ⊥ Null(Aᵀ) in ℝᵐ
//!
//! # References
//!
//! - Strang, "Linear Algebra and Its Applications", Chapter 3
//! - Strang, "The Fundamental Theorem of Linear Algebra"

use crate::core::traits::Float;
use crate::linear_algebra::spaces::gram_schmidt;
use crate::linear_algebra::rref::rref_with_pivots;
use crate::linear_algebra::LinearAlgebraError;
use crate::matrix::{transpose, Matrix};
use crate::vector::Vector;

/// Computes a basis for the column space of a matrix.
///
/// # Mathematical Definition
///
/// The column space (or range) of a matrix A is:
/// ```text
/// Col(A) = span{a₁, a₂, ..., aₙ}
/// ```
/// where aᵢ are the columns of A.
///
/// # Algorithm
///
/// 1. Extract columns from the matrix
/// 2. Apply Gram-Schmidt to find an orthonormal basis
/// 3. The resulting independent vectors form a basis
///
/// # Complexity
///
/// - Time: O(m²n) where m is rows, n is columns
/// - Space: O(mn)
///
/// # Arguments
///
/// * `matrix` - Matrix to compute column space for (m×n)
///
/// # Returns
///
/// A vector of linearly independent vectors that span Col(A).
/// The number of vectors equals rank(A).
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::spaces::column_space;
///
/// let a = Matrix::from_rows([
///     [1.0, 2.0, 3.0],
///     [4.0, 5.0, 6.0],
///     [7.0, 8.0, 9.0],
/// ]);
///
/// let col_basis = column_space(&a).unwrap();
/// // col_basis contains 2 vectors (rank = 2)
/// ```
pub fn column_space<T: Float, const M: usize, const N: usize>(
    matrix: &Matrix<T, M, N>,
) -> Result<Vec<Vector<T, M>>, LinearAlgebraError> {
    // Extract columns
    let mut columns = Vec::with_capacity(N);
    for j in 0..N {
        let mut col = Vector::<T, M>::zero();
        for i in 0..M {
            col[i] = matrix[(i, j)];
        }
        columns.push(col);
    }

    // Apply Gram-Schmidt to get orthonormal basis
    gram_schmidt(&columns)
}

/// Computes a basis for the row space of a matrix.
///
/// # Mathematical Definition
///
/// The row space of A is the column space of Aᵀ:
/// ```text
/// Row(A) = Col(Aᵀ) = span{r₁, r₂, ..., rₘ}
/// ```
/// where rᵢ are the rows of A.
///
/// # Algorithm
///
/// 1. Transpose the matrix
/// 2. Compute column space of the transpose
///
/// # Complexity
///
/// - Time: O(mn²)
/// - Space: O(mn)
///
/// # Arguments
///
/// * `matrix` - Matrix to compute row space for (m×n)
///
/// # Returns
///
/// A vector of linearly independent vectors that span Row(A).
/// The number of vectors equals rank(A).
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::spaces::row_space;
///
/// let a = Matrix::from_rows([
///     [1.0, 2.0],
///     [3.0, 4.0],
///     [5.0, 6.0],
/// ]);
///
/// let row_basis = row_space(&a).unwrap();
/// // row_basis contains rank(A) vectors
/// ```
pub fn row_space<T: Float, const M: usize, const N: usize>(
    matrix: &Matrix<T, M, N>,
) -> Result<Vec<Vector<T, N>>, LinearAlgebraError> {
    // Row space of A = column space of A^T
    let at = transpose(matrix);
    column_space(&at)
}

/// Computes a basis for the null space of a matrix.
///
/// # Mathematical Definition
///
/// The null space (or kernel) of a matrix A is:
/// ```text
/// Null(A) = {x ∈ ℝⁿ : Ax = 0}
/// ```
///
/// The dimension of the null space (nullity) is given by:
/// ```text
/// nullity(A) = n - rank(A)
/// ```
///
/// # Algorithm
///
/// 1. Compute the RREF (Reduced Row Echelon Form) of A
/// 2. Identify pivot columns and free columns
///    - Pivot columns contain the leading 1's in the RREF
///    - Free columns are all other columns
/// 3. For each free column j:
///    - Set x_j = 1 (the free variable)
///    - Set other free variables to 0
///    - Read off pivot variables from the RREF to satisfy Ax = 0
///    - This gives one basis vector
///
/// # Complexity
///
/// - Time: O(m²n) for RREF computation
/// - Space: O(mn)
///
/// # Arguments
///
/// * `matrix` - Matrix to compute null space for (m×n)
///
/// # Returns
///
/// A vector of linearly independent vectors that span Null(A).
/// The number of vectors equals n - rank(A) (nullity).
///
/// If the matrix has full column rank (rank = n), returns an empty vector.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::spaces::null_space;
/// use axiomath::matrix::matrix_vector_multiply;
/// use axiomath::vector::magnitude;
///
/// // Matrix with a null space
/// let a = Matrix::from_rows([
///     [1.0, 2.0, 3.0],
///     [4.0, 5.0, 6.0],
///     [7.0, 8.0, 9.0],
/// ]);
///
/// let null_basis = null_space(&a).unwrap();
///
/// // Verify each basis vector satisfies Ax = 0
/// for v in &null_basis {
///     let av = matrix_vector_multiply(&a, v);
///     let norm = magnitude(&av);
///     assert!(norm < 1e-10);  // av ≈ 0
/// }
/// ```
pub fn null_space<T: Float, const M: usize, const N: usize>(
    matrix: &Matrix<T, M, N>,
) -> Result<Vec<Vector<T, N>>, LinearAlgebraError> {
    // Compute RREF and get pivot columns
    let (rref, pivot_columns, rank) = rref_with_pivots(matrix);

    // Nullity = n - rank
    let nullity = N - rank;

    // If full column rank, null space is trivial (only zero vector)
    if nullity == 0 {
        return Ok(Vec::new());
    }

    // Identify free columns (not in pivot set)
    let mut free_columns = Vec::new();
    for j in 0..N {
        if !pivot_columns.contains(&j) {
            free_columns.push(j);
        }
    }

    assert_eq!(
        free_columns.len(),
        nullity,
        "Number of free columns should equal nullity"
    );

    // Build basis vectors
    let mut basis = Vec::with_capacity(nullity);
    let tolerance = T::epsilon() * T::from(N).unwrap() * T::from(100.0).unwrap();

    for &free_col in &free_columns {
        // Construct a basis vector with this free variable = 1
        let mut basis_vec = Vector::<T, N>::zero();
        basis_vec[free_col] = T::one();

        // For each pivot column, read the RREF to find what value to assign
        // We need to find which row contains the pivot for each pivot column
        for &pivot_col in &pivot_columns {
            // Find the row with the pivot in this column
            let mut pivot_row = None;
            for i in 0..M {
                if (rref[(i, pivot_col)] - T::one()).abs() < tolerance {
                    pivot_row = Some(i);
                    break;
                }
            }

            if let Some(prow) = pivot_row {
                // The entry in this row at the free column tells us the coefficient
                // In RREF form: x_pivot + c*x_free = 0, so x_pivot = -c*x_free
                // Since x_free = 1, we have x_pivot = -c
                let coeff = rref[(prow, free_col)];
                basis_vec[pivot_col] = -coeff;
            }
        }

        basis.push(basis_vec);
    }

    Ok(basis)
}

/// Computes a basis for the left null space of a matrix.
///
/// # Mathematical Definition
///
/// The left null space of A is the null space of Aᵀ:
/// ```text
/// Null(Aᵀ) = {y ∈ ℝᵐ : Aᵀy = 0} = {y ∈ ℝᵐ : yᵀA = 0}
/// ```
///
/// # Algorithm
///
/// Compute Null(Aᵀ) using the null space algorithm.
///
/// # Status
///
/// Currently unimplemented. Requires null space implementation.
///
/// # Arguments
///
/// * `matrix` - Matrix to compute left null space for (m×n)
///
/// # Returns
///
/// A vector of linearly independent vectors that span Null(Aᵀ).
/// The number of vectors equals m - rank(A).
///
/// # Examples
///
/// ```ignore
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::spaces::left_null_space;
///
/// let a = Matrix::from_rows([
///     [1.0, 2.0],
///     [3.0, 4.0],
///     [5.0, 6.0],
/// ]);
///
/// let left_null_basis = left_null_space(&a).unwrap();
/// // left_null_basis dimension = 3 - rank(A)
/// ```
pub fn left_null_space<T: Float, const M: usize, const N: usize>(
    matrix: &Matrix<T, M, N>,
) -> Result<Vec<Vector<T, M>>, LinearAlgebraError> {
    // Left null space of A = null space of A^T
    let at = transpose(matrix);
    null_space(&at)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::comparison::approximately_equal;
    use crate::vector::{dot_product, magnitude};

    #[test]
    fn test_column_space_full_rank() {
        // Full rank 3×3 matrix
        let a = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);

        let col_basis = column_space(&a).unwrap();

        // Should have 3 basis vectors (full rank)
        assert_eq!(col_basis.len(), 3);

        // All should be orthonormal
        for i in 0..3 {
            let norm = magnitude(&col_basis[i]);
            assert!(
                approximately_equal(norm, 1.0, 1e-10, 1e-10),
                "Basis vector {} has norm {}, expected 1.0",
                i,
                norm
            );

            for j in (i + 1)..3 {
                let dot = dot_product(&col_basis[i], &col_basis[j]);
                assert!(
                    approximately_equal(dot, 0.0, 1e-10, 1e-10),
                    "Basis vectors {} and {} have dot product {}, expected 0",
                    i,
                    j,
                    dot
                );
            }
        }
    }

    #[test]
    fn test_column_space_rank_deficient() {
        // Rank 2 matrix (third column = 2*first - second)
        let a = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 0.0, 2.0],
            [0.0, 1.0, -1.0],
            [0.0, 0.0, 0.0],
        ]);

        let col_basis = column_space(&a).unwrap();

        // Should have 2 basis vectors (rank = 2)
        assert_eq!(col_basis.len(), 2);

        // Should be orthonormal
        for vec in &col_basis {
            let norm = magnitude(vec);
            assert!(approximately_equal(norm, 1.0, 1e-10, 1e-10));
        }

        let dot = dot_product(&col_basis[0], &col_basis[1]);
        assert!(approximately_equal(dot, 0.0, 1e-10, 1e-10));
    }

    #[test]
    fn test_column_space_tall_matrix() {
        // 4×2 matrix
        let a = Matrix::<f64, 4, 2>::from_rows([
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [0.0, 0.0],
        ]);

        let col_basis = column_space(&a).unwrap();

        // Should have 2 basis vectors (full column rank)
        assert_eq!(col_basis.len(), 2);

        // Should be orthonormal
        for vec in &col_basis {
            let norm = magnitude(vec);
            assert!(approximately_equal(norm, 1.0, 1e-10, 1e-10));
        }
    }

    #[test]
    fn test_row_space_full_rank() {
        let a = Matrix::<f64, 3, 3>::identity();

        let row_basis = row_space(&a).unwrap();

        // Should have 3 basis vectors
        assert_eq!(row_basis.len(), 3);

        // Should be orthonormal
        for vec in &row_basis {
            let norm = magnitude(vec);
            assert!(approximately_equal(norm, 1.0, 1e-10, 1e-10));
        }
    }

    #[test]
    fn test_row_space_wide_matrix() {
        // 2×4 matrix with rank 2
        let a = Matrix::<f64, 2, 4>::from_rows([
            [1.0, 0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0, 1.0],
        ]);

        let row_basis = row_space(&a).unwrap();

        // Should have 2 basis vectors (full row rank)
        assert_eq!(row_basis.len(), 2);

        // Should be orthonormal
        for vec in &row_basis {
            let norm = magnitude(vec);
            assert!(approximately_equal(norm, 1.0, 1e-10, 1e-10));
        }

        let dot = dot_product(&row_basis[0], &row_basis[1]);
        assert!(approximately_equal(dot, 0.0, 1e-10, 1e-10));
    }

    #[test]
    fn test_column_and_row_space_same_rank() {
        // Any matrix should have same dimension for column and row space
        let a = Matrix::<f64, 3, 4>::from_rows([
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 4.0, 6.0, 8.0],  // 2× first row
            [0.0, 0.0, 1.0, 1.0],
        ]);

        let col_basis = column_space(&a).unwrap();
        let row_basis = row_space(&a).unwrap();

        // rank(A) should be same from column and row space
        assert_eq!(col_basis.len(), row_basis.len());
    }

    // ========== NULL SPACE TESTS ==========

    #[test]
    fn test_null_space_full_rank() {
        // Full rank square matrix has trivial null space
        let a = Matrix::<f64, 3, 3>::identity();
        let null_basis = null_space(&a).unwrap();

        // Null space should be empty (only zero vector)
        assert_eq!(null_basis.len(), 0);
    }

    #[test]
    fn test_null_space_full_column_rank() {
        // Full column rank tall matrix has trivial null space
        let a = Matrix::<f64, 3, 2>::from_rows([
            [1.0, 0.0],
            [0.0, 1.0],
            [0.0, 0.0],
        ]);

        let null_basis = null_space(&a).unwrap();

        // Null space should be empty
        assert_eq!(null_basis.len(), 0);
    }

    #[test]
    fn test_null_space_rank_one() {
        // Rank 1 matrix: all rows proportional
        let a = Matrix::<f64, 2, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [2.0, 4.0, 6.0],
        ]);

        let null_basis = null_space(&a).unwrap();

        // nullity = 3 - 1 = 2
        assert_eq!(null_basis.len(), 2);

        // Verify each basis vector satisfies Ax = 0
        use crate::matrix::matrix_vector_multiply;
        for v in &null_basis {
            let av = matrix_vector_multiply(&a, v);
            for i in 0..2 {
                assert!(
                    av[i].abs() < 1e-10,
                    "Ax should be zero, but component {} = {}",
                    i,
                    av[i]
                );
            }
        }
    }

    #[test]
    fn test_null_space_simple_3x3() {
        // Matrix with known null space
        // Third row = first + second, so rank = 2, nullity = 1
        let a = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 2.0],
        ]);

        let null_basis = null_space(&a).unwrap();

        // nullity = 3 - 2 = 1
        assert_eq!(null_basis.len(), 1);

        // Verify Ax = 0
        use crate::matrix::matrix_vector_multiply;
        let v = &null_basis[0];
        let av = matrix_vector_multiply(&a, v);

        for i in 0..3 {
            assert!(
                av[i].abs() < 1e-10,
                "Ax should be zero, but component {} = {}",
                i,
                av[i]
            );
        }

        // The null space should be span{[-1, -1, 1]} or a scalar multiple
        // Verify v is proportional to [-1, -1, 1]
        let expected_ratio = if v[2].abs() > 1e-10 {
            v[2] // Use last component as reference
        } else {
            panic!("Null space basis vector has unexpected form");
        };

        let normalized_v = Vector::<f64, 3>::from([
            v[0] / expected_ratio,
            v[1] / expected_ratio,
            v[2] / expected_ratio,
        ]);

        assert!(
            (normalized_v[0] + 1.0).abs() < 1e-10
                || (normalized_v[0] - (-1.0)).abs() < 1e-10
        );
        assert!(
            (normalized_v[1] + 1.0).abs() < 1e-10
                || (normalized_v[1] - (-1.0)).abs() < 1e-10
        );
        assert!((normalized_v[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_null_space_zero_matrix() {
        // Zero matrix has full null space
        let a = Matrix::<f64, 2, 3>::from_rows([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ]);

        let null_basis = null_space(&a).unwrap();

        // nullity = 3 - 0 = 3 (all columns are free)
        assert_eq!(null_basis.len(), 3);

        // Verify all basis vectors are linearly independent
        // (this is automatically guaranteed by the algorithm)
    }

    #[test]
    fn test_null_space_wide_matrix() {
        // Wide matrix (more columns than rows) always has nontrivial null space
        let a = Matrix::<f64, 2, 4>::from_rows([
            [1.0, 0.0, 1.0, 2.0],
            [0.0, 1.0, 1.0, 3.0],
        ]);

        let null_basis = null_space(&a).unwrap();

        // nullity = 4 - rank(A), and rank(A) ≤ min(2, 4) = 2
        // So nullity ≥ 2
        assert!(null_basis.len() >= 2);

        // Verify Ax = 0 for all basis vectors
        use crate::matrix::matrix_vector_multiply;
        for v in &null_basis {
            let av = matrix_vector_multiply(&a, v);
            for i in 0..2 {
                assert!(
                    av[i].abs() < 1e-10,
                    "Ax should be zero, but component {} = {}",
                    i,
                    av[i]
                );
            }
        }
    }

    #[test]
    fn test_null_space_basis_independence() {
        // Test that null space basis vectors are linearly independent
        let a = Matrix::<f64, 2, 5>::from_rows([
            [1.0, 2.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0, 2.0, 1.0],
        ]);

        let null_basis = null_space(&a).unwrap();

        // Should have 5 - 2 = 3 basis vectors
        assert_eq!(null_basis.len(), 3);

        // Verify linear independence using Gram-Schmidt test
        // (if vectors are dependent, Gram-Schmidt will reduce the count)
        let gs_basis = gram_schmidt(&null_basis).unwrap();
        assert_eq!(
            gs_basis.len(),
            null_basis.len(),
            "Null space basis vectors should be linearly independent"
        );
    }

    #[test]
    fn test_null_space_orthogonal_to_row_space() {
        // Fundamental theorem: Null(A) ⊥ Row(A)
        let a = Matrix::<f64, 3, 4>::from_rows([
            [1.0, 2.0, 3.0, 4.0],
            [0.0, 1.0, 2.0, 3.0],
            [0.0, 0.0, 0.0, 0.0],
        ]);

        let null_basis = null_space(&a).unwrap();
        let row_basis = row_space(&a).unwrap();

        // Every null space vector should be orthogonal to every row space vector
        for nv in &null_basis {
            for rv in &row_basis {
                let dot = dot_product(nv, rv);
                assert!(
                    dot.abs() < 1e-10,
                    "Null space and row space should be orthogonal, but dot = {}",
                    dot
                );
            }
        }
    }

    #[test]
    fn test_null_space_dimension_theorem_2x3() {
        // Rank-nullity theorem: rank(A) + nullity(A) = n
        let a = Matrix::<f64, 2, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
        ]);

        let null_basis = null_space(&a).unwrap();
        let col_basis = column_space(&a).unwrap();

        let nullity = null_basis.len();
        let rank = col_basis.len();

        // For 2×3 matrix: rank + nullity = 3
        assert_eq!(rank + nullity, 3);
    }

    #[test]
    fn test_null_space_dimension_theorem_3x4() {
        let a = Matrix::<f64, 3, 4>::from_rows([
            [1.0, 0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0, 1.0],
        ]);

        let null_basis = null_space(&a).unwrap();
        let col_basis = column_space(&a).unwrap();

        let nullity = null_basis.len();
        let rank = col_basis.len();

        // For 3×4 matrix: rank + nullity = 4
        assert_eq!(rank + nullity, 4);
    }

    #[test]
    fn test_null_space_dimension_theorem_4x3() {
        let a = Matrix::<f64, 4, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [0.0, 1.0, 2.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0],
        ]);

        let null_basis = null_space(&a).unwrap();
        let col_basis = column_space(&a).unwrap();

        let nullity = null_basis.len();
        let rank = col_basis.len();

        // For 4×3 matrix: rank + nullity = 3
        assert_eq!(rank + nullity, 3);
    }
}
