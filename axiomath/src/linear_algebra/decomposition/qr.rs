//! QR decomposition algorithms.
//!
//! QR decomposition factors a matrix A as A = QR where:
//! - Q is an orthogonal matrix (Q^T Q = I)
//! - R is an upper triangular matrix
//!
//! # Algorithms
//!
//! - **Classical Gram-Schmidt**: Simple but can be numerically unstable
//! - **Modified Gram-Schmidt**: More numerically stable variant
//! - **Householder Reflections**: Most stable for ill-conditioned matrices
//!
//! # Applications
//!
//! - Solving least squares problems
//! - Computing eigenvalues (QR algorithm)
//! - Orthogonal basis construction
//!
//! # References
//!
//! - Trefethen & Bau, "Numerical Linear Algebra", Lectures 7-10
//! - Golub & Van Loan, "Matrix Computations", Section 5.2

use crate::core::traits::Float;
use crate::linear_algebra::LinearAlgebraError;
use crate::matrix::Matrix;
use crate::vector::{dot_product, magnitude, scale, subtract, Vector};
use num_traits::{One, Zero};

/// Performs QR decomposition using Modified Gram-Schmidt process.
///
/// Factors a matrix A (m×n) into A = QR where:
/// - Q is orthogonal (m×m): Q^T Q = I
/// - R is upper triangular (m×n)
///
/// This is the recommended QR decomposition method for most applications
/// due to its good balance of simplicity and numerical stability.
///
/// # Mathematical Definition
///
/// ```text
/// A = QR
/// ```
///
/// where:
/// - Q^T Q = I (orthogonality)
/// - R_ij = 0 for i > j (upper triangular)
///
/// # Algorithm: Modified Gram-Schmidt
///
/// ```text
/// for j = 1 to n:
///     q_j = a_j
///     for i = 1 to j-1:
///         r_ij = q_i^T * q_j
///         q_j = q_j - r_ij * q_i
///     r_jj = ||q_j||
///     q_j = q_j / r_jj
/// ```
///
/// Modified Gram-Schmidt performs orthogonalization in a different order
/// than Classical GS, which improves numerical stability.
///
/// # Numerical Stability
///
/// Modified Gram-Schmidt is more stable than Classical Gram-Schmidt,
/// but still can have issues with nearly-dependent columns. For maximum
/// stability, use Householder reflections (not yet implemented).
///
/// # Complexity
///
/// - Time: O(mn²) for m×n matrix
/// - Space: O(mn) for Q and R
///
/// # Arguments
///
/// * `matrix` - Matrix to decompose (m×n)
///
/// # Returns
///
/// Returns `Ok((Q, R))` where Q is m×m orthogonal and R is m×n upper triangular.
///
/// Returns `Err(LinearAlgebraError::SingularMatrix)` if the matrix has
/// linearly dependent columns.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::decomposition::qr_decompose;
///
/// let a = Matrix::from_rows([
///     [12.0, -51.0, 4.0],
///     [6.0, 167.0, -68.0],
///     [-4.0, 24.0, -41.0],
/// ]);
///
/// let (q, r) = qr_decompose(&a).unwrap();
///
/// // Verify: Q^T Q = I and A = QR
/// ```
pub fn qr_decompose<T: Float, const M: usize, const N: usize>(
    matrix: &Matrix<T, M, N>,
) -> Result<(Matrix<T, M, M>, Matrix<T, M, N>), LinearAlgebraError> {
    qr_decompose_gram_schmidt(matrix)
}

/// QR decomposition using Modified Gram-Schmidt process.
///
/// This is the implementation used by [`qr_decompose`]. See that function
/// for detailed documentation.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::decomposition::qr_decompose_gram_schmidt;
///
/// let a = Matrix::from_rows([
///     [1.0, 1.0],
///     [1.0, 0.0],
///     [1.0, -1.0],
/// ]);
///
/// let (q, r) = qr_decompose_gram_schmidt(&a).unwrap();
/// ```
pub fn qr_decompose_gram_schmidt<T: Float, const M: usize, const N: usize>(
    matrix: &Matrix<T, M, N>,
) -> Result<(Matrix<T, M, M>, Matrix<T, M, N>), LinearAlgebraError> {
    let epsilon = T::epsilon() * T::from(M.max(N)).unwrap() * T::from(100.0).unwrap();

    // Extract columns as vectors
    let mut q_vecs: Vec<Vector<T, M>> = Vec::new();
    let mut r = Matrix::<T, M, N>::zero();

    // Process each column
    for j in 0..N {
        // Extract column j
        let mut v = Vector::<T, M>::zero();
        for i in 0..M {
            v[i] = matrix[(i, j)];
        }

        // Modified Gram-Schmidt: orthogonalize against all previous q vectors
        for i in 0..j {
            let r_ij = dot_product(&q_vecs[i], &v);
            r[(i, j)] = r_ij;
            v = subtract(&v, &scale(&q_vecs[i], r_ij));
        }

        // Compute norm
        let norm = magnitude(&v);

        // Check for linear dependence
        if norm < epsilon {
            return Err(LinearAlgebraError::SingularMatrix);
        }

        r[(j, j)] = norm;

        // Normalize
        let q_j = scale(&v, T::one() / norm);
        q_vecs.push(q_j);
    }

    // If N < M, we need to extend Q to a full orthonormal basis
    // For now, we'll just fill with zeros and add orthonormal vectors
    while q_vecs.len() < M {
        // Find a vector orthogonal to all existing q vectors
        let mut v = Vector::<T, M>::zero();

        // Try standard basis vectors
        for k in 0..M {
            v = Vector::<T, M>::zero();
            v[k] = T::one();

            // Orthogonalize against existing q vectors
            for q in &q_vecs {
                let proj = dot_product(q, &v);
                v = subtract(&v, &scale(q, proj));
            }

            let norm = magnitude(&v);
            if norm > epsilon {
                v = scale(&v, T::one() / norm);
                q_vecs.push(v);
                break;
            }
        }
    }

    // Build Q matrix from q vectors
    let mut q = Matrix::<T, M, M>::zero();
    for j in 0..M {
        for i in 0..M {
            q[(i, j)] = q_vecs[j][i];
        }
    }

    Ok((q, r))
}

/// QR decomposition using Householder reflections.
///
/// This method is more numerically stable than Gram-Schmidt, especially
/// for ill-conditioned matrices. It's the recommended method for production
/// use when numerical stability is critical.
///
/// # Mathematical Definition
///
/// Uses Householder reflections H_k to zero out elements below the diagonal:
/// ```text
/// H_n ... H_2 H_1 A = R
/// Q = H_1 H_2 ... H_n
/// ```
///
/// Each H_k is of the form: H = I - 2vv^T where v is a unit vector (Householder vector).
///
/// # Algorithm
///
/// ```text
/// for k = 1 to min(m-1, n):
///     x = A[k:m, k]  (column k from row k downward)
///     v = sign(x_1) * ||x|| * e_1 + x
///     v = v / ||v||
///     H_k = I - 2vv^T
///     A = H_k * A
/// R = A
/// Q = H_1 * H_2 * ... * H_k
/// ```
///
/// The Householder vector v is chosen to reflect column k below the diagonal
/// to align with the first coordinate axis, thus creating zeros.
///
/// # Numerical Stability
///
/// Householder reflections are backward stable and maintain orthogonality
/// better than Gram-Schmidt methods. The conditioning of Q is independent
/// of the conditioning of A.
///
/// # Complexity
///
/// - Time: O(mn² - n³/3) for m×n matrix with m ≥ n
/// - Space: O(mn) for Q and R
///
/// # Arguments
///
/// * `matrix` - Matrix to decompose (m×n)
///
/// # Returns
///
/// Returns `Ok((Q, R))` where Q is m×m orthogonal and R is m×n upper triangular.
///
/// Returns `Err(LinearAlgebraError::SingularMatrix)` if the matrix has
/// linearly dependent columns (column norm is effectively zero).
///
/// # References
///
/// - Golub & Van Loan, "Matrix Computations" (4th ed.), Algorithm 5.1.1
/// - Trefethen & Bau, "Numerical Linear Algebra", Lecture 10
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::decomposition::qr_decompose_householder;
///
/// let a = Matrix::from_rows([
///     [12.0, -51.0, 4.0],
///     [6.0, 167.0, -68.0],
///     [-4.0, 24.0, -41.0],
/// ]);
///
/// let (q, r) = qr_decompose_householder(&a).unwrap();
/// // Q is orthogonal and R is upper triangular
/// ```
pub fn qr_decompose_householder<T: Float, const M: usize, const N: usize>(
    matrix: &Matrix<T, M, N>,
) -> Result<(Matrix<T, M, M>, Matrix<T, M, N>), LinearAlgebraError> {
    let epsilon = T::epsilon() * T::from(M.max(N)).unwrap() * T::from(100.0).unwrap();

    // Initialize R as a copy of A
    let mut r = *matrix;

    // Store Householder vectors for computing Q later
    let mut householder_vecs: Vec<(usize, Vector<T, M>)> = Vec::new();

    // Apply Householder reflections to each column
    let k_max = M.min(N);
    for k in 0..k_max {
        // Extract column k from row k downward: x = A[k:m, k]
        let mut x = Vector::<T, M>::zero();
        let len = M - k;
        for i in k..M {
            x[i - k] = r[(i, k)];
        }

        // Compute norm of x (only non-zero part)
        let mut norm_sq = T::zero();
        for i in 0..len {
            norm_sq = norm_sq + x[i] * x[i];
        }
        let norm_x = norm_sq.sqrt();

        // Check for near-zero column (singular matrix)
        if norm_x < epsilon {
            return Err(LinearAlgebraError::SingularMatrix);
        }

        // Compute Householder vector v = sign(x_1) * ||x|| * e_1 + x
        // Using sign(x_1) instead of -sign(x_1) for numerical stability
        let sign = if x[0] >= T::zero() { T::one() } else { -T::one() };
        let alpha = sign * norm_x;

        let mut v = x;
        v[0] = v[0] + alpha;

        // Normalize v to unit vector
        let mut v_norm_sq = T::zero();
        for i in 0..len {
            v_norm_sq = v_norm_sq + v[i] * v[i];
        }
        let v_norm = v_norm_sq.sqrt();

        if v_norm < epsilon {
            // This shouldn't happen if norm_x > epsilon, but handle defensively
            continue;
        }

        for i in 0..len {
            v[i] = v[i] / v_norm;
        }

        // Store the Householder vector (with its starting index k)
        let mut v_full = Vector::<T, M>::zero();
        for i in 0..len {
            v_full[i] = v[i];
        }
        householder_vecs.push((k, v_full));

        // Apply Householder reflection to R: R = (I - 2vv^T) R
        // This is equivalent to: R = R - 2v(v^T R)
        // We only need to update rows k through M and columns k through N

        for j in k..N {
            // Compute v^T * R[:, j] (only the part from row k onward)
            let mut vt_r_col = T::zero();
            for i in 0..len {
                vt_r_col = vt_r_col + v[i] * r[(k + i, j)];
            }

            // Update column j: R[:, j] = R[:, j] - 2 * v * (v^T * R[:, j])
            let factor = T::from(2.0).unwrap() * vt_r_col;
            for i in 0..len {
                r[(k + i, j)] = r[(k + i, j)] - factor * v[i];
            }
        }
    }

    // Now construct Q^T = H_k * ... * H_2 * H_1 (applied left to right)
    // Then transpose to get Q
    // Start with identity matrix for Q^T
    let mut qt = Matrix::<T, M, M>::identity();

    // Apply Householder reflections in reverse order (H_k * ... * H_1)
    for (k, v) in householder_vecs.iter().rev() {
        let k = *k;
        let len = M - k;

        // Apply H_k to qt: qt = H_k * qt = (I - 2vv^T) * qt
        // This is: qt = qt - 2v(v^T * qt)
        // We process each row of qt
        for i in 0..M {
            // Compute (v^T * qt[i, :]) for the portion affected (columns k onward)
            let mut vt_row = T::zero();
            for j in 0..len {
                vt_row = vt_row + v[j] * qt[(i, k + j)];
            }

            // Update row i: qt[i, :] = qt[i, :] - 2 * vt_row * v
            let factor = T::from(2.0).unwrap() * vt_row;
            for j in 0..len {
                qt[(i, k + j)] = qt[(i, k + j)] - factor * v[j];
            }
        }
    }

    // Transpose to get Q
    let q = crate::matrix::transpose(&qt);

    Ok((q, r))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::{multiply, transpose};
    use crate::utils::comparison::approximately_equal;

    #[test]
    fn test_qr_decompose_3x3() {
        // Example from Trefethen & Bau
        let a = Matrix::<f64, 3, 3>::from_rows([
            [12.0, -51.0, 4.0],
            [6.0, 167.0, -68.0],
            [-4.0, 24.0, -41.0],
        ]);

        let (q, r) = qr_decompose(&a).unwrap();

        // Check Q is orthogonal: Q^T Q = I
        let qt = transpose(&q);
        let qtq = multiply(&qt, &q);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    approximately_equal(qtq[(i, j)], expected, 1e-10, 1e-10),
                    "Q^T Q[{},{}] = {}, expected {}",
                    i,
                    j,
                    qtq[(i, j)],
                    expected
                );
            }
        }

        // Check R is upper triangular
        for i in 1..3 {
            for j in 0..i {
                assert!(
                    r[(i, j)].abs() < 1e-10_f64,
                    "R[{},{}] = {}, should be zero",
                    i,
                    j,
                    r[(i, j)]
                );
            }
        }

        // Check A = QR
        let qr = multiply(&q, &r);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    approximately_equal(a[(i, j)], qr[(i, j)], 1e-9, 1e-9),
                    "A[{},{}] = {}, QR[{},{}] = {}",
                    i,
                    j,
                    a[(i, j)],
                    i,
                    j,
                    qr[(i, j)]
                );
            }
        }
    }

    #[test]
    fn test_qr_decompose_3x2() {
        // Tall matrix (overdetermined system)
        let a = Matrix::<f64, 3, 2>::from_rows([[1.0, 1.0], [1.0, 0.0], [1.0, -1.0]]);

        let (q, r) = qr_decompose(&a).unwrap();

        // Check Q is orthogonal: Q^T Q = I
        let qt = transpose(&q);
        let qtq = multiply(&qt, &q);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    approximately_equal(qtq[(i, j)], expected, 1e-10, 1e-10),
                    "Q^T Q[{},{}] = {}",
                    i,
                    j,
                    qtq[(i, j)]
                );
            }
        }

        // Check R is upper triangular (first 2 columns)
        for j in 0..2 {
            for i in (j + 1)..3 {
                assert!(r[(i, j)].abs() < 1e-10_f64, "R[{},{}] = {}", i, j, r[(i, j)]);
            }
        }

        // Check A = QR (first 2 columns of Q times R)
        let qr = multiply(&q, &r);
        for i in 0..3 {
            for j in 0..2 {
                assert!(
                    approximately_equal(a[(i, j)], qr[(i, j)], 1e-9, 1e-9),
                    "A[{},{}] = {}, QR[{},{}] = {}",
                    i,
                    j,
                    a[(i, j)],
                    i,
                    j,
                    qr[(i, j)]
                );
            }
        }
    }

    #[test]
    fn test_qr_decompose_identity() {
        let a = Matrix::<f64, 3, 3>::identity();
        let (q, r) = qr_decompose(&a).unwrap();

        // Q should be identity (or a permutation with sign changes)
        // R should be identity (or diagonal with ±1)
        let qr = multiply(&q, &r);

        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(approximately_equal(qr[(i, j)], expected, 1e-10, 1e-10));
            }
        }
    }

    #[test]
    fn test_qr_decompose_singular() {
        // Matrix with linearly dependent columns
        let a = Matrix::from_rows([[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [1.0, 2.0, 3.0]]);

        let result = qr_decompose(&a);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), LinearAlgebraError::SingularMatrix);
    }

    #[test]
    fn test_qr_orthogonality_property() {
        let a = Matrix::<f64, 3, 2>::from_rows([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]);

        let (q, _) = qr_decompose(&a).unwrap();

        // Check that columns of Q are orthonormal
        for i in 0..3 {
            // Extract column i
            let col_i = Vector::from([q[(0, i)], q[(1, i)], q[(2, i)]]);

            for j in 0..3 {
                // Extract column j
                let col_j = Vector::from([q[(0, j)], q[(1, j)], q[(2, j)]]);

                let dot = dot_product(&col_i, &col_j);
                let expected = if i == j { 1.0 } else { 0.0 };

                assert!(
                    approximately_equal(dot, expected, 1e-10, 1e-10),
                    "q_{} · q_{} = {}, expected {}",
                    i,
                    j,
                    dot,
                    expected
                );
            }
        }
    }

    // ========== Householder QR Tests ==========

    #[test]
    fn test_householder_qr_basic_3x3() {
        // Same test matrix as Gram-Schmidt test
        let a = Matrix::<f64, 3, 3>::from_rows([
            [12.0, -51.0, 4.0],
            [6.0, 167.0, -68.0],
            [-4.0, 24.0, -41.0],
        ]);

        let (q, r) = qr_decompose_householder(&a).unwrap();

        // Check Q is orthogonal: Q^T Q = I
        let qt = transpose(&q);
        let qtq = multiply(&qt, &q);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    approximately_equal(qtq[(i, j)], expected, 1e-10, 1e-10),
                    "Q^T Q[{},{}] = {}, expected {}",
                    i,
                    j,
                    qtq[(i, j)],
                    expected
                );
            }
        }

        // Check R is upper triangular
        for i in 1..3 {
            for j in 0..i {
                assert!(
                    r[(i, j)].abs() < 1e-10_f64,
                    "R[{},{}] = {}, should be zero",
                    i,
                    j,
                    r[(i, j)]
                );
            }
        }

        // Check A = QR
        let qr = multiply(&q, &r);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    approximately_equal(a[(i, j)], qr[(i, j)], 1e-9, 1e-9),
                    "A[{},{}] = {}, QR[{},{}] = {}",
                    i,
                    j,
                    a[(i, j)],
                    i,
                    j,
                    qr[(i, j)]
                );
            }
        }
    }

    #[test]
    fn test_householder_qr_tall_matrix() {
        // Tall matrix (m > n)
        let a = Matrix::<f64, 4, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.1], // Slightly perturbed to avoid singularity
            [10.0, 11.0, 12.0],
        ]);

        let (q, r) = qr_decompose_householder(&a).unwrap();

        // Property 1: Q^T Q = I (orthogonality)
        let qt = transpose(&q);
        let qtq = multiply(&qt, &q);
        for i in 0..4 {
            for j in 0..4 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    approximately_equal(qtq[(i, j)], expected, 1e-9, 1e-9),
                    "Q^T Q[{},{}] = {}",
                    i,
                    j,
                    qtq[(i, j)]
                );
            }
        }

        // Property 2: R is upper triangular
        for i in 0..4 {
            for j in 0..3 {
                if i > j {
                    assert!(
                        r[(i, j)].abs() < 1e-9_f64,
                        "R[{},{}] = {}, should be zero",
                        i,
                        j,
                        r[(i, j)]
                    );
                }
            }
        }

        // Property 3: A = QR (reconstruction)
        let qr = multiply(&q, &r);
        for i in 0..4 {
            for j in 0..3 {
                assert!(
                    approximately_equal(a[(i, j)], qr[(i, j)], 1e-8, 1e-8),
                    "A[{},{}] = {}, QR[{},{}] = {}",
                    i,
                    j,
                    a[(i, j)],
                    i,
                    j,
                    qr[(i, j)]
                );
            }
        }
    }

    #[test]
    fn test_householder_qr_wide_matrix() {
        // Wide matrix (m < n) - should still work
        let a = Matrix::<f64, 2, 3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);

        let (q, r) = qr_decompose_householder(&a).unwrap();

        // Q should be 2x2 and orthogonal
        let qt = transpose(&q);
        let qtq = multiply(&qt, &q);
        for i in 0..2 {
            for j in 0..2 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    approximately_equal(qtq[(i, j)], expected, 1e-10, 1e-10),
                    "Q^T Q[{},{}] = {}",
                    i,
                    j,
                    qtq[(i, j)]
                );
            }
        }

        // R should be 2x3 upper triangular (only first 2 rows matter)
        for i in 1..2 {
            for j in 0..i {
                assert!(r[(i, j)].abs() < 1e-10_f64);
            }
        }

        // Reconstruction
        let qr = multiply(&q, &r);
        for i in 0..2 {
            for j in 0..3 {
                assert!(approximately_equal(a[(i, j)], qr[(i, j)], 1e-9, 1e-9));
            }
        }
    }

    #[test]
    fn test_householder_vs_gram_schmidt_consistency() {
        // Both methods should give similar results (up to sign differences)
        let a = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.1], // Slightly perturbed
            [7.0, 8.1, 9.0],
        ]);

        let (q_gs, r_gs) = qr_decompose_gram_schmidt(&a).unwrap();
        let (q_hh, r_hh) = qr_decompose_householder(&a).unwrap();

        // Both should reconstruct A
        let qr_gs = multiply(&q_gs, &r_gs);
        let qr_hh = multiply(&q_hh, &r_hh);

        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    approximately_equal(qr_gs[(i, j)], a[(i, j)], 1e-8, 1e-8),
                    "GS reconstruction failed"
                );
                assert!(
                    approximately_equal(qr_hh[(i, j)], a[(i, j)], 1e-8, 1e-8),
                    "Householder reconstruction failed"
                );
            }
        }

        // R matrices should have same absolute values on diagonal (up to sign)
        for i in 0..3 {
            assert!(
                approximately_equal(r_gs[(i, i)].abs(), r_hh[(i, i)].abs(), 1e-8, 1e-8),
                "Diagonal elements differ: GS={}, HH={}",
                r_gs[(i, i)],
                r_hh[(i, i)]
            );
        }
    }

    #[test]
    fn test_householder_qr_identity() {
        let a = Matrix::<f64, 3, 3>::identity();
        let (q, r) = qr_decompose_householder(&a).unwrap();

        // Q and R should both be close to identity (up to sign)
        let qr = multiply(&q, &r);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    approximately_equal(qr[(i, j)], expected, 1e-10, 1e-10),
                    "Identity reconstruction failed"
                );
            }
        }
    }

    #[test]
    fn test_householder_qr_diagonal() {
        // Diagonal matrix - should be trivial decomposition
        let a = Matrix::<f64, 3, 3>::from_rows([
            [3.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);

        let (q, r) = qr_decompose_householder(&a).unwrap();

        // Q^T Q = I
        let qt = transpose(&q);
        let qtq = multiply(&qt, &q);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(approximately_equal(qtq[(i, j)], expected, 1e-10, 1e-10));
            }
        }

        // R should be upper triangular (diagonal in this case)
        for i in 1..3 {
            for j in 0..i {
                assert!(r[(i, j)].abs() < 1e-10_f64);
            }
        }

        // Reconstruction
        let qr = multiply(&q, &r);
        for i in 0..3 {
            for j in 0..3 {
                assert!(approximately_equal(a[(i, j)], qr[(i, j)], 1e-9, 1e-9));
            }
        }
    }

    #[test]
    fn test_householder_qr_singular_matrix() {
        // Matrix with linearly dependent columns
        let a = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [2.0, 4.0, 6.0],
            [1.0, 2.0, 3.0],
        ]);

        let result = qr_decompose_householder(&a);
        assert!(
            result.is_err(),
            "Expected error for singular matrix, got Ok"
        );
        assert_eq!(result.unwrap_err(), LinearAlgebraError::SingularMatrix);
    }

    #[test]
    fn test_householder_numerical_stability() {
        // Ill-conditioned matrix where Householder should be more stable
        // Hilbert matrix is notoriously ill-conditioned
        let a = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 1.0 / 2.0, 1.0 / 3.0],
            [1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0],
            [1.0 / 3.0, 1.0 / 4.0, 1.0 / 5.0],
        ]);

        let (q, r) = qr_decompose_householder(&a).unwrap();

        // Property: Orthogonality should be preserved even for ill-conditioned matrices
        let qt = transpose(&q);
        let qtq = multiply(&qt, &q);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    approximately_equal(qtq[(i, j)], expected, 1e-8, 1e-8),
                    "Orthogonality lost for ill-conditioned matrix"
                );
            }
        }

        // Property: Reconstruction should still work
        let qr = multiply(&q, &r);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    approximately_equal(a[(i, j)], qr[(i, j)], 1e-7, 1e-7),
                    "Reconstruction failed for Hilbert matrix"
                );
            }
        }
    }
}
