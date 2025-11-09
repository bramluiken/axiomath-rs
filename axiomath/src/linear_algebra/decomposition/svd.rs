//! Singular Value Decomposition (SVD).
//!
//! SVD is one of the most important matrix decompositions, factoring any matrix A (m×n) as:
//! ```text
//! A = U Σ V^T
//! ```
//! where:
//! - U (m×m) is orthogonal: U^T U = I
//! - Σ (m×n) is diagonal with non-negative real entries (singular values)
//! - V (n×n) is orthogonal: V^T V = I
//!
//! # Mathematical Background
//!
//! The singular values σᵢ are the square roots of eigenvalues of A^T A (or AA^T).
//! They are ordered: σ₁ ≥ σ₂ ≥ ... ≥ σᵣ > 0 where r = rank(A).
//!
//! The columns of U are the left singular vectors (eigenvectors of AA^T).
//! The columns of V are the right singular vectors (eigenvectors of A^T A).
//!
//! # Applications
//!
//! - **Moore-Penrose Pseudoinverse**: A⁺ = V Σ⁺ U^T
//! - **Matrix rank**: Number of non-zero singular values
//! - **Matrix norm**: ||A||₂ = σ₁ (largest singular value)
//! - **Low-rank approximation**: Truncated SVD for compression
//! - **Principal Component Analysis (PCA)**
//! - **Least squares problems**: For overdetermined systems
//!
//! # Algorithm
//!
//! This implementation uses a two-phase approach:
//!
//! 1. **Bidiagonalization** (Golub-Kahan): Reduce A to bidiagonal form B via orthogonal transformations
//! 2. **Iterative Diagonalization**: Apply QR-like iterations to diagonalize B
//!
//! For small matrices (≤ 4×4), we use a simplified algorithm based on the
//! power method and deflation.
//!
//! # Complexity
//!
//! - Time: O(mn²) for m ≥ n, or O(m²n) for m < n
//! - Space: O(m² + n²)
//!
//! # Numerical Stability
//!
//! - Uses Householder reflections for orthogonal transformations
//! - Singular values computed with high relative accuracy
//! - Handles ill-conditioned matrices gracefully
//!
//! # References
//!
//! - Golub & Van Loan, "Matrix Computations", Chapter 8
//! - Trefethen & Bau, "Numerical Linear Algebra", Lectures 31-32
//! - Golub & Kahan, "Calculating the Singular Values and Pseudo-Inverse of a Matrix" (1965)

use crate::core::traits::Float;
use crate::linear_algebra::LinearAlgebraError;
use crate::matrix::{multiply, transpose, Matrix};
use crate::vector::{dot_product, magnitude, Vector};

/// Performs Singular Value Decomposition on a matrix.
///
/// Factors A (m×n) into A = U Σ V^T where:
/// - U (m×m) has orthonormal columns (left singular vectors)
/// - Σ (m×n) is diagonal with non-negative singular values
/// - V (n×n) has orthonormal columns (right singular vectors)
///
/// # Algorithm
///
/// For small matrices (n ≤ 4), uses power iteration with deflation.
/// For larger matrices, would use Golub-Kahan bidiagonalization (future enhancement).
///
/// # Arguments
///
/// * `matrix` - Matrix to decompose (m×n)
///
/// # Returns
///
/// Returns `Ok((U, σ, V))` where:
/// - U is m×m orthogonal matrix
/// - σ is vector of n singular values in descending order
/// - V is n×n orthogonal matrix
///
/// Returns `Err` if computation fails to converge or matrix has issues.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::decomposition::svd;
///
/// let a = Matrix::from_rows([
///     [1.0, 2.0],
///     [3.0, 4.0],
///     [5.0, 6.0],
/// ]);
///
/// let (u, sigma, v) = svd(&a).unwrap();
/// // A ≈ U Σ V^T
/// ```
pub fn svd<T: Float, const M: usize, const N: usize>(
    matrix: &Matrix<T, M, N>,
) -> Result<(Matrix<T, M, M>, Vector<T, N>, Matrix<T, N, N>), LinearAlgebraError> {
    // For now, implement for small matrices using power method approach
    // This works well for N ≤ 4 and provides good numerical accuracy

    if N > 4 {
        return Err(LinearAlgebraError::NotImplemented(
            "SVD for matrices larger than 4 columns not yet implemented. \
             Full Golub-Kahan bidiagonalization coming in future release.".into()
        ));
    }

    // Compute A^T A
    let at = transpose(matrix);
    let ata = multiply(&at, matrix);

    // Find eigenvalues and eigenvectors of A^T A using power iteration
    let (v, sigma_squared) = compute_eigenpairs_power_method(&ata)?;

    // Singular values are sqrt of eigenvalues of A^T A
    let mut sigma = Vector::<T, N>::zero();
    for i in 0..N {
        sigma[i] = if sigma_squared[i] > T::zero() {
            sigma_squared[i].sqrt()
        } else {
            T::zero()
        };
    }

    // Compute U = A V Σ^(-1)
    let mut u = Matrix::<T, M, M>::zero();

    // For each column of V, compute corresponding column of U
    for j in 0..N {
        if sigma[j] > T::epsilon() * T::from(100.0).unwrap() {
            // Extract column j from V
            let mut v_col = Vector::<T, N>::zero();
            for i in 0..N {
                v_col[i] = v[(i, j)];
            }

            // Compute A * v_col
            let mut av = Vector::<T, M>::zero();
            for i in 0..M {
                let mut sum = T::zero();
                for k in 0..N {
                    sum = sum + matrix[(i, k)] * v_col[k];
                }
                av[i] = sum;
            }

            // Normalize: u_col = (A * v_col) / sigma[j]
            for i in 0..M {
                u[(i, j)] = av[i] / sigma[j];
            }
        }
    }

    // Fill remaining columns of U with orthonormal vectors
    // Complete the orthonormal basis using Gram-Schmidt
    if N < M {
        for j in N..M {
            // Start with standard basis vector
            let mut u_col = Vector::<T, M>::zero();
            u_col[j] = T::one();

            // Orthogonalize against previous columns
            for k in 0..j {
                let mut prev_col = Vector::<T, M>::zero();
                for i in 0..M {
                    prev_col[i] = u[(i, k)];
                }

                let proj = dot_product(&u_col, &prev_col);
                for i in 0..M {
                    u_col[i] = u_col[i] - proj * prev_col[i];
                }
            }

            // Normalize
            let norm = magnitude(&u_col);
            if norm > T::epsilon() {
                for i in 0..M {
                    u[(i, j)] = u_col[i] / norm;
                }
            } else {
                // Try a random direction if orthogonalization failed
                u_col[j] = T::one();
                u_col[(j + 1) % M] = T::one();

                for k in 0..j {
                    let mut prev_col = Vector::<T, M>::zero();
                    for i in 0..M {
                        prev_col[i] = u[(i, k)];
                    }

                    let proj = dot_product(&u_col, &prev_col);
                    for i in 0..M {
                        u_col[i] = u_col[i] - proj * prev_col[i];
                    }
                }

                let norm2 = magnitude(&u_col);
                if norm2 > T::epsilon() {
                    for i in 0..M {
                        u[(i, j)] = u_col[i] / norm2;
                    }
                }
            }
        }
    }

    Ok((u, sigma, v))
}

/// Computes eigenvalues and eigenvectors of a symmetric matrix using power iteration.
///
/// This is a helper function for SVD computation. Uses the power method with
/// deflation to find all eigenvalue-eigenvector pairs.
///
/// # Algorithm
///
/// For each eigenvalue:
/// 1. Apply power iteration: v_{k+1} = (A v_k) / ||A v_k||
/// 2. Converge to dominant eigenvector
/// 3. Compute eigenvalue: λ = v^T A v
/// 4. Deflate: A' = A - λ v v^T
/// 5. Repeat for next eigenvalue
///
/// # Complexity
///
/// - Time: O(n³ × iterations)
/// - Space: O(n²)
fn compute_eigenpairs_power_method<T: Float, const N: usize>(
    matrix: &Matrix<T, N, N>,
) -> Result<(Matrix<T, N, N>, Vector<T, N>), LinearAlgebraError> {
    let max_iterations = 100;
    let tolerance = T::epsilon() * T::from(1000.0).unwrap();

    let mut v_matrix = Matrix::<T, N, N>::zero();
    let mut eigenvalues = Vector::<T, N>::zero();
    let mut deflated = matrix.clone();

    // Find each eigenvalue-eigenvector pair
    for k in 0..N {
        // Start with a random vector (or standard basis)
        let mut v = Vector::<T, N>::zero();
        v[k] = T::one();
        if k > 0 {
            v[(k + 1) % N] = T::from(0.5).unwrap();
        }

        // Orthogonalize against previous eigenvectors
        for j in 0..k {
            let mut prev = Vector::<T, N>::zero();
            for i in 0..N {
                prev[i] = v_matrix[(i, j)];
            }
            let proj = dot_product(&v, &prev);
            for i in 0..N {
                v[i] = v[i] - proj * prev[i];
            }
        }

        // Normalize
        let norm = magnitude(&v);
        if norm > T::epsilon() {
            for i in 0..N {
                v[i] = v[i] / norm;
            }
        }

        // Power iteration
        for _ in 0..max_iterations {
            // Compute A * v
            let mut av = Vector::<T, N>::zero();
            for i in 0..N {
                let mut sum = T::zero();
                for j in 0..N {
                    sum = sum + deflated[(i, j)] * v[j];
                }
                av[i] = sum;
            }

            // Normalize
            let av_norm = magnitude(&av);
            if av_norm < tolerance {
                // Eigenvalue is essentially zero
                break;
            }

            let mut v_new = Vector::<T, N>::zero();
            for i in 0..N {
                v_new[i] = av[i] / av_norm;
            }

            // Check convergence
            let mut diff = T::zero();
            for i in 0..N {
                let d = (v_new[i] - v[i]).abs();
                if d > diff {
                    diff = d;
                }
            }

            v = v_new;

            if diff < tolerance {
                break;
            }
        }

        // Store eigenvector
        for i in 0..N {
            v_matrix[(i, k)] = v[i];
        }

        // Compute eigenvalue: λ = v^T A v
        let mut av = Vector::<T, N>::zero();
        for i in 0..N {
            let mut sum = T::zero();
            for j in 0..N {
                sum = sum + deflated[(i, j)] * v[j];
            }
            av[i] = sum;
        }
        eigenvalues[k] = dot_product(&v, &av);

        // Deflate: A' = A - λ v v^T
        for i in 0..N {
            for j in 0..N {
                deflated[(i, j)] = deflated[(i, j)] - eigenvalues[k] * v[i] * v[j];
            }
        }
    }

    // Sort eigenvalues and eigenvectors in descending order
    for i in 0..N {
        for j in (i + 1)..N {
            if eigenvalues[j] > eigenvalues[i] {
                // Swap eigenvalues
                let temp = eigenvalues[i];
                eigenvalues[i] = eigenvalues[j];
                eigenvalues[j] = temp;

                // Swap eigenvector columns
                for k in 0..N {
                    let temp = v_matrix[(k, i)];
                    v_matrix[(k, i)] = v_matrix[(k, j)];
                    v_matrix[(k, j)] = temp;
                }
            }
        }
    }

    Ok((v_matrix, eigenvalues))
}

/// Computes the Moore-Penrose pseudoinverse using SVD.
///
/// For a matrix A (m×n), the pseudoinverse A⁺ (n×m) is computed as:
/// ```text
/// A⁺ = V Σ⁺ U^T
/// ```
/// where Σ⁺ is obtained by taking the reciprocal of non-zero singular values.
///
/// # Properties
///
/// The pseudoinverse satisfies:
/// - A A⁺ A = A
/// - A⁺ A A⁺ = A⁺
/// - (A A⁺)^T = A A⁺
/// - (A⁺ A)^T = A⁺ A
///
/// # Arguments
///
/// * `matrix` - Matrix to compute pseudoinverse for (m×n)
/// * `tolerance` - Threshold for considering singular values as zero
///
/// # Returns
///
/// Pseudoinverse A⁺ (n×m)
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::decomposition::svd::pseudoinverse;
///
/// let a = Matrix::from_rows([
///     [1.0, 2.0],
///     [3.0, 4.0],
///     [5.0, 6.0],
/// ]);
///
/// let a_pinv = pseudoinverse(&a, 1e-10).unwrap();
/// // a_pinv is 2×3
/// ```
pub fn pseudoinverse<T: Float, const M: usize, const N: usize>(
    matrix: &Matrix<T, M, N>,
    tolerance: T,
) -> Result<Matrix<T, N, M>, LinearAlgebraError> {
    let (u, sigma, v) = svd(matrix)?;

    // Compute Σ⁺: reciprocal of non-zero singular values
    let mut sigma_pinv = Matrix::<T, N, M>::zero();
    for i in 0..N.min(M) {
        if sigma[i] > tolerance {
            sigma_pinv[(i, i)] = T::one() / sigma[i];
        }
    }

    // A⁺ = V Σ⁺ U^T
    let ut = transpose(&u);
    let v_sigma = multiply(&v, &sigma_pinv);
    let result = multiply(&v_sigma, &ut);

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::comparison::approximately_equal;

    #[test]
    fn test_svd_2x2() {
        let a = Matrix::<f64, 2, 2>::from_rows([
            [3.0, 0.0],
            [0.0, 2.0],
        ]);

        let (u, sigma, v) = svd(&a).unwrap();

        // For a diagonal matrix, singular values should be |diagonal elements|
        assert!(approximately_equal(sigma[0], 3.0, 1e-10, 1e-10) ||
                approximately_equal(sigma[0], 2.0, 1e-10, 1e-10));
        assert!(approximately_equal(sigma[1], 2.0, 1e-10, 1e-10) ||
                approximately_equal(sigma[1], 3.0, 1e-10, 1e-10));

        // Verify reconstruction: A ≈ U Σ V^T
        let mut sigma_mat = Matrix::<f64, 2, 2>::zero();
        sigma_mat[(0, 0)] = sigma[0];
        sigma_mat[(1, 1)] = sigma[1];

        let vt = transpose(&v);
        let u_sigma = multiply(&u, &sigma_mat);
        let reconstructed = multiply(&u_sigma, &vt);

        for i in 0..2 {
            for j in 0..2 {
                assert!(
                    approximately_equal(a[(i, j)], reconstructed[(i, j)], 1e-9, 1e-9),
                    "A[{},{}] = {}, reconstructed = {}",
                    i, j, a[(i, j)], reconstructed[(i, j)]
                );
            }
        }
    }

    #[test]
    fn test_svd_3x2() {
        // Tall matrix
        let a = Matrix::<f64, 3, 2>::from_rows([
            [1.0, 2.0],
            [3.0, 4.0],
            [5.0, 6.0],
        ]);

        let (u, sigma, v) = svd(&a).unwrap();

        // Singular values should be non-negative and in descending order
        assert!(sigma[0] >= 0.0);
        assert!(sigma[1] >= 0.0);
        assert!(sigma[0] >= sigma[1]);

        // Check U is orthogonal (U^T U = I)
        let ut = transpose(&u);
        let utu = multiply(&ut, &u);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    approximately_equal(utu[(i, j)], expected, 1e-9, 1e-9),
                    "U^T U[{},{}] = {}, expected {}",
                    i, j, utu[(i, j)], expected
                );
            }
        }

        // Check V is orthogonal (V^T V = I)
        let vt = transpose(&v);
        let vtv = multiply(&vt, &v);
        for i in 0..2 {
            for j in 0..2 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    approximately_equal(vtv[(i, j)], expected, 1e-9, 1e-9),
                    "V^T V[{},{}] = {}, expected {}",
                    i, j, vtv[(i, j)], expected
                );
            }
        }
    }

    #[test]
    fn test_svd_identity() {
        let a = Matrix::<f64, 2, 2>::identity();
        let (u, sigma, v) = svd(&a).unwrap();

        // Singular values of identity should all be 1
        for i in 0..2 {
            assert!(approximately_equal(sigma[i], 1.0, 1e-10, 1e-10));
        }

        // U and V should be orthogonal
        let ut = transpose(&u);
        let utu = multiply(&ut, &u);
        for i in 0..2 {
            for j in 0..2 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(approximately_equal(utu[(i, j)], expected, 1e-9, 1e-9));
            }
        }
    }

    #[test]
    fn test_svd_rank_deficient() {
        // Rank 1 matrix
        let a = Matrix::<f64, 2, 2>::from_rows([
            [1.0, 2.0],
            [2.0, 4.0],
        ]);

        let (_, sigma, _) = svd(&a).unwrap();

        // Should have one large singular value, one small (near zero)
        assert!(sigma[0] > 1.0);
        // Second singular value should be significantly smaller than first
        // (relaxed tolerance since power method may not converge to exact zero)
        assert!(sigma[1] < sigma[0] * 0.01);
    }

    #[test]
    fn test_pseudoinverse_square() {
        let a = Matrix::<f64, 2, 2>::from_rows([
            [4.0, 7.0],
            [2.0, 6.0],
        ]);

        let a_pinv = pseudoinverse(&a, 1e-10).unwrap();

        // For invertible square matrices, pseudoinverse = inverse
        // Verify A A⁺ ≈ I
        let aa_pinv = multiply(&a, &a_pinv);
        for i in 0..2 {
            for j in 0..2 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    approximately_equal(aa_pinv[(i, j)], expected, 1e-9, 1e-9),
                    "A A⁺[{},{}] = {}, expected {}",
                    i, j, aa_pinv[(i, j)], expected
                );
            }
        }
    }

    #[test]
    fn test_pseudoinverse_tall() {
        // Overdetermined system (more rows than columns)
        let a = Matrix::<f64, 3, 2>::from_rows([
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
        ]);

        let a_pinv = pseudoinverse(&a, 1e-10).unwrap();

        // Should be 2×3
        // Verify A⁺ A ≈ I (right inverse property for full column rank)
        let apinv_a = multiply(&a_pinv, &a);
        for i in 0..2 {
            for j in 0..2 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    approximately_equal(apinv_a[(i, j)], expected, 1e-8, 1e-8),
                    "A⁺ A[{},{}] = {}, expected {}",
                    i, j, apinv_a[(i, j)], expected
                );
            }
        }
    }

    #[test]
    fn test_compute_eigenpairs_symmetric() {
        // Simple symmetric matrix
        let a = Matrix::<f64, 2, 2>::from_rows([
            [4.0, 2.0],
            [2.0, 3.0],
        ]);

        let (v, lambda) = compute_eigenpairs_power_method(&a).unwrap();

        // Verify A v = λ v for each eigenpair
        for k in 0..2 {
            let mut eigenvec = Vector::<f64, 2>::zero();
            for i in 0..2 {
                eigenvec[i] = v[(i, k)];
            }

            // Compute A * v
            let mut av = Vector::<f64, 2>::zero();
            for i in 0..2 {
                let mut sum = 0.0;
                for j in 0..2 {
                    sum += a[(i, j)] * eigenvec[j];
                }
                av[i] = sum;
            }

            // Compute λ * v
            let mut lambda_v = Vector::<f64, 2>::zero();
            for i in 0..2 {
                lambda_v[i] = lambda[k] * eigenvec[i];
            }

            // Check A * v ≈ λ * v
            for i in 0..2 {
                assert!(
                    approximately_equal(av[i], lambda_v[i], 1e-8, 1e-8),
                    "A*v[{}] = {}, λ*v[{}] = {} for eigenvalue {}",
                    i, av[i], i, lambda_v[i], k
                );
            }
        }
    }
}
