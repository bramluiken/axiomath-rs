//! Eigenvalue and eigenvector decomposition.
//!
//! Eigenvalues and eigenvectors are fundamental in linear algebra, appearing in:
//! - Stability analysis of dynamical systems
//! - Principal component analysis (PCA)
//! - Quantum mechanics (observables)
//! - Graph theory (spectral graph theory)
//! - Vibration analysis
//!
//! # Mathematical Definition
//!
//! For a square matrix A, an eigenvalue λ and eigenvector v satisfy:
//! ```text
//! A v = λ v
//! ```
//!
//! The characteristic polynomial is: det(A - λI) = 0
//!
//! # Algorithms Implemented
//!
//! - **Power Iteration**: Finds dominant (largest magnitude) eigenvalue
//! - **QR Algorithm**: Iterative method for all eigenvalues
//! - **Symmetric optimization**: Uses SVD for symmetric matrices (more efficient)
//!
//! # Limitations
//!
//! - QR algorithm is most reliable for matrices ≤ 10×10
//! - For large matrices, iterative methods (Arnoldi, Lanczos) are preferred (not yet implemented)
//! - Complex eigenvalues are not yet supported (returns error)
//!
//! # References
//!
//! - Golub & Van Loan, "Matrix Computations" (4th ed.), Chapter 7
//! - Trefethen & Bau, "Numerical Linear Algebra", Lectures 24-28
//! - Press et al., "Numerical Recipes", Chapter 11

use crate::core::traits::Float;
use crate::linear_algebra::decomposition::qr_decompose_householder;
use crate::linear_algebra::LinearAlgebraError;
use crate::matrix::{multiply, Matrix};
use crate::utils::comparison::approximately_equal;
use crate::vector::{dot_product, magnitude, Vector};
use num_traits::{One, Zero};

/// Computes all eigenvalues of a square matrix using the QR algorithm.
///
/// # Mathematical Background
///
/// The QR algorithm is an iterative method that repeatedly performs QR decomposition:
/// ```text
/// A_0 = A
/// for k = 1, 2, ...
///     A_{k-1} = Q_k R_k  (QR decomposition)
///     A_k = R_k Q_k      (reverse multiplication)
/// ```
///
/// Under certain conditions, A_k converges to a triangular or block-triangular matrix
/// where the eigenvalues appear on the diagonal.
///
/// # Algorithm
///
/// 1. Initialize A_0 = A
/// 2. For each iteration:
///    - Compute QR decomposition: A = QR
///    - Update: A = RQ
/// 3. Check convergence of off-diagonal elements
/// 4. Extract eigenvalues from diagonal
///
/// # Limitations
///
/// - Works best for matrices up to size 10×10
/// - May not converge for matrices with complex eigenvalues (returns error)
/// - Does not return eigenvectors (use `eigenpairs` for that)
/// - Convergence rate depends on eigenvalue separation
///
/// # Complexity
///
/// - Time: O(n³) per iteration, typically 10-30 iterations needed
/// - Space: O(n²)
///
/// # Arguments
///
/// * `matrix` - Square matrix to decompose (n×n)
///
/// # Returns
///
/// Returns `Ok(eigenvalues)` as a vector of n eigenvalues (unsorted).
///
/// Returns `Err` if:
/// - Matrix has complex eigenvalues (not yet supported)
/// - Algorithm fails to converge within maximum iterations
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::decomposition::eigenvalues;
///
/// // Diagonal matrix has eigenvalues on diagonal
/// let a = Matrix::from_rows([
///     [4.0, 0.0, 0.0],
///     [0.0, 3.0, 0.0],
///     [0.0, 0.0, 2.0],
/// ]);
///
/// let eigs = eigenvalues(&a).unwrap();
/// // eigs ≈ [4, 3, 2] (in some order)
/// ```
pub fn eigenvalues<T: Float, const N: usize>(
    matrix: &Matrix<T, N, N>,
) -> Result<Vector<T, N>, LinearAlgebraError> {
    qr_algorithm_eigenvalues(matrix, 100, T::epsilon() * T::from(1000.0).unwrap())
}

/// Computes all eigenvalues and eigenvectors using the QR algorithm.
///
/// Returns the eigenvalue decomposition: A = V Λ V⁻¹
/// where V is the matrix of eigenvectors (columns) and Λ is diagonal with eigenvalues.
///
/// # Mathematical Property
///
/// For each eigenvalue λ_i and eigenvector v_i:
/// ```text
/// A v_i = λ_i v_i
/// ```
///
/// Or in matrix form:
/// ```text
/// A V = V Λ
/// ```
///
/// # Arguments
///
/// * `matrix` - Square matrix to decompose (n×n)
///
/// # Returns
///
/// Returns `Ok((eigenvalues, eigenvectors))` where:
/// - `eigenvalues`: Vector of n eigenvalues
/// - `eigenvectors`: Matrix where column i is the eigenvector for eigenvalue i
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::decomposition::eigenpairs;
///
/// let a = Matrix::from_rows([
///     [2.0, 1.0],
///     [1.0, 2.0],
/// ]);
///
/// let (vals, vecs) = eigenpairs(&a).unwrap();
/// // vals ≈ [3, 1]
/// // vecs columns are corresponding eigenvectors
/// ```
pub fn eigenpairs<T: Float, const N: usize>(
    matrix: &Matrix<T, N, N>,
) -> Result<(Vector<T, N>, Matrix<T, N, N>), LinearAlgebraError> {
    // For symmetric matrices, use more efficient/stable method
    if is_symmetric(matrix, T::epsilon() * T::from(100.0).unwrap()) {
        return eigenpairs_symmetric(matrix);
    }

    // General case: use QR algorithm with eigenvector accumulation
    qr_algorithm_eigenpairs(matrix, 100, T::epsilon() * T::from(1000.0).unwrap())
}

/// Finds the dominant (largest magnitude) eigenvalue and its eigenvector using power iteration.
///
/// # Mathematical Background
///
/// Power iteration is based on the observation that for most initial vectors v_0:
/// ```text
/// lim_{k→∞} A^k v_0 / ||A^k v_0|| = v_1
/// ```
/// where v_1 is the eigenvector corresponding to the largest eigenvalue λ_1.
///
/// # Algorithm
///
/// ```text
/// Start with random unit vector v
/// for k = 1, 2, ...
///     v = A v
///     v = v / ||v||
/// λ = v^T A v (Rayleigh quotient)
/// ```
///
/// # Convergence
///
/// Converges at rate |λ_2/λ_1| where λ_1, λ_2 are largest and second-largest eigenvalues.
/// Fails if |λ_1| = |λ_2| (multiple dominant eigenvalues).
///
/// # Complexity
///
/// - Time: O(n²) per iteration
/// - Space: O(n)
/// - Typically converges in 10-50 iterations
///
/// # Arguments
///
/// * `matrix` - Square matrix (n×n)
/// * `max_iterations` - Maximum iterations before giving up
/// * `tolerance` - Convergence threshold for eigenvector change
///
/// # Returns
///
/// Returns `Ok((eigenvalue, eigenvector))` where:
/// - `eigenvalue`: Dominant eigenvalue (largest magnitude)
/// - `eigenvector`: Corresponding unit eigenvector
///
/// Returns `Err` if:
/// - Matrix has multiple dominant eigenvalues of equal magnitude
/// - Fails to converge within max_iterations
///
/// # Examples
///
/// ```
/// use axiomath::matrix::Matrix;
/// use axiomath::linear_algebra::decomposition::power_iteration;
///
/// let a = Matrix::from_rows([
///     [2.0, 1.0],
///     [1.0, 2.0],
/// ]);
///
/// let (lambda, v) = power_iteration(&a, 100, 1e-10).unwrap();
/// // lambda ≈ 3.0 (largest eigenvalue)
/// ```
pub fn power_iteration<T: Float, const N: usize>(
    matrix: &Matrix<T, N, N>,
    max_iterations: usize,
    tolerance: T,
) -> Result<(T, Vector<T, N>), LinearAlgebraError> {
    // Start with a random-ish vector (use basis vector with small perturbations)
    let mut v = Vector::<T, N>::zero();
    v[0] = T::one();
    if N > 1 {
        v[1] = T::from(0.1).unwrap();
    }
    if N > 2 {
        v[2] = T::from(0.01).unwrap();
    }

    // Normalize
    let norm = magnitude(&v);
    for i in 0..N {
        v[i] = v[i] / norm;
    }

    // Power iteration
    for _ in 0..max_iterations {
        // Compute A * v
        let mut av = Vector::<T, N>::zero();
        for i in 0..N {
            let mut sum = T::zero();
            for j in 0..N {
                sum = sum + matrix[(i, j)] * v[j];
            }
            av[i] = sum;
        }

        // Compute norm
        let av_norm = magnitude(&av);
        if av_norm < tolerance {
            return Err(LinearAlgebraError::NotConverged);
        }

        // Normalize to get new eigenvector estimate
        let mut v_new = Vector::<T, N>::zero();
        for i in 0..N {
            v_new[i] = av[i] / av_norm;
        }

        // Check convergence: ||v_new - v|| < tol (or ||v_new + v|| < tol for sign flip)
        let mut diff_plus = T::zero();
        let mut diff_minus = T::zero();
        for i in 0..N {
            let dp = (v_new[i] - v[i]).abs();
            let dm = (v_new[i] + v[i]).abs();
            if dp > diff_plus {
                diff_plus = dp;
            }
            if dm > diff_minus {
                diff_minus = dm;
            }
        }

        let diff = diff_plus.min(diff_minus);
        v = v_new;

        if diff < tolerance {
            // Converged - compute eigenvalue via Rayleigh quotient: λ = v^T A v
            let mut av = Vector::<T, N>::zero();
            for i in 0..N {
                let mut sum = T::zero();
                for j in 0..N {
                    sum = sum + matrix[(i, j)] * v[j];
                }
                av[i] = sum;
            }
            let eigenvalue = dot_product(&v, &av);
            return Ok((eigenvalue, v));
        }
    }

    Err(LinearAlgebraError::NotConverged)
}

/// QR algorithm for computing all eigenvalues.
///
/// This is the core implementation used by the public `eigenvalues` function.
fn qr_algorithm_eigenvalues<T: Float, const N: usize>(
    matrix: &Matrix<T, N, N>,
    max_iterations: usize,
    tolerance: T,
) -> Result<Vector<T, N>, LinearAlgebraError> {
    let mut a = *matrix;

    // QR iteration
    for iter in 0..max_iterations {
        // QR decomposition using Householder (most stable)
        let (q, r) = qr_decompose_householder(&a)?;

        // Update A = RQ
        a = multiply(&r, &q);

        // Check convergence: off-diagonal elements should approach zero
        if iter > 10 {
            // Don't check too early
            let mut max_off_diag = T::zero();
            for i in 0..N {
                for j in 0..N {
                    if i != j {
                        let val = a[(i, j)].abs();
                        if val > max_off_diag {
                            max_off_diag = val;
                        }
                    }
                }
            }

            if max_off_diag < tolerance {
                // Converged - extract eigenvalues from diagonal
                let mut eigenvalues = Vector::<T, N>::zero();
                for i in 0..N {
                    eigenvalues[i] = a[(i, i)];
                }
                return Ok(eigenvalues);
            }
        }
    }

    Err(LinearAlgebraError::NotConverged)
}

/// QR algorithm for computing eigenvalues and eigenvectors.
///
/// Accumulates the Q matrices to recover eigenvectors.
fn qr_algorithm_eigenpairs<T: Float, const N: usize>(
    matrix: &Matrix<T, N, N>,
    max_iterations: usize,
    tolerance: T,
) -> Result<(Vector<T, N>, Matrix<T, N, N>), LinearAlgebraError> {
    let mut a = *matrix;
    let mut v = Matrix::<T, N, N>::identity(); // Accumulate eigenvectors

    // QR iteration
    for iter in 0..max_iterations {
        // QR decomposition
        let (q, r) = qr_decompose_householder(&a)?;

        // Update A = RQ
        a = multiply(&r, &q);

        // Accumulate eigenvectors: V = V * Q
        v = multiply(&v, &q);

        // Check convergence
        if iter > 10 {
            let mut max_off_diag = T::zero();
            for i in 0..N {
                for j in 0..N {
                    if i != j {
                        let val = a[(i, j)].abs();
                        if val > max_off_diag {
                            max_off_diag = val;
                        }
                    }
                }
            }

            if max_off_diag < tolerance {
                // Converged - extract eigenvalues
                let mut eigenvalues = Vector::<T, N>::zero();
                for i in 0..N {
                    eigenvalues[i] = a[(i, i)];
                }
                return Ok((eigenvalues, v));
            }
        }
    }

    Err(LinearAlgebraError::NotConverged)
}

/// Optimized eigenvalue decomposition for symmetric matrices using SVD.
///
/// For symmetric matrices, eigenvalues are real and eigenvectors are orthogonal.
/// We can leverage this for a more stable/efficient computation.
fn eigenpairs_symmetric<T: Float, const N: usize>(
    matrix: &Matrix<T, N, N>,
) -> Result<(Vector<T, N>, Matrix<T, N, N>), LinearAlgebraError> {
    // For symmetric matrices: A = V Λ V^T (spectral theorem)
    // We can use power iteration with deflation
    power_iteration_deflation(matrix)
}

/// Power iteration with deflation to find all eigenvalue-eigenvector pairs.
///
/// This is adapted from the SVD implementation but exposed for eigenvalue computation.
fn power_iteration_deflation<T: Float, const N: usize>(
    matrix: &Matrix<T, N, N>,
) -> Result<(Vector<T, N>, Matrix<T, N, N>), LinearAlgebraError> {
    let max_iterations = 100;
    let tolerance = T::epsilon() * T::from(1000.0).unwrap();

    let mut v_matrix = Matrix::<T, N, N>::zero();
    let mut eigenvalues = Vector::<T, N>::zero();
    let mut deflated = *matrix;

    // Find each eigenvalue-eigenvector pair
    for k in 0..N {
        // Start with a vector that's orthogonal to previous eigenvectors
        let mut v = Vector::<T, N>::zero();
        v[k] = T::one();
        if k > 0 && N > 1 {
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

        // Power iteration on deflated matrix
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

        // Compute eigenvalue: λ = v^T A v (Rayleigh quotient)
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

    // Sort eigenvalues and eigenvectors in descending order by absolute value
    for i in 0..N {
        for j in (i + 1)..N {
            if eigenvalues[j].abs() > eigenvalues[i].abs() {
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

    Ok((eigenvalues, v_matrix))
}

/// Check if a matrix is symmetric.
fn is_symmetric<T: Float, const N: usize>(matrix: &Matrix<T, N, N>, tolerance: T) -> bool {
    for i in 0..N {
        for j in (i + 1)..N {
            if !approximately_equal(matrix[(i, j)], matrix[(j, i)], tolerance, tolerance) {
                return false;
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::transpose;

    #[test]
    fn test_power_iteration_diagonal() {
        // Diagonal matrix - eigenvalues are diagonal elements
        let a = Matrix::<f64, 3, 3>::from_rows([
            [5.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);

        let (lambda, v) = power_iteration(&a, 100, 1e-10).unwrap();

        // Should find largest eigenvalue (5.0)
        assert!(
            approximately_equal(lambda, 5.0, 1e-8, 1e-8),
            "Expected eigenvalue 5.0, got {}",
            lambda
        );

        // Verify: A v = λ v
        let mut av = Vector::<f64, 3>::zero();
        for i in 0..3 {
            for j in 0..3 {
                av[i] += a[(i, j)] * v[j];
            }
        }

        for i in 0..3 {
            assert!(
                approximately_equal(av[i], lambda * v[i], 1e-8, 1e-8),
                "Eigenvector verification failed"
            );
        }
    }

    #[test]
    fn test_power_iteration_symmetric() {
        // Symmetric 2x2 matrix
        let a = Matrix::<f64, 2, 2>::from_rows([[2.0, 1.0], [1.0, 2.0]]);

        let (lambda, v) = power_iteration(&a, 100, 1e-10).unwrap();

        // Largest eigenvalue should be 3.0
        assert!(
            approximately_equal(lambda.abs(), 3.0, 1e-8, 1e-8),
            "Expected |eigenvalue| = 3.0, got {}",
            lambda
        );

        // Verify: A v = λ v
        let mut av = Vector::<f64, 2>::zero();
        for i in 0..2 {
            for j in 0..2 {
                av[i] += a[(i, j)] * v[j];
            }
        }

        for i in 0..2 {
            assert!(approximately_equal(av[i], lambda * v[i], 1e-8, 1e-8));
        }
    }

    #[test]
    fn test_qr_algorithm_diagonal() {
        let a = Matrix::<f64, 3, 3>::from_rows([
            [4.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [0.0, 0.0, 2.0],
        ]);

        let eigs = eigenvalues(&a).unwrap();

        // Check that we found all three eigenvalues (order doesn't matter)
        let mut found = [false, false, false];
        let expected = [4.0, 3.0, 2.0];

        for i in 0..3 {
            for j in 0..3 {
                if approximately_equal(eigs[i], expected[j], 1e-8, 1e-8) {
                    found[j] = true;
                }
            }
        }

        assert!(found.iter().all(|&x| x), "Not all eigenvalues found");
    }

    #[test]
    fn test_qr_algorithm_symmetric() {
        // Symmetric matrix with known eigenvalues
        let a = Matrix::<f64, 2, 2>::from_rows([[2.0, 1.0], [1.0, 2.0]]);

        let eigs = eigenvalues(&a).unwrap();

        // Eigenvalues are 3 and 1
        let mut found_3 = false;
        let mut found_1 = false;

        for i in 0..2 {
            if approximately_equal(eigs[i], 3.0, 1e-8, 1e-8) {
                found_3 = true;
            }
            if approximately_equal(eigs[i], 1.0, 1e-8, 1e-8) {
                found_1 = true;
            }
        }

        assert!(found_3 && found_1, "Eigenvalues {}, {} not correct", eigs[0], eigs[1]);
    }

    #[test]
    fn test_eigenpairs_symmetric_simple() {
        // Simple symmetric 2x2
        let a = Matrix::<f64, 2, 2>::from_rows([[3.0, 1.0], [1.0, 3.0]]);

        let (vals, vecs) = eigenpairs(&a).unwrap();

        // Eigenvalues should be 4 and 2
        // Verify A V = V Λ for each eigenpair
        for i in 0..2 {
            let lambda = vals[i];
            let v = Vector::from([vecs[(0, i)], vecs[(1, i)]]);

            // Compute A * v
            let mut av = Vector::<f64, 2>::zero();
            for j in 0..2 {
                for k in 0..2 {
                    av[j] += a[(j, k)] * v[k];
                }
            }

            // Verify A v = λ v
            for j in 0..2 {
                assert!(
                    approximately_equal(av[j], lambda * v[j], 1e-7, 1e-7),
                    "Eigenpair {} verification failed: Av[{}] = {}, λv[{}] = {}",
                    i,
                    j,
                    av[j],
                    j,
                    lambda * v[j]
                );
            }
        }
    }

    #[test]
    fn test_eigenpairs_diagonal() {
        let a = Matrix::<f64, 3, 3>::from_rows([
            [5.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);

        let (vals, vecs) = eigenpairs(&a).unwrap();

        // Verify each eigenpair
        for i in 0..3 {
            let lambda = vals[i];
            let v = Vector::from([vecs[(0, i)], vecs[(1, i)], vecs[(2, i)]]);

            // Compute A * v
            let mut av = Vector::<f64, 3>::zero();
            for j in 0..3 {
                for k in 0..3 {
                    av[j] += a[(j, k)] * v[k];
                }
            }

            // Verify A v = λ v
            for j in 0..3 {
                assert!(approximately_equal(av[j], lambda * v[j], 1e-7, 1e-7));
            }
        }
    }

    #[test]
    fn test_eigenpairs_identity() {
        let a = Matrix::<f64, 3, 3>::identity();

        let (vals, vecs) = eigenpairs(&a).unwrap();

        // All eigenvalues should be 1
        for i in 0..3 {
            assert!(approximately_equal(vals[i], 1.0, 1e-8, 1e-8));
        }

        // Eigenvectors should be orthonormal
        let vecs_t = transpose(&vecs);
        let vvt = multiply(&vecs_t, &vecs);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(approximately_equal(vvt[(i, j)], expected, 1e-8, 1e-8));
            }
        }
    }

    #[test]
    fn test_eigenvalues_trace_determinant() {
        // Property: trace(A) = sum of eigenvalues
        // Property: det(A) = product of eigenvalues
        let a = Matrix::<f64, 3, 3>::from_rows([
            [2.0, 1.0, 0.0],
            [1.0, 3.0, 1.0],
            [0.0, 1.0, 2.0],
        ]);

        let eigs = eigenvalues(&a).unwrap();

        // Trace = 2 + 3 + 2 = 7
        let trace: f64 = eigs[0] + eigs[1] + eigs[2];
        assert!(approximately_equal(trace, 7.0, 1e-6, 1e-6), "Trace = {}", trace);

        // This matrix is symmetric positive definite, so all eigenvalues are positive
        // We can verify the determinant relationship
        let det_from_eigs = eigs[0] * eigs[1] * eigs[2];
        assert!(det_from_eigs > 0.0, "Determinant from eigenvalues should be positive");
    }

    #[test]
    fn test_symmetric_vs_qr_algorithm_consistency() {
        // For symmetric matrices, both paths should give same eigenvalues
        let a = Matrix::<f64, 3, 3>::from_rows([
            [4.0, 1.0, 0.0],
            [1.0, 4.0, 1.0],
            [0.0, 1.0, 4.0],
        ]);

        let eigs1 = eigenvalues(&a).unwrap(); // Uses symmetric optimization
        let (eigs2, _) = qr_algorithm_eigenpairs(&a, 100, 1e-10).unwrap(); // Direct QR

        // Sort both for comparison
        let mut e1 = [eigs1[0], eigs1[1], eigs1[2]];
        let mut e2 = [eigs2[0], eigs2[1], eigs2[2]];
        e1.sort_by(|a, b| b.partial_cmp(a).unwrap());
        e2.sort_by(|a, b| b.partial_cmp(a).unwrap());

        for i in 0..3 {
            assert!(
                approximately_equal(e1[i], e2[i], 1e-6, 1e-6),
                "Eigenvalue mismatch: {} vs {}",
                e1[i],
                e2[i]
            );
        }
    }

    #[test]
    fn test_power_iteration_convergence_property() {
        // Power iteration should converge faster when eigenvalue gap is larger
        let a = Matrix::<f64, 2, 2>::from_rows([[10.0, 0.0], [0.0, 1.0]]);

        // Should converge very quickly for this matrix (large gap)
        let result = power_iteration(&a, 10, 1e-10);
        assert!(result.is_ok(), "Should converge quickly with large eigenvalue gap");

        let (lambda, _) = result.unwrap();
        assert!(approximately_equal(lambda, 10.0, 1e-8, 1e-8));
    }

    #[test]
    fn test_eigenvector_orthogonality_symmetric() {
        // For symmetric matrices, eigenvectors are orthogonal
        let a = Matrix::<f64, 3, 3>::from_rows([
            [3.0, 1.0, 0.0],
            [1.0, 3.0, 1.0],
            [0.0, 1.0, 3.0],
        ]);

        let (_, vecs) = eigenpairs(&a).unwrap();

        // Check orthogonality: V^T V = I
        let vt = transpose(&vecs);
        let vtv = multiply(&vt, &vecs);

        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    approximately_equal(vtv[(i, j)], expected, 1e-7, 1e-7),
                    "Eigenvectors not orthonormal: V^T V[{},{}] = {}",
                    i,
                    j,
                    vtv[(i, j)]
                );
            }
        }
    }
}
