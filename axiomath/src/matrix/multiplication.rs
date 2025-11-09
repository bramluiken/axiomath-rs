//! Matrix multiplication operations.
//!
//! This module implements standard matrix multiplication, which is fundamentally
//! different from element-wise (Hadamard) multiplication.
//!
//! # Mathematical Background
//!
//! For matrices A (M×N) and B (N×P), the matrix product AB is an M×P matrix where:
//!
//! ```text
//! (AB)_ij = Σ(k=1 to N) A_ik · B_kj
//! ```
//!
//! The key constraint is that the number of columns in A must equal the number
//! of rows in B. This is enforced at compile-time using const generics.
//!
//! # Type Safety
//!
//! ```compile_fail
//! use axiomath::matrix::{Matrix, multiply};
//!
//! let a = Matrix::<f64, 2, 3>::zero();  // 2×3
//! let b = Matrix::<f64, 2, 4>::zero();  // 2×4
//!
//! // This will NOT compile - dimensions incompatible!
//! let c = multiply(&a, &b);
//! ```
//!
//! # Properties
//!
//! - Associative: (AB)C = A(BC)
//! - Distributive: A(B + C) = AB + AC
//! - NOT commutative: AB ≠ BA in general
//! - Identity: AI = IA = A
//! - Transpose: (AB)^T = B^T A^T
//!
//! # Examples
//!
//! ```
//! use axiomath::matrix::{Matrix, multiply};
//!
//! let a = Matrix::<f64, 2, 3>::from_rows([
//!     [1.0, 2.0, 3.0],
//!     [4.0, 5.0, 6.0],
//! ]);
//!
//! let b = Matrix::<f64, 3, 2>::from_rows([
//!     [7.0, 8.0],
//!     [9.0, 10.0],
//!     [11.0, 12.0],
//! ]);
//!
//! let c = multiply(&a, &b);  // Results in 2×2 matrix
//! assert_eq!(c[(0, 0)], 58.0);  // 1*7 + 2*9 + 3*11
//! ```

use crate::core::Scalar;
use crate::matrix::Matrix;
use crate::vector::Vector;
use num_traits::Zero;

/// Multiplies two matrices using standard matrix multiplication.
///
/// # Type Signature
///
/// ```text
/// Matrix<T, M, N> × Matrix<T, N, P> → Matrix<T, M, P>
/// ```
///
/// # Mathematical Definition
///
/// For matrices A (M×N) and B (N×P), the product C = AB is defined as:
/// ```text
/// C_ij = Σ(k=1 to N) A_ik · B_kj
/// ```
///
/// # Complexity
///
/// - Time: O(M·N·P)
/// - Space: O(M·P)
///
/// # Algorithm
///
/// Uses the standard naive triply-nested loop algorithm. For large matrices,
/// more sophisticated algorithms (Strassen, etc.) may be more efficient,
/// but this implementation prioritizes clarity and correctness.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, multiply};
///
/// let a = Matrix::<f64, 2, 3>::from_rows([
///     [1.0, 2.0, 3.0],
///     [4.0, 5.0, 6.0],
/// ]);
///
/// let b = Matrix::<f64, 3, 2>::from_rows([
///     [7.0, 8.0],
///     [9.0, 10.0],
///     [11.0, 12.0],
/// ]);
///
/// let c = multiply(&a, &b);
///
/// assert_eq!(c.dimensions(), (2, 2));
/// assert_eq!(c[(0, 0)], 58.0);   // 1*7 + 2*9 + 3*11
/// assert_eq!(c[(0, 1)], 64.0);   // 1*8 + 2*10 + 3*12
/// assert_eq!(c[(1, 0)], 139.0);  // 4*7 + 5*9 + 6*11
/// assert_eq!(c[(1, 1)], 154.0);  // 4*8 + 5*10 + 6*12
/// ```
pub fn multiply<T: Scalar, const M: usize, const N: usize, const P: usize>(
    a: &Matrix<T, M, N>,
    b: &Matrix<T, N, P>,
) -> Matrix<T, M, P> {
    let mut result = Matrix::<T, M, P>::zero();

    for i in 0..M {
        for j in 0..P {
            let mut sum = T::zero();
            for k in 0..N {
                sum += a[(i, k)] * b[(k, j)];
            }
            result[(i, j)] = sum;
        }
    }

    result
}

/// Multiplies a matrix by a column vector.
///
/// # Type Signature
///
/// ```text
/// Matrix<T, M, N> × Vector<T, N> → Vector<T, M>
/// ```
///
/// # Mathematical Definition
///
/// For matrix A (M×N) and vector v (N×1), the product is a vector w (M×1):
/// ```text
/// w_i = Σ(j=1 to N) A_ij · v_j
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, matrix_vector_multiply};
/// use axiomath::vector::Vector;
///
/// let m = Matrix::<f64, 2, 3>::from_rows([
///     [1.0, 2.0, 3.0],
///     [4.0, 5.0, 6.0],
/// ]);
///
/// let v = Vector::<f64, 3>::from([7.0, 8.0, 9.0]);
///
/// let result = matrix_vector_multiply(&m, &v);
///
/// assert_eq!(result[0], 50.0);  // 1*7 + 2*8 + 3*9
/// assert_eq!(result[1], 122.0); // 4*7 + 5*8 + 6*9
/// ```
pub fn matrix_vector_multiply<T: Scalar, const M: usize, const N: usize>(
    matrix: &Matrix<T, M, N>,
    vector: &Vector<T, N>,
) -> Vector<T, M> {
    let mut result = Vector::<T, M>::zero();

    for i in 0..M {
        let mut sum = T::zero();
        for j in 0..N {
            sum += matrix[(i, j)] * vector[j];
        }
        result[i] = sum;
    }

    result
}

/// Multiplies a row vector by a matrix.
///
/// # Type Signature
///
/// ```text
/// Vector<T, M> × Matrix<T, M, N> → Vector<T, N>
/// ```
///
/// # Mathematical Definition
///
/// For row vector v (1×M) and matrix A (M×N), the product is a row vector w (1×N):
/// ```text
/// w_j = Σ(i=1 to M) v_i · A_ij
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, vector_matrix_multiply};
/// use axiomath::vector::Vector;
///
/// let v = Vector::<f64, 2>::from([1.0, 2.0]);
///
/// let m = Matrix::<f64, 2, 3>::from_rows([
///     [3.0, 4.0, 5.0],
///     [6.0, 7.0, 8.0],
/// ]);
///
/// let result = vector_matrix_multiply(&v, &m);
///
/// assert_eq!(result[0], 15.0);  // 1*3 + 2*6
/// assert_eq!(result[1], 18.0);  // 1*4 + 2*7
/// assert_eq!(result[2], 21.0);  // 1*5 + 2*8
/// ```
pub fn vector_matrix_multiply<T: Scalar, const M: usize, const N: usize>(
    vector: &Vector<T, M>,
    matrix: &Matrix<T, M, N>,
) -> Vector<T, N> {
    let mut result = Vector::<T, N>::zero();

    for j in 0..N {
        let mut sum = T::zero();
        for i in 0..M {
            sum += vector[i] * matrix[(i, j)];
        }
        result[j] = sum;
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multiply_2x3_3x2() {
        let a = Matrix::<f64, 2, 3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);

        let b = Matrix::<f64, 3, 2>::from_rows([[7.0, 8.0], [9.0, 10.0], [11.0, 12.0]]);

        let c = multiply(&a, &b);

        assert_eq!(c.rows(), 2);
        assert_eq!(c.cols(), 2);

        // [1*7 + 2*9 + 3*11,  1*8 + 2*10 + 3*12]
        // [4*7 + 5*9 + 6*11,  4*8 + 5*10 + 6*12]
        assert_eq!(c[(0, 0)], 58.0);
        assert_eq!(c[(0, 1)], 64.0);
        assert_eq!(c[(1, 0)], 139.0);
        assert_eq!(c[(1, 1)], 154.0);
    }

    #[test]
    fn test_multiply_3x2_2x3() {
        let a = Matrix::<f64, 3, 2>::from_rows([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]);

        let b = Matrix::<f64, 2, 3>::from_rows([[7.0, 8.0, 9.0], [10.0, 11.0, 12.0]]);

        let c = multiply(&a, &b);

        assert_eq!(c.rows(), 3);
        assert_eq!(c.cols(), 3);

        assert_eq!(c[(0, 0)], 27.0); // 1*7 + 2*10
        assert_eq!(c[(0, 1)], 30.0); // 1*8 + 2*11
        assert_eq!(c[(0, 2)], 33.0); // 1*9 + 2*12
        assert_eq!(c[(1, 0)], 61.0); // 3*7 + 4*10
        assert_eq!(c[(2, 2)], 117.0); // 5*9 + 6*12
    }

    #[test]
    fn test_multiply_square_2x2() {
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);

        let b = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);

        let c = multiply(&a, &b);

        assert_eq!(c[(0, 0)], 19.0); // 1*5 + 2*7
        assert_eq!(c[(0, 1)], 22.0); // 1*6 + 2*8
        assert_eq!(c[(1, 0)], 43.0); // 3*5 + 4*7
        assert_eq!(c[(1, 1)], 50.0); // 3*6 + 4*8
    }

    #[test]
    fn test_multiply_identity_left() {
        let identity = Matrix::<f64, 3, 3>::identity();
        let a = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ]);

        let result = multiply(&identity, &a);

        assert_eq!(result, a);
    }

    #[test]
    fn test_multiply_identity_right() {
        let a = Matrix::<f64, 3, 3>::from_rows([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ]);
        let identity = Matrix::<f64, 3, 3>::identity();

        let result = multiply(&a, &identity);

        assert_eq!(result, a);
    }

    #[test]
    fn test_multiply_by_zero() {
        let a = Matrix::<f64, 2, 3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);

        let zero = Matrix::<f64, 3, 2>::zero();

        let result = multiply(&a, &zero);

        assert!(result.is_zero());
    }

    #[test]
    fn test_multiply_associative() {
        // (AB)C = A(BC)
        let a = Matrix::<f64, 2, 3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);

        let b = Matrix::<f64, 3, 2>::from_rows([[7.0, 8.0], [9.0, 10.0], [11.0, 12.0]]);

        let c = Matrix::<f64, 2, 2>::from_rows([[13.0, 14.0], [15.0, 16.0]]);

        let left = multiply(&multiply(&a, &b), &c);
        let right = multiply(&a, &multiply(&b, &c));

        assert_eq!(left, right);
    }

    #[test]
    fn test_multiply_distributive_left() {
        // A(B + C) = AB + AC
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);

        let b = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);

        let c = Matrix::<f64, 2, 2>::from_rows([[9.0, 10.0], [11.0, 12.0]]);

        let left = multiply(&a, &(b + c));
        let right = multiply(&a, &b) + multiply(&a, &c);

        assert_eq!(left, right);
    }

    #[test]
    fn test_multiply_distributive_right() {
        // (A + B)C = AC + BC
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);

        let b = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);

        let c = Matrix::<f64, 2, 2>::from_rows([[9.0, 10.0], [11.0, 12.0]]);

        let left = multiply(&(a + b), &c);
        let right = multiply(&a, &c) + multiply(&b, &c);

        assert_eq!(left, right);
    }

    #[test]
    fn test_multiply_not_commutative() {
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);

        let b = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);

        let ab = multiply(&a, &b);
        let ba = multiply(&b, &a);

        assert_ne!(ab, ba);
    }

    #[test]
    fn test_matrix_vector_multiply() {
        let m = Matrix::<f64, 2, 3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);

        let v = Vector::<f64, 3>::from([7.0, 8.0, 9.0]);

        let result = matrix_vector_multiply(&m, &v);

        assert_eq!(result[0], 50.0); // 1*7 + 2*8 + 3*9
        assert_eq!(result[1], 122.0); // 4*7 + 5*8 + 6*9
    }

    #[test]
    fn test_matrix_vector_multiply_identity() {
        let identity = Matrix::<f64, 3, 3>::identity();
        let v = Vector::<f64, 3>::from([1.0, 2.0, 3.0]);

        let result = matrix_vector_multiply(&identity, &v);

        assert_eq!(result, v);
    }

    #[test]
    fn test_matrix_vector_multiply_zero() {
        let m = Matrix::<f64, 2, 3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);

        let zero = Vector::<f64, 3>::zero();

        let result = matrix_vector_multiply(&m, &zero);

        assert!(result.is_zero());
    }

    #[test]
    fn test_vector_matrix_multiply() {
        let v = Vector::<f64, 2>::from([1.0, 2.0]);

        let m = Matrix::<f64, 2, 3>::from_rows([[3.0, 4.0, 5.0], [6.0, 7.0, 8.0]]);

        let result = vector_matrix_multiply(&v, &m);

        assert_eq!(result[0], 15.0); // 1*3 + 2*6
        assert_eq!(result[1], 18.0); // 1*4 + 2*7
        assert_eq!(result[2], 21.0); // 1*5 + 2*8
    }

    #[test]
    fn test_vector_matrix_multiply_identity() {
        let v = Vector::<f64, 3>::from([1.0, 2.0, 3.0]);
        let identity = Matrix::<f64, 3, 3>::identity();

        let result = vector_matrix_multiply(&v, &identity);

        assert_eq!(result, v);
    }

    #[test]
    fn test_multiply_scalar_compatibility() {
        // (λA)B = λ(AB) = A(λB)
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);

        let b = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);

        let scalar = 2.5;

        let left = multiply(&(a * scalar), &b);
        let middle = multiply(&a, &b) * scalar;
        let right = multiply(&a, &(b * scalar));

        assert_eq!(left, middle);
        assert_eq!(left, right);
    }
}
