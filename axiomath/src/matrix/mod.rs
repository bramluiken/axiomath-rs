//! Matrix types and operations.
//!
//! This module provides generic matrix types using const generics for compile-time
//! dimension safety, along with all fundamental matrix operations implemented from
//! first principles.
//!
//! # Core Type
//!
//! The [`Matrix<T, ROWS, COLS>`](Matrix) type uses const generics to enforce dimension
//! correctness at compile time. This prevents common errors like multiplying
//! incompatible matrices.
//!
//! # Type Aliases
//!
//! Common matrix sizes have convenient type aliases:
//! - `Mat2<T>`, `Mat3<T>`, `Mat4<T>` for square matrices
//! - `Mat2x3<T>`, `Mat3x4<T>`, etc. for rectangular matrices
//!
//! # Operations
//!
//! - **Arithmetic**: Element-wise addition, subtraction, scalar multiplication
//! - **Matrix multiplication**: Type-safe standard matrix multiplication
//! - **Transpose**: Convert M×N to N×M
//! - **Properties**: Trace, determinant, symmetry checks, etc.
//!
//! # Examples
//!
//! ```
//! use axiomath::matrix::{Matrix, Mat2, Mat3, multiply, transpose, determinant};
//!
//! // Create matrices
//! let a = Mat2::<f64>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
//! let b = Mat2::<f64>::identity();
//!
//! // Arithmetic operations
//! let sum = a + b;
//! let scaled = a * 2.0;
//!
//! // Matrix multiplication
//! let product = multiply(&a, &b);
//!
//! // Properties
//! let det = determinant(&a);  // -2.0
//! let at = transpose(&a);
//! ```
//!
//! # Design Philosophy
//!
//! All operations are implemented from first principles using mathematical
//! definitions. The type system ensures dimension safety at compile time,
//! preventing many common linear algebra errors.

mod arithmetic;
mod constructors;
mod matrix_type;
mod multiplication;
mod properties;
mod transpose;

// Re-export core type
pub use matrix_type::Matrix;

// Re-export arithmetic operations
pub use arithmetic::{
    add, divide_by_scalar, divide_elementwise, multiply_elementwise, negate, scale, subtract,
};

// Re-export multiplication operations
pub use multiplication::{matrix_vector_multiply, multiply, vector_matrix_multiply};

// Re-export transpose operations
pub use transpose::transpose;

// Re-export properties
pub use properties::{
    determinant, determinant_2x2, determinant_3x3, determinant_4x4, is_diagonal, is_invertible,
    is_lower_triangular, is_singular, is_symmetric, is_upper_triangular, trace,
};

// Type aliases for square matrices
/// 2×2 square matrix.
pub type Mat2<T> = Matrix<T, 2, 2>;

/// 3×3 square matrix.
pub type Mat3<T> = Matrix<T, 3, 3>;

/// 4×4 square matrix.
pub type Mat4<T> = Matrix<T, 4, 4>;

// Type aliases for common rectangular matrices
/// 2×3 matrix (2 rows, 3 columns).
pub type Mat2x3<T> = Matrix<T, 2, 3>;

/// 2×4 matrix (2 rows, 4 columns).
pub type Mat2x4<T> = Matrix<T, 2, 4>;

/// 3×2 matrix (3 rows, 2 columns).
pub type Mat3x2<T> = Matrix<T, 3, 2>;

/// 3×4 matrix (3 rows, 4 columns).
pub type Mat3x4<T> = Matrix<T, 3, 4>;

/// 4×2 matrix (4 rows, 2 columns).
pub type Mat4x2<T> = Matrix<T, 4, 2>;

/// 4×3 matrix (4 rows, 3 columns).
pub type Mat4x3<T> = Matrix<T, 4, 3>;

#[cfg(test)]
mod integration_tests {
    use super::*;

    #[test]
    fn test_matrix_workflow() {
        // Create matrices
        let a = Mat3::<f64>::from_rows([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ]);

        let b = Mat3::<f64>::identity();

        // Arithmetic
        let sum = add(&a, &b);
        assert_eq!(sum[(0, 0)], 2.0); // 1 + 1
        assert_eq!(sum[(0, 1)], 2.0); // 2 + 0

        // Multiplication
        let product = multiply(&a, &b);
        assert_eq!(product, a); // A * I = A

        // Transpose
        let at = transpose(&a);
        assert_eq!(at[(0, 1)], 4.0); // Was a[(1, 0)]

        // Properties
        let tr = trace(&a);
        assert_eq!(tr, 15.0); // 1 + 5 + 9

        let det = determinant(&a);
        assert_eq!(det, 0.0); // Singular matrix
    }

    #[test]
    fn test_type_aliases() {
        let m2 = Mat2::<f64>::identity();
        assert_eq!(m2.rows(), 2);
        assert_eq!(m2.cols(), 2);

        let m3 = Mat3::<f64>::zero();
        assert_eq!(m3.rows(), 3);
        assert_eq!(m3.cols(), 3);

        let m2x3 = Mat2x3::<f64>::zero();
        assert_eq!(m2x3.rows(), 2);
        assert_eq!(m2x3.cols(), 3);
    }

    #[test]
    fn test_operator_overloading() {
        let a = Mat2::<f64>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let b = Mat2::<f64>::from_rows([[5.0, 6.0], [7.0, 8.0]]);

        // Addition
        let sum = a + b;
        assert_eq!(sum[(0, 0)], 6.0);

        // Subtraction
        let diff = a - b;
        assert_eq!(diff[(0, 0)], -4.0);

        // Negation
        let neg = -a;
        assert_eq!(neg[(0, 0)], -1.0);

        // Scalar multiplication
        let scaled = a * 2.0;
        assert_eq!(scaled[(0, 0)], 2.0);

        // Scalar division
        let divided = a / 2.0;
        assert_eq!(divided[(0, 0)], 0.5);
    }

    #[test]
    fn test_chained_operations() {
        let a = Mat2::<f64>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let b = Mat2::<f64>::from_rows([[5.0, 6.0], [7.0, 8.0]]);

        // (2A + B) * 0.5
        let result = (a * 2.0 + b) / 2.0;

        assert_eq!(result[(0, 0)], 3.5); // (2*1 + 5) / 2
        assert_eq!(result[(0, 1)], 5.0); // (2*2 + 6) / 2
        assert_eq!(result[(1, 0)], 6.5); // (2*3 + 7) / 2
        assert_eq!(result[(1, 1)], 8.0); // (2*4 + 8) / 2
    }

    #[test]
    fn test_matrix_vector_integration() {
        let m = Mat2x3::<f64>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);

        let v = crate::vector::Vector::<f64, 3>::from([1.0, 0.0, 1.0]);

        let result = matrix_vector_multiply(&m, &v);

        assert_eq!(result[0], 4.0); // 1*1 + 2*0 + 3*1
        assert_eq!(result[1], 10.0); // 4*1 + 5*0 + 6*1
    }

    #[test]
    fn test_mathematical_properties() {
        let a = Mat2::<f64>::from_rows([[1.0, 2.0], [3.0, 4.0]]);

        // Transpose twice returns original
        let at = transpose(&a);
        let att = transpose(&at);
        assert_eq!(att, a);

        // Determinant of transpose equals determinant
        let det_a = determinant(&a);
        let det_at = determinant(&at);
        assert_eq!(det_a, det_at);

        // Identity properties
        let identity = Mat2::<f64>::identity();
        assert_eq!(determinant(&identity), 1.0);
        assert!(is_symmetric(&identity));
        assert!(is_diagonal(&identity));
        assert!(is_upper_triangular(&identity));
        assert!(is_lower_triangular(&identity));
    }
}
