//! Prelude module providing commonly used types and traits.
//!
//! This module re-exports the most commonly used items from axiomath,
//! allowing users to quickly import everything they need with a single use statement.
//!
//! # Examples
//!
//! ```
//! use axiomath::prelude::*;
//!
//! // Core traits and constants
//! let x: f64 = 2.0;
//! let y = x.sqrt();
//!
//! // Vector operations
//! let v1 = vec3(1.0, 2.0, 3.0);
//! let v2 = vec3(4.0, 5.0, 6.0);
//! let dot = dot_product(&v1, &v2);
//! let mag = magnitude(&v1);
//!
//! // Matrix operations
//! let m1 = Mat2::<f64>::identity();
//! let m2 = Mat2::<f64>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
//! let product = matrix::multiply(&m1, &m2);
//! let det = matrix::determinant(&m2);
//! ```

// Core traits
pub use crate::core::{Float, Real, Scalar};
pub use crate::core::constants::*;

// Vector types and constructors
pub use crate::vector::{Vector, Vec2, Vec3, Vec4, vec2, vec3, vec4};

// Common vector operations
pub use crate::vector::{
    // Arithmetic
    add as vector_add,
    subtract as vector_subtract,
    multiply as vector_multiply,
    divide as vector_divide,
    negate as vector_negate,
    reciprocal as vector_reciprocal,
    // Scalar operations
    scale as vector_scale,
    divide_by_scalar as vector_divide_by_scalar,
    // Norms
    magnitude, magnitude_squared, normalize, distance, euclidean_norm,
    // Products
    dot_product, cross_product,
    // Angles
    angle_between,
};

// Matrix types and constructors
pub use crate::matrix::{Matrix, Mat2, Mat3, Mat4, Mat2x3, Mat3x2};

// Matrix operations (re-exported in matrix namespace to avoid naming conflicts)
pub mod matrix {
    //! Matrix operations namespace.
    //!
    //! This submodule groups matrix operations to avoid naming conflicts
    //! with vector operations.

    pub use crate::matrix::{
        // Arithmetic
        add, subtract, scale, divide_by_scalar, negate,
        multiply_elementwise, divide_elementwise,
        // Multiplication
        multiply, matrix_vector_multiply, vector_matrix_multiply,
        // Transpose
        transpose,
        // Properties
        trace, determinant, determinant_2x2, determinant_3x3, determinant_4x4,
        is_symmetric, is_diagonal, is_upper_triangular, is_lower_triangular,
        is_singular, is_invertible,
    };
}

// Re-export num-traits for convenience
pub use num_traits::{One, Zero};
