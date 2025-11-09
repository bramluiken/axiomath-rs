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
//!
//! // Linear algebra
//! let a = Mat3::<f64>::from_rows([
//!     [1.0, 2.0, 3.0],
//!     [0.0, 1.0, 4.0],
//!     [5.0, 6.0, 0.0],
//! ]);
//! let b = vec3(1.0, 2.0, 3.0);
//! let x = linear_algebra::solve(&a, &b).unwrap();
//! let a_inv = linear_algebra::invert(&a).unwrap();
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

// Linear algebra operations
pub mod linear_algebra {
    //! Linear algebra operations namespace.
    //!
    //! This submodule provides decompositions, solvers, and vector space operations.

    // Error type
    pub use crate::linear_algebra::LinearAlgebraError;

    // Decompositions
    pub mod decomposition {
        //! Matrix decomposition algorithms.
        pub use crate::linear_algebra::decomposition::{
            lu_decompose, qr_decompose, qr_decompose_householder,
            cholesky::{cholesky_decompose, cholesky_solve},
            svd::{svd, pseudoinverse},
            eigen::{eigenvalues, eigenpairs, power_iteration},
        };
    }

    // System solvers
    pub mod systems {
        //! Linear system solvers and matrix inversion.
        pub use crate::linear_algebra::systems::{
            // Triangular systems
            triangular::{solve_lower_triangular, solve_upper_triangular},
            // Direct solver
            direct::{solve_linear_system, gaussian_elimination},
            // Matrix inversion
            inverse::invert,
            // Iterative solvers
            iterative::{jacobi, gauss_seidel, conjugate_gradient},
        };
    }

    // Vector spaces
    pub mod spaces {
        //! Vector space and subspace operations.
        pub use crate::linear_algebra::spaces::{
            // Basis operations
            basis::{gram_schmidt, is_linearly_independent, span_dimension, find_basis},
            // Subspace operations
            subspaces::{column_space, row_space, null_space, left_null_space},
        };
    }

    // Row operations and RREF
    pub mod rref {
        //! Reduced Row Echelon Form operations.
        pub use crate::linear_algebra::rref::{rref, rref_with_pivots};
        pub use crate::linear_algebra::row_operations::{swap_rows, scale_row, add_scaled_row};
    }

    // Flat re-exports for common operations
    pub use crate::linear_algebra::decomposition::{lu_decompose, qr_decompose, svd::svd};
    pub use crate::linear_algebra::systems::direct::solve_linear_system as solve;
    pub use crate::linear_algebra::systems::inverse::invert;
}

// Re-export num-traits for convenience
pub use num_traits::{One, Zero};
