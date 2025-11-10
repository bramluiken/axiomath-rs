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
//!
//! // Calculus operations
//! let f = |x: f64| x.sin();
//! let derivative = calculus::central_difference(f, 0.0, 1e-5);
//! let integral = calculus::simpsons_rule(f, 0.0, std::f64::consts::PI, 100);
//!
//! // ODE solving
//! let ode = |_t: f64, y: f64| y;
//! let solution = calculus::runge_kutta_4(ode, 1.0, (0.0, 1.0), 0.01);
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

// Calculus operations
pub mod calculus {
    //! Calculus operations namespace.
    //!
    //! This submodule provides numerical differentiation, integration, and ODE solvers.

    // Differentiation
    pub mod differentiation {
        //! Numerical differentiation methods.
        pub use crate::calculus::differentiation::{
            forward_difference, backward_difference, central_difference,
            five_point_stencil, derivative_at_point, second_derivative, nth_derivative,
        };
    }

    // Integration
    pub mod integration {
        //! Numerical integration (quadrature) methods.
        pub use crate::calculus::integration::{
            // Newton-Cotes
            trapezoidal_rule, simpsons_rule, simpsons_3_8_rule,
            // Adaptive
            adaptive_simpson, romberg_integration,
            // Gaussian quadrature
            gauss_legendre_quadrature, gauss_laguerre_quadrature, gauss_hermite_quadrature,
        };
    }

    // ODE solvers
    pub mod ode {
        //! Ordinary differential equation solvers.
        pub use crate::calculus::ode::{
            euler_method, runge_kutta_4, adaptive_rk45,
        };
    }

    // Flat re-exports for common operations
    pub use crate::calculus::differentiation::central_difference;
    pub use crate::calculus::integration::{simpsons_rule, adaptive_simpson, gauss_legendre_quadrature};
    pub use crate::calculus::ode::{runge_kutta_4, adaptive_rk45};
}

// Complex numbers
pub mod complex {
    //! Complex number operations namespace.
    //!
    //! This submodule provides complex number arithmetic, polar form operations,
    //! and complex-valued functions.

    // Type
    pub use crate::complex::Complex;

    // Arithmetic
    pub use crate::complex::{
        add, subtract, multiply, divide, negate, conjugate, reciprocal,
    };

    // Polar form
    pub use crate::complex::{
        magnitude, argument, to_polar, from_polar, multiply_polar,
    };

    // Complex functions
    pub use crate::complex::{
        complex_exp, complex_ln, complex_sqrt, complex_pow, complex_sin, complex_cos,
    };
}

// Polynomials
pub mod polynomial {
    //! Polynomial operations namespace.
    //!
    //! This submodule provides polynomial arithmetic, evaluation, calculus operations,
    //! and root-finding algorithms.

    // Type
    pub use crate::polynomial::Polynomial;

    // Arithmetic
    pub use crate::polynomial::{
        add, subtract, multiply, negate,
    };

    // Evaluation
    pub use crate::polynomial::{
        evaluate, evaluate_complex,
    };

    // Calculus
    pub use crate::polynomial::{
        derivative, integral, definite_integral,
    };

    // Properties
    pub use crate::polynomial::{
        discriminant,
    };

    // Root finding
    pub use crate::polynomial::{
        linear_root, quadratic_roots, cubic_root,
        newton_raphson_root, durand_kerner,
    };
}

// Statistics operations
pub mod statistics {
    //! Statistics and probability distributions namespace.
    //!
    //! This submodule provides descriptive statistics and probability distributions.

    // Error type
    pub use crate::statistics::StatisticsError;

    // Descriptive statistics
    pub mod descriptive {
        //! Descriptive statistics functions.
        pub use crate::statistics::descriptive::{
            // Central tendency
            mean, weighted_mean, median, mode, geometric_mean,
            // Dispersion
            variance, standard_deviation, mean_absolute_deviation, range, coefficient_of_variation,
            // Shape
            skewness, kurtosis,
            // Order statistics
            percentile, quantile, quartiles, interquartile_range, five_number_summary,
            // Multivariate
            covariance, correlation,
        };
    }

    // Probability distributions
    pub mod distributions {
        //! Probability distribution functions.

        // Discrete distributions
        pub mod discrete {
            //! Discrete probability distributions.
            pub use crate::statistics::distributions::discrete::{
                bernoulli_pmf, bernoulli_cdf, bernoulli_mean, bernoulli_variance,
                binomial_pmf, binomial_cdf, binomial_mean, binomial_variance,
                poisson_pmf, poisson_cdf, poisson_mean, poisson_variance,
                geometric_pmf, geometric_cdf, geometric_mean, geometric_variance,
            };
        }

        // Continuous distributions
        pub mod continuous {
            //! Continuous probability distributions.
            pub use crate::statistics::distributions::continuous::{
                normal_pdf, normal_cdf, normal_mean, normal_variance,
                uniform_pdf, uniform_cdf, uniform_mean, uniform_variance,
                exponential_pdf, exponential_cdf, exponential_mean, exponential_variance,
                gamma_pdf, gamma_cdf, gamma_mean, gamma_variance,
                beta_pdf, beta_cdf, beta_mean, beta_variance,
                chi_squared_pdf, chi_squared_cdf, chi_squared_mean, chi_squared_variance,
                students_t_pdf, students_t_cdf, students_t_mean, students_t_variance,
                f_pdf, f_cdf, f_mean, f_variance,
            };
        }

        // Flat re-exports for common distributions
        pub use crate::statistics::distributions::{
            binomial_pmf, poisson_pmf, normal_pdf, normal_cdf, exponential_pdf,
        };
    }

    // Flat re-exports for common operations
    pub use crate::statistics::descriptive::{mean, standard_deviation, variance, correlation};
    pub use crate::statistics::distributions::{binomial_pmf, normal_pdf, normal_cdf};
}

// Re-export num-traits for convenience
pub use num_traits::{One, Zero};
