//! Special mathematical functions.
//!
//! This module implements special functions that appear frequently in mathematics,
//! physics, statistics, and engineering. These functions extend beyond elementary
//! functions and often arise as solutions to differential equations or as integrals
//! that cannot be expressed in terms of elementary functions.
//!
//! # Contents
//!
//! - **Error Function Family**: erf, erfc, and their inverses
//! - **Gamma Function Family**: Γ(x), ln Γ(x), and incomplete gamma functions
//! - **Beta Function Family**: B(a,b), ln B(a,b), and incomplete beta functions
//!
//! # Applications
//!
//! These special functions are fundamental to:
//! - Probability distributions (normal, gamma, beta, chi-squared, t, F)
//! - Statistical inference and hypothesis testing
//! - Quantum mechanics and thermodynamics
//! - Signal processing and communications
//! - Numerical analysis and approximation theory
//!
//! # Implementation Notes
//!
//! All functions are implemented from first principles using:
//! - Series expansions (Taylor, Laurent)
//! - Continued fractions
//! - Asymptotic expansions
//! - Rational approximations (Lanczos, Chebyshev)
//!
//! Implementations prioritize:
//! 1. **Numerical stability** - avoiding overflow/underflow
//! 2. **Accuracy** - typically 10-15 digits of precision
//! 3. **Efficiency** - fast convergence and minimal iterations
//!
//! # Examples
//!
//! ```
//! use axiomath::special::{error_function, gamma_function, beta_function};
//!
//! // Error function at x=1
//! let erf_1 = error_function(1.0);
//! assert!((erf_1 - 0.8427007929497149).abs() < 1e-10);
//!
//! // Gamma function: Γ(5) = 4! = 24
//! let gamma_5 = gamma_function(5.0);
//! assert!((gamma_5 - 24.0).abs() < 1e-10);
//!
//! // Beta function: B(2,3) = 1/12
//! let beta_23 = beta_function(2.0, 3.0);
//! assert!((beta_23 - 1.0/12.0).abs() < 1e-10);
//! ```
//!
//! # References
//!
//! - Abramowitz, M., and Stegun, I. A. (1964). "Handbook of Mathematical Functions"
//! - Press, W. H., et al. (2007). "Numerical Recipes: The Art of Scientific Computing"
//! - NIST Digital Library of Mathematical Functions: https://dlmf.nist.gov/

// Error function family
pub mod error_function;
pub use error_function::{
    complementary_error_function,
    error_function,
    inverse_error_function,
    scaled_complementary_error_function,
};

// Gamma function family
pub mod gamma;
pub use gamma::{
    gamma_function,
    incomplete_gamma_lower,
    log_gamma_function,
    regularized_incomplete_gamma,
};

// Beta function family
pub mod beta;
pub use beta::{
    beta_function,
    incomplete_beta,
    log_beta_function,
    regularized_incomplete_beta,
};
