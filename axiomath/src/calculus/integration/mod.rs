//! Numerical integration (quadrature) methods.
//!
//! This module implements various techniques for approximating definite integrals.
//!
//! # Theory
//!
//! Numerical integration approximates the definite integral:
//!
//! ```text
//! ∫[a,b] f(x)dx ≈ Σᵢ wᵢ f(xᵢ)
//! ```
//!
//! where xᵢ are evaluation points and wᵢ are weights.
//!
//! # Methods
//!
//! ## Newton-Cotes Formulas
//!
//! Based on equally-spaced points:
//! - **Trapezoidal Rule**: Linear interpolation, O(h²) accuracy
//! - **Simpson's 1/3 Rule**: Quadratic interpolation, O(h⁴) accuracy
//! - **Simpson's 3/8 Rule**: Cubic interpolation, O(h⁴) accuracy
//!
//! ## Adaptive Methods
//!
//! Automatically adjust the integration intervals:
//! - **Adaptive Simpson**: Recursively subdivides intervals
//! - **Romberg Integration**: Richardson extrapolation for high accuracy
//!
//! ## Gaussian Quadrature
//!
//! Optimal placement of evaluation points:
//! - **Gauss-Legendre**: For finite intervals [a,b]
//! - **Gauss-Laguerre**: For semi-infinite [0,∞) with weight e^(-x)
//! - **Gauss-Hermite**: For infinite (-∞,∞) with weight e^(-x²)
//!
//! # Examples
//!
//! ```
//! use axiomath::calculus::integration::{simpsons_rule, trapezoidal_rule};
//!
//! // Integrate x² from 0 to 1, exact value is 1/3
//! let f = |x: f64| x * x;
//! let result_simpson = simpsons_rule(f, 0.0, 1.0, 100);
//! let result_trap = trapezoidal_rule(f, 0.0, 1.0, 100);
//!
//! assert!((result_simpson - 1.0/3.0).abs() < 1e-8);
//! assert!((result_trap - 1.0/3.0).abs() < 1e-5);
//! ```

mod adaptive;
mod gaussian;
mod newton_cotes;

// Re-export all public functions
pub use adaptive::{adaptive_simpson, romberg_integration};
pub use gaussian::{gauss_hermite_quadrature, gauss_laguerre_quadrature, gauss_legendre_quadrature};
pub use newton_cotes::{simpsons_3_8_rule, simpsons_rule, trapezoidal_rule};
