//! Calculus operations.
//!
//! This module implements numerical differentiation, integration, and
//! ordinary differential equation (ODE) solvers from first principles.
//!
//! # Theory
//!
//! Calculus is the mathematical study of continuous change. This module provides
//! numerical methods for computing derivatives, integrals, and solutions to
//! differential equations when analytical solutions are impractical or impossible.
//!
//! ## Differentiation
//!
//! Numerical differentiation approximates derivatives using finite differences:
//!
//! ```text
//! f'(x) ≈ [f(x+h) - f(x)] / h     (forward difference, O(h))
//! f'(x) ≈ [f(x+h) - f(x-h)] / 2h  (central difference, O(h²))
//! ```
//!
//! Higher-order methods provide better accuracy at the cost of more function evaluations.
//!
//! ## Integration
//!
//! Numerical integration (quadrature) approximates definite integrals:
//!
//! ```text
//! ∫[a,b] f(x)dx ≈ Σ wᵢ f(xᵢ)
//! ```
//!
//! Methods include:
//! - Newton-Cotes formulas (trapezoidal rule, Simpson's rule)
//! - Adaptive methods (adaptive Simpson, Romberg integration)
//! - Gaussian quadrature (optimal node placement)
//!
//! ## Ordinary Differential Equations
//!
//! ODE solvers find solutions to equations of the form:
//!
//! ```text
//! dy/dt = f(t, y),  y(t₀) = y₀
//! ```
//!
//! Methods include explicit methods (Euler, Runge-Kutta) and adaptive
//! step-size control for efficiency and accuracy.
//!
//! # Examples
//!
//! ## Differentiation
//!
//! ```
//! use axiomath::calculus::differentiation::central_difference;
//!
//! // Derivative of sin(x) at x=0 is cos(0)=1
//! let f = |x: f64| x.sin();
//! let derivative = central_difference(f, 0.0, 1e-5);
//! assert!((derivative - 1.0).abs() < 1e-8);
//! ```
//!
//! ## Integration
//!
//! ```
//! use axiomath::calculus::integration::simpsons_rule;
//!
//! // Integral of x² from 0 to 1 is 1/3
//! let f = |x: f64| x * x;
//! let integral = simpsons_rule(f, 0.0, 1.0, 100);
//! assert!((integral - 1.0/3.0).abs() < 1e-6);
//! ```
//!
//! ## ODE Solving
//!
//! ```
//! use axiomath::calculus::ode::runge_kutta_4;
//!
//! // Solve dy/dt = y with y(0) = 1, solution is y(t) = e^t
//! let f = |_t: f64, y: f64| y;
//! let solution = runge_kutta_4(f, 1.0, (0.0, 1.0), 0.01);
//! let y_at_1 = solution.last().unwrap().1;
//! assert!((y_at_1 - std::f64::consts::E).abs() < 1e-3);
//! ```
//!
//! # References
//!
//! - Press, William H., et al. "Numerical Recipes: The Art of Scientific Computing." 3rd ed., 2007.
//! - Burden, Richard L., and J. Douglas Faires. "Numerical Analysis." 9th ed., 2010.
//! - Hairer, Ernst, Syvert P. Nørsett, and Gerhard Wanner. "Solving Ordinary Differential Equations I." 2nd ed., 1993.

pub mod differentiation;
pub mod integration;
pub mod ode;

// Re-export commonly used functions
pub use differentiation::{
    backward_difference, central_difference, derivative_at_point, five_point_stencil,
    forward_difference, nth_derivative, second_derivative,
};

pub use integration::{
    adaptive_simpson, gauss_hermite_quadrature, gauss_laguerre_quadrature,
    gauss_legendre_quadrature, romberg_integration, simpsons_3_8_rule, simpsons_rule,
    trapezoidal_rule,
};

pub use ode::{adaptive_rk45, euler_method, runge_kutta_4};
