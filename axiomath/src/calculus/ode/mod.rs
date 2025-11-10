//! Ordinary Differential Equation (ODE) solvers.
//!
//! This module implements numerical methods for solving initial value problems
//! of the form:
//!
//! ```text
//! dy/dt = f(t, y),  y(t₀) = y₀
//! ```
//!
//! # Theory
//!
//! ODE solvers compute approximate solutions by discretizing the time domain
//! and updating the solution stepwise. Methods differ in:
//!
//! - **Accuracy order**: How fast the error decreases with step size
//! - **Stability**: Ability to handle stiff equations
//! - **Efficiency**: Computational cost per step
//!
//! ## Methods
//!
//! ### Euler Method
//!
//! The simplest method with first-order accuracy:
//! ```text
//! yₙ₊₁ = yₙ + h·f(tₙ, yₙ)
//! ```
//!
//! ### Runge-Kutta Methods
//!
//! Higher-order methods that evaluate f at multiple points:
//! - **RK4**: Fourth-order accuracy, widely used
//! - **RK45**: Adaptive step size with error control
//!
//! # Examples
//!
//! ## Exponential Growth
//!
//! ```
//! use axiomath::calculus::ode::runge_kutta_4;
//!
//! // Solve dy/dt = y with y(0) = 1, solution is y(t) = e^t
//! let f = |_t: f64, y: f64| y;
//! let solution = runge_kutta_4(f, 1.0, (0.0, 1.0), 0.01);
//!
//! let y_at_1 = solution.last().unwrap().1;
//! assert!((y_at_1 - std::f64::consts::E).abs() < 1e-6);
//! ```
//!
//! ## Harmonic Oscillator
//!
//! For systems, use multiple calls or extend to vector-valued functions:
//! ```text
//! d²x/dt² = -x  →  { dx/dt = v, dv/dt = -x }
//! ```
//!
//! # References
//!
//! - Hairer, Ernst, et al. "Solving Ordinary Differential Equations I." 2nd ed., 1993.
//! - Burden & Faires, "Numerical Analysis", Chapter 5.
//! - Press et al., "Numerical Recipes", Chapter 16.

mod adaptive_rk;
mod euler;
mod rk4;

// Re-export public functions
pub use adaptive_rk::adaptive_rk45;
pub use euler::euler_method;
pub use rk4::runge_kutta_4;
