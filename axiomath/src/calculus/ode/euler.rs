//! Euler's method for solving ODEs.
//!
//! The simplest explicit method for solving ordinary differential equations.

use crate::core::Float;

/// Solves an ODE using Euler's method.
///
/// # Algorithm
///
/// Euler's method is a first-order numerical procedure for solving ODEs with
/// a given initial value:
///
/// ```text
/// yₙ₊₁ = yₙ + h·f(tₙ, yₙ)
/// ```
///
/// where h is the step size and f(t,y) is the derivative dy/dt.
///
/// # Accuracy
///
/// - Local truncation error: O(h²)
/// - Global error: O(h)
///
/// This means doubling the number of steps (halving h) roughly halves the error.
///
/// # Stability
///
/// Euler's method is conditionally stable. For the test equation dy/dt = λy,
/// it is stable when |1 + hλ| ≤ 1.
///
/// # Arguments
///
/// * `f` - The derivative function f(t, y) where dy/dt = f(t, y)
/// * `y0` - Initial condition y(t₀) = y₀
/// * `t_range` - Time interval (t₀, t_final)
/// * `dt` - Step size (must be positive)
///
/// # Returns
///
/// A vector of (time, value) pairs representing the solution trajectory.
///
/// # Panics
///
/// Panics if `dt` is not positive or if `t_range` is invalid.
///
/// # Examples
///
/// ```
/// use axiomath::calculus::ode::euler_method;
///
/// // Solve dy/dt = -y with y(0) = 1, exact solution: y(t) = e^(-t)
/// let f = |_t: f64, y: f64| -y;
/// let solution = euler_method(f, 1.0, (0.0, 2.0), 0.01);
///
/// let y_at_2 = solution.last().unwrap().1;
/// let exact = (-2.0_f64).exp();
/// assert!((y_at_2 - exact).abs() < 0.01);
/// ```
///
/// ```
/// use axiomath::calculus::ode::euler_method;
///
/// // Solve dy/dt = t with y(0) = 0, exact solution: y(t) = t²/2
/// let f = |t: f64, _y: f64| t;
/// let solution = euler_method(f, 0.0, (0.0, 1.0), 0.001);
///
/// let y_at_1 = solution.last().unwrap().1;
/// assert!((y_at_1 - 0.5).abs() < 0.001);
/// ```
///
/// # References
///
/// - Euler, Leonhard. "Institutionum calculi integralis." 1768-1770.
/// - Burden & Faires, "Numerical Analysis", Section 5.2.
pub fn euler_method<T, F>(f: F, y0: T, t_range: (T, T), dt: T) -> Vec<(T, T)>
where
    T: Float,
    F: Fn(T, T) -> T,
{
    let (t0, t_final) = t_range;

    assert!(dt > T::zero(), "Step size must be positive");
    assert!(t_final >= t0, "Final time must be >= initial time");

    // Calculate number of steps
    let n_steps = ((t_final - t0) / dt).ceil().to_usize().unwrap();

    // Initialize solution vector
    let mut solution = Vec::with_capacity(n_steps + 1);
    solution.push((t0, y0));

    let mut t = t0;
    let mut y = y0;

    for _ in 0..n_steps {
        // Euler step: y_{n+1} = y_n + h * f(t_n, y_n)
        let dy = f(t, y);
        y = y + dt * dy;
        t = t + dt;

        // Store the solution point
        solution.push((t, y));

        // Stop if we've reached or passed t_final
        if t >= t_final {
            break;
        }
    }

    solution
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_euler_exponential_growth() {
        // dy/dt = y, y(0) = 1
        // Exact solution: y(t) = e^t
        let f = |_t: f64, y: f64| y;
        let solution = euler_method(f, 1.0, (0.0, 1.0), 0.001);

        let y_at_1 = solution.last().unwrap().1;
        let exact = std::f64::consts::E;

        // Euler has O(h) global error
        assert_abs_diff_eq!(y_at_1, exact, epsilon = 0.002);
    }

    #[test]
    fn test_euler_exponential_decay() {
        // dy/dt = -y, y(0) = 1
        // Exact solution: y(t) = e^(-t)
        let f = |_t: f64, y: f64| -y;
        let solution = euler_method(f, 1.0, (0.0, 2.0), 0.01);

        let y_at_2 = solution.last().unwrap().1;
        let exact = (-2.0_f64).exp();

        assert_abs_diff_eq!(y_at_2, exact, epsilon = 0.02);
    }

    #[test]
    fn test_euler_linear_ode() {
        // dy/dt = t, y(0) = 0
        // Exact solution: y(t) = t²/2
        let f = |t: f64, _y: f64| t;
        let solution = euler_method(f, 0.0, (0.0, 1.0), 0.001);

        let y_at_1 = solution.last().unwrap().1;
        let exact = 0.5;

        assert_abs_diff_eq!(y_at_1, exact, epsilon = 0.001);
    }

    #[test]
    fn test_euler_constant_solution() {
        // dy/dt = 0, y(0) = 5
        // Exact solution: y(t) = 5
        let f = |_t: f64, _y: f64| 0.0;
        let solution = euler_method(f, 5.0, (0.0, 1.0), 0.1);

        for &(_t, y) in &solution {
            assert_abs_diff_eq!(y, 5.0, epsilon = 1e-15);
        }
    }

    #[test]
    fn test_euler_linear_growth() {
        // dy/dt = 2, y(0) = 1
        // Exact solution: y(t) = 1 + 2t
        let f = |_t: f64, _y: f64| 2.0;
        let solution = euler_method(f, 1.0, (0.0, 1.0), 0.01);

        let y_at_1 = solution.last().unwrap().1;
        let exact = 3.0;

        assert_abs_diff_eq!(y_at_1, exact, epsilon = 1e-10);
    }

    #[test]
    fn test_euler_solution_length() {
        let f = |_t: f64, y: f64| y;
        let dt = 0.1;
        let solution = euler_method(f, 1.0, (0.0, 1.0), dt);

        // Should have approximately (1.0 - 0.0) / 0.1 + 1 = 11 points
        assert!(solution.len() >= 11);
        assert!(solution.len() <= 12); // Allow for floating-point rounding
    }

    #[test]
    fn test_euler_time_values() {
        let f = |_t: f64, _y: f64| 0.0;
        let dt = 0.1;
        let solution = euler_method(f, 0.0, (0.0, 1.0), dt);

        // Check that time values are correct
        for (i, &(t, _y)) in solution.iter().enumerate() {
            let expected_t = (i as f64) * dt;
            assert_abs_diff_eq!(t, expected_t, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_euler_convergence() {
        // Verify first-order convergence: halving h should approximately halve error
        let f = |_t: f64, y: f64| y;
        let exact = std::f64::consts::E;

        let sol_h1 = euler_method(f, 1.0, (0.0, 1.0), 0.1);
        let sol_h2 = euler_method(f, 1.0, (0.0, 1.0), 0.05);

        let error_h1 = (sol_h1.last().unwrap().1 - exact).abs();
        let error_h2 = (sol_h2.last().unwrap().1 - exact).abs();

        // For first-order method, error_h2 should be roughly error_h1/2
        let ratio = error_h1 / error_h2;
        assert!(ratio > 1.8 && ratio < 2.2);
    }

    #[test]
    fn test_euler_zero_interval() {
        let f = |_t: f64, y: f64| y;
        let solution = euler_method(f, 1.0, (0.0, 0.0), 0.1);

        // Should return just the initial condition
        assert_eq!(solution.len(), 1);
        assert_eq!(solution[0], (0.0, 1.0));
    }

    #[test]
    fn test_euler_f32_support() {
        // Test that method works with f32
        let f = |_t: f32, y: f32| y;
        let solution = euler_method(f, 1.0_f32, (0.0_f32, 1.0_f32), 0.01_f32);

        let y_at_1 = solution.last().unwrap().1;
        let exact = std::f32::consts::E;

        assert_abs_diff_eq!(y_at_1, exact, epsilon = 0.02_f32);
    }

    #[test]
    #[should_panic(expected = "Step size must be positive")]
    fn test_euler_negative_step() {
        let f = |_t: f64, y: f64| y;
        euler_method(f, 1.0, (0.0, 1.0), -0.1);
    }

    #[test]
    #[should_panic(expected = "Step size must be positive")]
    fn test_euler_zero_step() {
        let f = |_t: f64, y: f64| y;
        euler_method(f, 1.0, (0.0, 1.0), 0.0);
    }

    #[test]
    #[should_panic(expected = "Final time must be >= initial time")]
    fn test_euler_invalid_time_range() {
        let f = |_t: f64, y: f64| y;
        euler_method(f, 1.0, (1.0, 0.0), 0.1);
    }
}
