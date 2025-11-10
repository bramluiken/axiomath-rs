//! Runge-Kutta 4th order method for solving ODEs.
//!
//! The classic RK4 method is one of the most widely used ODE solvers, offering
//! a good balance between accuracy and computational cost.

use crate::core::Float;

/// Solves an ODE using the classical Runge-Kutta 4th order method (RK4).
///
/// # Algorithm
///
/// RK4 is a fourth-order accurate method that uses four function evaluations
/// per step to achieve high accuracy:
///
/// ```text
/// k₁ = f(tₙ, yₙ)
/// k₂ = f(tₙ + h/2, yₙ + h·k₁/2)
/// k₃ = f(tₙ + h/2, yₙ + h·k₂/2)
/// k₄ = f(tₙ + h, yₙ + h·k₃)
///
/// yₙ₊₁ = yₙ + (h/6)·(k₁ + 2k₂ + 2k₃ + k₄)
/// ```
///
/// # Accuracy
///
/// - Local truncation error: O(h⁵)
/// - Global error: O(h⁴)
///
/// This means halving the step size reduces the error by a factor of ~16.
///
/// # Stability
///
/// RK4 has a larger stability region than Euler's method, making it suitable
/// for a wider class of problems. For the test equation dy/dt = λy, the
/// stability region extends to approximately |hλ| < 2.78.
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
/// use axiomath::calculus::ode::runge_kutta_4;
///
/// // Solve dy/dt = y with y(0) = 1, exact solution: y(t) = e^t
/// let f = |_t: f64, y: f64| y;
/// let solution = runge_kutta_4(f, 1.0, (0.0, 1.0), 0.1);
///
/// let y_at_1 = solution.last().unwrap().1;
/// let exact = std::f64::consts::E;
/// assert!((y_at_1 - exact).abs() < 1e-6);
/// ```
///
/// ```
/// use axiomath::calculus::ode::runge_kutta_4;
///
/// // Solve dy/dt = -2ty with y(0) = 1, exact solution: y(t) = e^(-t²)
/// let f = |t: f64, y: f64| -2.0 * t * y;
/// let solution = runge_kutta_4(f, 1.0, (0.0, 1.0), 0.01);
///
/// let y_at_1 = solution.last().unwrap().1;
/// let exact = (-1.0_f64).exp();
/// assert!((y_at_1 - exact).abs() < 1e-8);
/// ```
///
/// # References
///
/// - Runge, Carl. "Über die numerische Auflösung von Differentialgleichungen." 1895.
/// - Kutta, Martin. "Beitrag zur näherungsweisen Integration totaler Differentialgleichungen." 1901.
/// - Burden & Faires, "Numerical Analysis", Section 5.4.
pub fn runge_kutta_4<T, F>(f: F, y0: T, t_range: (T, T), dt: T) -> Vec<(T, T)>
where
    T: Float,
    F: Fn(T, T) -> T,
{
    let (t0, t_final) = t_range;

    assert!(dt > T::zero(), "Step size must be positive");
    assert!(t_final >= t0, "Final time must be >= initial time");

    let two = T::one() + T::one();
    let three = two + T::one();
    let four = two * two;
    let six = three * two;

    // Calculate number of steps
    let n_steps = ((t_final - t0) / dt).ceil().to_usize().unwrap();

    // Initialize solution vector
    let mut solution = Vec::with_capacity(n_steps + 1);
    solution.push((t0, y0));

    let mut t = t0;
    let mut y = y0;

    for _ in 0..n_steps {
        // RK4 stages
        let k1 = f(t, y);
        let k2 = f(t + dt / two, y + dt * k1 / two);
        let k3 = f(t + dt / two, y + dt * k2 / two);
        let k4 = f(t + dt, y + dt * k3);

        // Weighted average: (k1 + 2k2 + 2k3 + k4) / 6
        let dy = (k1 + two * k2 + two * k3 + k4) / six;
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
    fn test_rk4_exponential_growth() {
        // dy/dt = y, y(0) = 1
        // Exact solution: y(t) = e^t
        let f = |_t: f64, y: f64| y;
        let solution = runge_kutta_4(f, 1.0, (0.0, 1.0), 0.1);

        let y_at_1 = solution.last().unwrap().1;
        let exact = std::f64::consts::E;

        // RK4 should be very accurate
        assert_abs_diff_eq!(y_at_1, exact, epsilon = 1e-5);
    }

    #[test]
    fn test_rk4_exponential_decay() {
        // dy/dt = -y, y(0) = 1
        // Exact solution: y(t) = e^(-t)
        let f = |_t: f64, y: f64| -y;
        let solution = runge_kutta_4(f, 1.0, (0.0, 2.0), 0.1);

        let y_at_2 = solution.last().unwrap().1;
        let exact = (-2.0_f64).exp();

        assert_abs_diff_eq!(y_at_2, exact, epsilon = 1e-6);
    }

    #[test]
    fn test_rk4_linear_ode() {
        // dy/dt = t, y(0) = 0
        // Exact solution: y(t) = t²/2
        let f = |t: f64, _y: f64| t;
        let solution = runge_kutta_4(f, 0.0, (0.0, 1.0), 0.1);

        let y_at_1 = solution.last().unwrap().1;
        let exact = 0.5;

        assert_abs_diff_eq!(y_at_1, exact, epsilon = 1e-10);
    }

    #[test]
    fn test_rk4_gaussian() {
        // dy/dt = -2ty, y(0) = 1
        // Exact solution: y(t) = e^(-t²)
        let f = |t: f64, y: f64| -2.0 * t * y;
        let solution = runge_kutta_4(f, 1.0, (0.0, 1.0), 0.01);

        let y_at_1 = solution.last().unwrap().1;
        let exact = (-1.0_f64).exp();

        assert_abs_diff_eq!(y_at_1, exact, epsilon = 1e-8);
    }

    #[test]
    fn test_rk4_sine_cosine() {
        // dy/dt = cos(t), y(0) = 0
        // Exact solution: y(t) = sin(t)
        let f = |t: f64, _y: f64| t.cos();
        let solution = runge_kutta_4(f, 0.0, (0.0, std::f64::consts::PI), 0.01);

        let y_at_pi = solution.last().unwrap().1;
        let exact = std::f64::consts::PI.sin(); // Should be ~0

        // At π, sin(π) is very close to 0 (within machine precision)
        assert_abs_diff_eq!(y_at_pi, exact, epsilon = 1e-2);
    }

    #[test]
    fn test_rk4_constant_solution() {
        // dy/dt = 0, y(0) = 7
        // Exact solution: y(t) = 7
        let f = |_t: f64, _y: f64| 0.0;
        let solution = runge_kutta_4(f, 7.0, (0.0, 1.0), 0.1);

        for &(_t, y) in &solution {
            assert_abs_diff_eq!(y, 7.0, epsilon = 1e-15);
        }
    }

    #[test]
    fn test_rk4_linear_growth() {
        // dy/dt = 3, y(0) = 2
        // Exact solution: y(t) = 2 + 3t
        let f = |_t: f64, _y: f64| 3.0;
        let solution = runge_kutta_4(f, 2.0, (0.0, 1.0), 0.1);

        let y_at_1 = solution.last().unwrap().1;
        let exact = 5.0;

        assert_abs_diff_eq!(y_at_1, exact, epsilon = 1e-14);
    }

    #[test]
    fn test_rk4_convergence() {
        // Verify fourth-order convergence: halving h should reduce error by ~16×
        let f = |_t: f64, y: f64| y;
        let exact = std::f64::consts::E;

        let sol_h1 = runge_kutta_4(f, 1.0, (0.0, 1.0), 0.2);
        let sol_h2 = runge_kutta_4(f, 1.0, (0.0, 1.0), 0.1);

        let error_h1 = (sol_h1.last().unwrap().1 - exact).abs();
        let error_h2 = (sol_h2.last().unwrap().1 - exact).abs();

        // For fourth-order method, error_h2 should be roughly error_h1/16
        let ratio = error_h1 / error_h2;
        assert!(ratio > 12.0 && ratio < 20.0);
    }

    #[test]
    fn test_rk4_better_than_euler() {
        // RK4 should be significantly more accurate than Euler for the same step size
        use crate::calculus::ode::euler_method;

        let f = |_t: f64, y: f64| y;
        let dt = 0.1;
        let exact = std::f64::consts::E;

        let rk4_sol = runge_kutta_4(f, 1.0, (0.0, 1.0), dt);
        let euler_sol = euler_method(f, 1.0, (0.0, 1.0), dt);

        let rk4_error = (rk4_sol.last().unwrap().1 - exact).abs();
        let euler_error = (euler_sol.last().unwrap().1 - exact).abs();

        assert!(rk4_error < euler_error / 100.0);
    }

    #[test]
    fn test_rk4_solution_length() {
        let f = |_t: f64, y: f64| y;
        let dt = 0.1;
        let solution = runge_kutta_4(f, 1.0, (0.0, 1.0), dt);

        // Should have approximately (1.0 - 0.0) / 0.1 + 1 = 11 points
        assert!(solution.len() >= 11);
        assert!(solution.len() <= 12);
    }

    #[test]
    fn test_rk4_time_values() {
        let f = |_t: f64, _y: f64| 0.0;
        let dt = 0.1;
        let solution = runge_kutta_4(f, 0.0, (0.0, 1.0), dt);

        // Check that time values are correct
        for (i, &(t, _y)) in solution.iter().enumerate() {
            let expected_t = (i as f64) * dt;
            assert_abs_diff_eq!(t, expected_t, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_rk4_zero_interval() {
        let f = |_t: f64, y: f64| y;
        let solution = runge_kutta_4(f, 1.0, (0.0, 0.0), 0.1);

        // Should return just the initial condition
        assert_eq!(solution.len(), 1);
        assert_eq!(solution[0], (0.0, 1.0));
    }

    #[test]
    fn test_rk4_high_accuracy() {
        // Test with very small step size for high accuracy
        let f = |_t: f64, y: f64| y;
        let solution = runge_kutta_4(f, 1.0, (0.0, 1.0), 0.001);

        let y_at_1 = solution.last().unwrap().1;
        let exact = std::f64::consts::E;

        // With dt=0.001, RK4 should be extremely accurate
        assert_abs_diff_eq!(y_at_1, exact, epsilon = 1e-12);
    }

    #[test]
    fn test_rk4_polynomial() {
        // dy/dt = 3t², y(0) = 0
        // Exact solution: y(t) = t³
        let f = |t: f64, _y: f64| 3.0 * t * t;
        let solution = runge_kutta_4(f, 0.0, (0.0, 2.0), 0.1);

        let y_at_2 = solution.last().unwrap().1;
        let exact = 8.0;

        assert_abs_diff_eq!(y_at_2, exact, epsilon = 1e-10);
    }

    #[test]
    fn test_rk4_f32_support() {
        // Test that method works with f32
        let f = |_t: f32, y: f32| y;
        let solution = runge_kutta_4(f, 1.0_f32, (0.0_f32, 1.0_f32), 0.1_f32);

        let y_at_1 = solution.last().unwrap().1;
        let exact = std::f32::consts::E;

        assert_abs_diff_eq!(y_at_1, exact, epsilon = 1e-5_f32);
    }

    #[test]
    #[should_panic(expected = "Step size must be positive")]
    fn test_rk4_negative_step() {
        let f = |_t: f64, y: f64| y;
        runge_kutta_4(f, 1.0, (0.0, 1.0), -0.1);
    }

    #[test]
    #[should_panic(expected = "Step size must be positive")]
    fn test_rk4_zero_step() {
        let f = |_t: f64, y: f64| y;
        runge_kutta_4(f, 1.0, (0.0, 1.0), 0.0);
    }

    #[test]
    #[should_panic(expected = "Final time must be >= initial time")]
    fn test_rk4_invalid_time_range() {
        let f = |_t: f64, y: f64| y;
        runge_kutta_4(f, 1.0, (1.0, 0.0), 0.1);
    }
}
