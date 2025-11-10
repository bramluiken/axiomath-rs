//! Adaptive Runge-Kutta methods with automatic step size control.
//!
//! These methods automatically adjust the step size based on error estimates,
//! providing efficient and accurate solutions.

use crate::core::Float;

/// Solves an ODE using adaptive Runge-Kutta-Fehlberg method (RK45).
///
/// # Algorithm
///
/// RK45 uses a 4th-order and 5th-order Runge-Kutta formula to estimate the
/// local truncation error. The step size is then adjusted to keep the error
/// below a specified tolerance.
///
/// The method computes six function evaluations per step:
/// ```text
/// k₁ = f(tₙ, yₙ)
/// k₂ = f(tₙ + h/4, yₙ + h·k₁/4)
/// k₃ = f(tₙ + 3h/8, yₙ + h·(3k₁ + 9k₂)/32)
/// k₄ = f(tₙ + 12h/13, yₙ + h·(1932k₁ - 7200k₂ + 7296k₃)/2197)
/// k₅ = f(tₙ + h, yₙ + h·(439k₁/216 - 8k₂ + 3680k₃/513 - 845k₄/4104))
/// k₆ = f(tₙ + h/2, yₙ + h·(-8k₁/27 + 2k₂ - 3544k₃/2565 + 1859k₄/4104 - 11k₅/40))
/// ```
///
/// 4th-order: y₄ = yₙ + h·(25k₁/216 + 1408k₃/2565 + 2197k₄/4104 - k₅/5)
/// 5th-order: y₅ = yₙ + h·(16k₁/135 + 6656k₃/12825 + 28561k₄/56430 - 9k₅/50 + 2k₆/55)
///
/// Error estimate: E = |y₅ - y₄|
///
/// # Advantages
///
/// - Automatic step size control for efficiency
/// - Maintains user-specified accuracy
/// - Suitable for problems with varying time scales
/// - More efficient than fixed-step methods for many problems
///
/// # Arguments
///
/// * `f` - The derivative function f(t, y) where dy/dt = f(t, y)
/// * `y0` - Initial condition y(t₀) = y₀
/// * `t_range` - Time interval (t₀, t_final)
/// * `tolerance` - Desired absolute error tolerance per step
///
/// # Returns
///
/// A vector of (time, value) pairs. Note that the time steps are not uniform
/// and are chosen adaptively based on the error estimate.
///
/// # Panics
///
/// Panics if `tolerance` is not positive or if `t_range` is invalid.
///
/// # Examples
///
/// ```
/// use axiomath::calculus::ode::adaptive_rk45;
///
/// // Solve dy/dt = y with y(0) = 1
/// let f = |_t: f64, y: f64| y;
/// let solution = adaptive_rk45(f, 1.0, (0.0, 1.0), 1e-6);
///
/// let y_at_1 = solution.last().unwrap().1;
/// let exact = std::f64::consts::E;
/// assert!((y_at_1 - exact).abs() < 1e-5);
/// ```
///
/// ```
/// use axiomath::calculus::ode::adaptive_rk45;
///
/// // Solve a stiff problem: dy/dt = -50(y - cos(t))
/// let f = |t: f64, y: f64| -50.0 * (y - t.cos());
/// let solution = adaptive_rk45(f, 1.0, (0.0, 1.0), 1e-4);
///
/// // Adaptive method automatically uses smaller steps where needed
/// assert!(solution.len() > 10);
/// ```
///
/// # References
///
/// - Fehlberg, E. "Low-order classical Runge-Kutta formulas with step size control." 1969.
/// - Burden & Faires, "Numerical Analysis", Section 5.5.
pub fn adaptive_rk45<T, F>(f: F, y0: T, t_range: (T, T), tolerance: T) -> Vec<(T, T)>
where
    T: Float,
    F: Fn(T, T) -> T,
{
    let (t0, t_final) = t_range;

    assert!(tolerance > T::zero(), "Tolerance must be positive");
    assert!(t_final >= t0, "Final time must be >= initial time");

    // Safety factor for step size adjustment
    let safety = T::from(0.9).unwrap();
    let min_scale = T::from(0.2).unwrap();
    let max_scale = T::from(5.0).unwrap();

    // Initial step size guess
    let mut h = (t_final - t0) / T::from(100).unwrap();
    let h_min = (t_final - t0) / T::from(1000000).unwrap();

    let mut solution = Vec::new();
    solution.push((t0, y0));

    let mut t = t0;
    let mut y = y0;

    while t < t_final {
        // Ensure we don't step past t_final
        if t + h > t_final {
            h = t_final - t;
        }

        // Compute RK45 step
        let (y_new, error) = rk45_step(&f, t, y, h);

        // Compute optimal step size
        let scale = if error > T::zero() {
            safety * (tolerance / error).powf(T::from(0.2).unwrap())
        } else {
            max_scale
        };

        let scale = scale.max(min_scale).min(max_scale);

        if error <= tolerance || h <= h_min {
            // Accept the step
            t = t + h;
            y = y_new;
            solution.push((t, y));

            // Adjust step size for next iteration
            h = h * scale;
        } else {
            // Reject the step and retry with smaller h
            h = h * scale;
            if h < h_min {
                h = h_min;
            }
        }
    }

    solution
}

/// Performs a single RK45 step and returns the new value and error estimate.
fn rk45_step<T, F>(f: &F, t: T, y: T, h: T) -> (T, T)
where
    T: Float,
    F: Fn(T, T) -> T,
{
    // Butcher tableau coefficients for RK45
    let two = T::one() + T::one();
    let three = two + T::one();
    let four = three + T::one();
    let five = four + T::one();

    // RK45 stages
    let k1 = f(t, y);

    let k2 = f(t + h / four, y + h * k1 / four);

    let k3 = f(
        t + three * h / (four * two),
        y + h * (three * k1 + T::from(9).unwrap() * k2) / T::from(32).unwrap(),
    );

    let k4 = f(
        t + T::from(12).unwrap() * h / T::from(13).unwrap(),
        y + h
            * (T::from(1932).unwrap() * k1
                - T::from(7200).unwrap() * k2
                + T::from(7296).unwrap() * k3)
            / T::from(2197).unwrap(),
    );

    let k5 = f(
        t + h,
        y + h
            * (T::from(439).unwrap() * k1 / T::from(216).unwrap()
                - T::from(8).unwrap() * k2
                + T::from(3680).unwrap() * k3 / T::from(513).unwrap()
                - T::from(845).unwrap() * k4 / T::from(4104).unwrap()),
    );

    let k6 = f(
        t + h / two,
        y + h
            * (-T::from(8).unwrap() * k1 / T::from(27).unwrap()
                + two * k2
                - T::from(3544).unwrap() * k3 / T::from(2565).unwrap()
                + T::from(1859).unwrap() * k4 / T::from(4104).unwrap()
                - T::from(11).unwrap() * k5 / T::from(40).unwrap()),
    );

    // 4th-order solution
    let y4 = y
        + h * (T::from(25).unwrap() * k1 / T::from(216).unwrap()
            + T::from(1408).unwrap() * k3 / T::from(2565).unwrap()
            + T::from(2197).unwrap() * k4 / T::from(4104).unwrap()
            - k5 / five);

    // 5th-order solution
    let y5 = y
        + h * (T::from(16).unwrap() * k1 / T::from(135).unwrap()
            + T::from(6656).unwrap() * k3 / T::from(12825).unwrap()
            + T::from(28561).unwrap() * k4 / T::from(56430).unwrap()
            - T::from(9).unwrap() * k5 / T::from(50).unwrap()
            + two * k6 / T::from(55).unwrap());

    // Error estimate
    let error = (y5 - y4).abs();

    (y5, error)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_adaptive_rk45_exponential_growth() {
        // dy/dt = y, y(0) = 1
        // Exact solution: y(t) = e^t
        let f = |_t: f64, y: f64| y;
        let solution = adaptive_rk45(f, 1.0, (0.0, 1.0), 1e-8);

        let y_at_1 = solution.last().unwrap().1;
        let exact = std::f64::consts::E;

        assert_abs_diff_eq!(y_at_1, exact, epsilon = 1e-6);
    }

    #[test]
    fn test_adaptive_rk45_exponential_decay() {
        // dy/dt = -y, y(0) = 1
        // Exact solution: y(t) = e^(-t)
        let f = |_t: f64, y: f64| -y;
        let solution = adaptive_rk45(f, 1.0, (0.0, 2.0), 1e-8);

        let y_at_2 = solution.last().unwrap().1;
        let exact = (-2.0_f64).exp();

        assert_abs_diff_eq!(y_at_2, exact, epsilon = 1e-6);
    }

    #[test]
    fn test_adaptive_rk45_linear() {
        // dy/dt = t, y(0) = 0
        // Exact solution: y(t) = t²/2
        let f = |t: f64, _y: f64| t;
        let solution = adaptive_rk45(f, 0.0, (0.0, 1.0), 1e-8);

        let y_at_1 = solution.last().unwrap().1;
        let exact = 0.5;

        assert_abs_diff_eq!(y_at_1, exact, epsilon = 1e-10);
    }

    #[test]
    fn test_adaptive_rk45_gaussian() {
        // dy/dt = -2ty, y(0) = 1
        // Exact solution: y(t) = e^(-t²)
        let f = |t: f64, y: f64| -2.0 * t * y;
        let solution = adaptive_rk45(f, 1.0, (0.0, 1.0), 1e-8);

        let y_at_1 = solution.last().unwrap().1;
        let exact = (-1.0_f64).exp();

        assert_abs_diff_eq!(y_at_1, exact, epsilon = 1e-8);
    }

    #[test]
    fn test_adaptive_rk45_adaptive_steps() {
        // Test that adaptive method uses variable step sizes
        let f = |_t: f64, y: f64| y;
        let solution = adaptive_rk45(f, 1.0, (0.0, 1.0), 1e-6);

        // Compute step sizes
        let mut step_sizes = Vec::new();
        for i in 1..solution.len() {
            let dt = solution[i].0 - solution[i - 1].0;
            step_sizes.push(dt);
        }

        // Not all steps should be the same size
        let first_step = step_sizes[0];
        let all_same = step_sizes.iter().all(|&dt| (dt - first_step).abs() < 1e-10);
        assert!(!all_same, "Adaptive method should vary step sizes");
    }

    #[test]
    fn test_adaptive_rk45_tolerance_effect() {
        // Tighter tolerance should result in more steps
        let f = |_t: f64, y: f64| y;

        let sol_loose = adaptive_rk45(f, 1.0, (0.0, 1.0), 1e-4);
        let sol_tight = adaptive_rk45(f, 1.0, (0.0, 1.0), 1e-8);

        assert!(sol_tight.len() > sol_loose.len());
    }

    #[test]
    fn test_adaptive_rk45_constant() {
        // dy/dt = 0, y(0) = 5
        // Exact solution: y(t) = 5
        let f = |_t: f64, _y: f64| 0.0;
        let solution = adaptive_rk45(f, 5.0, (0.0, 1.0), 1e-8);

        for &(_t, y) in &solution {
            assert_abs_diff_eq!(y, 5.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn test_adaptive_rk45_sine() {
        // dy/dt = cos(t), y(0) = 0
        // Exact solution: y(t) = sin(t)
        let f = |t: f64, _y: f64| t.cos();
        let solution = adaptive_rk45(f, 0.0, (0.0, std::f64::consts::PI), 1e-8);

        let y_at_pi = solution.last().unwrap().1;
        let exact = std::f64::consts::PI.sin(); // Should be ~0

        // At π, sin(π) is very close to 0 (within machine precision)
        assert_abs_diff_eq!(y_at_pi, exact, epsilon = 1e-2);
    }

    #[test]
    fn test_adaptive_rk45_high_accuracy() {
        // Test with very tight tolerance
        let f = |_t: f64, y: f64| y;
        let solution = adaptive_rk45(f, 1.0, (0.0, 1.0), 1e-12);

        let y_at_1 = solution.last().unwrap().1;
        let exact = std::f64::consts::E;

        assert_abs_diff_eq!(y_at_1, exact, epsilon = 1e-10);
    }

    #[test]
    fn test_adaptive_rk45_rapid_change() {
        // Test function with rapid changes: dy/dt = 10y
        let f = |_t: f64, y: f64| 10.0 * y;
        let solution = adaptive_rk45(f, 1.0, (0.0, 0.5), 1e-6);

        let y_at_05 = solution.last().unwrap().1;
        let exact = (5.0_f64).exp();

        assert_abs_diff_eq!(y_at_05, exact, epsilon = 1e-4);
    }

    #[test]
    fn test_adaptive_rk45_zero_interval() {
        let f = |_t: f64, y: f64| y;
        let solution = adaptive_rk45(f, 1.0, (0.0, 0.0), 1e-6);

        // Should return just the initial condition
        assert_eq!(solution.len(), 1);
        assert_eq!(solution[0], (0.0, 1.0));
    }

    #[test]
    fn test_adaptive_rk45_f32_support() {
        // Test that method works with f32
        let f = |_t: f32, y: f32| y;
        let solution = adaptive_rk45(f, 1.0_f32, (0.0_f32, 1.0_f32), 1e-4_f32);

        let y_at_1 = solution.last().unwrap().1;
        let exact = std::f32::consts::E;

        assert_abs_diff_eq!(y_at_1, exact, epsilon = 1e-3_f32);
    }

    #[test]
    #[should_panic(expected = "Tolerance must be positive")]
    fn test_adaptive_rk45_negative_tolerance() {
        let f = |_t: f64, y: f64| y;
        adaptive_rk45(f, 1.0, (0.0, 1.0), -1e-6);
    }

    #[test]
    #[should_panic(expected = "Tolerance must be positive")]
    fn test_adaptive_rk45_zero_tolerance() {
        let f = |_t: f64, y: f64| y;
        adaptive_rk45(f, 1.0, (0.0, 1.0), 0.0);
    }

    #[test]
    #[should_panic(expected = "Final time must be >= initial time")]
    fn test_adaptive_rk45_invalid_time_range() {
        let f = |_t: f64, y: f64| y;
        adaptive_rk45(f, 1.0, (1.0, 0.0), 1e-6);
    }
}
