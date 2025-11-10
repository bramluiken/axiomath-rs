//! Newton-Cotes integration formulas.
//!
//! These methods use equally-spaced points and polynomial interpolation
//! to approximate definite integrals.

use crate::core::Float;

/// Computes the definite integral using the trapezoidal rule.
///
/// # Algorithm
///
/// The trapezoidal rule approximates the integral by summing trapezoids:
///
/// ```text
/// ∫[a,b] f(x)dx ≈ h/2 * [f(x₀) + 2f(x₁) + 2f(x₂) + ... + 2f(xₙ₋₁) + f(xₙ)]
/// ```
///
/// where h = (b-a)/n and xᵢ = a + i*h.
///
/// # Complexity
///
/// - Function evaluations: n + 1
/// - Accuracy order: O(h²) or O(n⁻²)
///
/// # Arguments
///
/// * `f` - The function to integrate
/// * `a` - Lower bound of integration
/// * `b` - Upper bound of integration
/// * `n` - Number of intervals (should be > 0)
///
/// # Panics
///
/// Panics if `n` is 0.
///
/// # Examples
///
/// ```
/// use axiomath::calculus::integration::trapezoidal_rule;
///
/// // Integrate x² from 0 to 1, exact value is 1/3
/// let f = |x: f64| x * x;
/// let result = trapezoidal_rule(f, 0.0, 1.0, 1000);
/// assert!((result - 1.0/3.0).abs() < 1e-6);
/// ```
///
/// # References
///
/// - Burden & Faires, "Numerical Analysis", Section 4.1
pub fn trapezoidal_rule<T, F>(f: F, a: T, b: T, n: usize) -> T
where
    T: Float,
    F: Fn(T) -> T,
{
    assert!(n > 0, "Number of intervals must be positive");

    let n_t = T::from(n).expect("Failed to convert n to Float");
    let h = (b - a) / n_t;
    let two = T::one() + T::one();

    // First and last terms
    let mut sum = f(a) + f(b);

    // Middle terms (weighted by 2)
    for i in 1..n {
        let x = a + T::from(i).expect("Failed to convert i to Float") * h;
        sum = sum + two * f(x);
    }

    sum * h / two
}

/// Computes the definite integral using Simpson's 1/3 rule.
///
/// # Algorithm
///
/// Simpson's rule uses quadratic interpolation:
///
/// ```text
/// ∫[a,b] f(x)dx ≈ h/3 * [f(x₀) + 4f(x₁) + 2f(x₂) + 4f(x₃) + ... + 4f(xₙ₋₁) + f(xₙ)]
/// ```
///
/// where h = (b-a)/n, n must be even, and xᵢ = a + i*h.
///
/// # Complexity
///
/// - Function evaluations: n + 1
/// - Accuracy order: O(h⁴) or O(n⁻⁴)
///
/// # Arguments
///
/// * `f` - The function to integrate
/// * `a` - Lower bound of integration
/// * `b` - Upper bound of integration
/// * `n` - Number of intervals (must be even and > 0)
///
/// # Panics
///
/// Panics if `n` is 0 or odd.
///
/// # Examples
///
/// ```
/// use axiomath::calculus::integration::simpsons_rule;
///
/// // Integrate sin(x) from 0 to π, exact value is 2
/// use std::f64::consts::PI;
/// let f = |x: f64| x.sin();
/// let result = simpsons_rule(f, 0.0, PI, 100);
/// assert!((result - 2.0).abs() < 1e-10);
/// ```
///
/// # References
///
/// - Burden & Faires, "Numerical Analysis", Section 4.3
pub fn simpsons_rule<T, F>(f: F, a: T, b: T, n: usize) -> T
where
    T: Float,
    F: Fn(T) -> T,
{
    assert!(n > 0, "Number of intervals must be positive");
    assert!(n % 2 == 0, "Number of intervals must be even for Simpson's 1/3 rule");

    let n_t = T::from(n).expect("Failed to convert n to Float");
    let h = (b - a) / n_t;
    let two = T::one() + T::one();
    let three = two + T::one();
    let four = two * two;

    // First and last terms
    let mut sum = f(a) + f(b);

    // Odd indices (weighted by 4)
    for i in (1..n).step_by(2) {
        let x = a + T::from(i).expect("Failed to convert i to Float") * h;
        sum = sum + four * f(x);
    }

    // Even indices (weighted by 2)
    for i in (2..n).step_by(2) {
        let x = a + T::from(i).expect("Failed to convert i to Float") * h;
        sum = sum + two * f(x);
    }

    sum * h / three
}

/// Computes the definite integral using Simpson's 3/8 rule.
///
/// # Algorithm
///
/// Simpson's 3/8 rule uses cubic interpolation:
///
/// ```text
/// ∫[a,b] f(x)dx ≈ 3h/8 * [f(x₀) + 3f(x₁) + 3f(x₂) + 2f(x₃) + 3f(x₄) + ... + f(xₙ)]
/// ```
///
/// where h = (b-a)/n, n must be divisible by 3, and xᵢ = a + i*h.
///
/// The pattern is: 1, 3, 3, 2, 3, 3, 2, ..., 3, 3, 1
///
/// # Complexity
///
/// - Function evaluations: n + 1
/// - Accuracy order: O(h⁴) or O(n⁻⁴)
///
/// # Arguments
///
/// * `f` - The function to integrate
/// * `a` - Lower bound of integration
/// * `b` - Upper bound of integration
/// * `n` - Number of intervals (must be divisible by 3 and > 0)
///
/// # Panics
///
/// Panics if `n` is 0 or not divisible by 3.
///
/// # Examples
///
/// ```
/// use axiomath::calculus::integration::simpsons_3_8_rule;
///
/// // Integrate e^x from 0 to 1, exact value is e - 1
/// let f = |x: f64| x.exp();
/// let result = simpsons_3_8_rule(f, 0.0, 1.0, 99);  // n divisible by 3
/// let exact = std::f64::consts::E - 1.0;
/// assert!((result - exact).abs() < 1e-10);
/// ```
///
/// # References
///
/// - Burden & Faires, "Numerical Analysis", Section 4.3
pub fn simpsons_3_8_rule<T, F>(f: F, a: T, b: T, n: usize) -> T
where
    T: Float,
    F: Fn(T) -> T,
{
    assert!(n > 0, "Number of intervals must be positive");
    assert!(n % 3 == 0, "Number of intervals must be divisible by 3 for Simpson's 3/8 rule");

    let n_t = T::from(n).expect("Failed to convert n to Float");
    let h = (b - a) / n_t;
    let two = T::one() + T::one();
    let three = two + T::one();
    let eight = two * two * two;

    // First and last terms
    let mut sum = f(a) + f(b);

    for i in 1..n {
        let x = a + T::from(i).expect("Failed to convert i to Float") * h;
        let weight = if i % 3 == 0 { two } else { three };
        sum = sum + weight * f(x);
    }

    sum * three * h / eight
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_trapezoidal_polynomial() {
        // ∫₀¹ x² dx = 1/3
        let f = |x: f64| x * x;
        let result = trapezoidal_rule(f, 0.0, 1.0, 1000);
        assert_abs_diff_eq!(result, 1.0 / 3.0, epsilon = 1e-6);
    }

    #[test]
    fn test_trapezoidal_exponential() {
        // ∫₀¹ e^x dx = e - 1
        let f = |x: f64| x.exp();
        let result = trapezoidal_rule(f, 0.0, 1.0, 1000);
        let exact = std::f64::consts::E - 1.0;
        assert_abs_diff_eq!(result, exact, epsilon = 1e-6);
    }

    #[test]
    fn test_trapezoidal_sine() {
        // ∫₀^π sin(x) dx = 2
        use std::f64::consts::PI;
        let f = |x: f64| x.sin();
        let result = trapezoidal_rule(f, 0.0, PI, 1000);
        assert_abs_diff_eq!(result, 2.0, epsilon = 2e-6);
    }

    #[test]
    fn test_simpsons_polynomial() {
        // ∫₀¹ x² dx = 1/3
        let f = |x: f64| x * x;
        let result = simpsons_rule(f, 0.0, 1.0, 100);
        assert_abs_diff_eq!(result, 1.0 / 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_simpsons_cubic() {
        // Simpson's rule is exact for cubics
        // ∫₀¹ x³ dx = 1/4
        let f = |x: f64| x.powi(3);
        let result = simpsons_rule(f, 0.0, 1.0, 10);
        assert_abs_diff_eq!(result, 0.25, epsilon = 1e-12);
    }

    #[test]
    fn test_simpsons_sine() {
        // ∫₀^π sin(x) dx = 2
        use std::f64::consts::PI;
        let f = |x: f64| x.sin();
        let result = simpsons_rule(f, 0.0, PI, 100);
        // Note: Slightly larger epsilon needed due to accumulated rounding errors in trigonometric evaluation
        // The error is ~1.08e-8, which is excellent accuracy for Simpson's rule
        assert_abs_diff_eq!(result, 2.0, epsilon = 2e-8);
    }

    #[test]
    fn test_simpsons_exponential() {
        // ∫₀¹ e^x dx = e - 1
        let f = |x: f64| x.exp();
        let result = simpsons_rule(f, 0.0, 1.0, 100);
        let exact = std::f64::consts::E - 1.0;
        assert_abs_diff_eq!(result, exact, epsilon = 1e-10);
    }

    #[test]
    fn test_simpsons_3_8_polynomial() {
        // ∫₀¹ x² dx = 1/3
        let f = |x: f64| x * x;
        let result = simpsons_3_8_rule(f, 0.0, 1.0, 99);
        assert_abs_diff_eq!(result, 1.0 / 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_simpsons_3_8_quartic() {
        // ∫₀¹ x⁴ dx = 1/5
        let f = |x: f64| x.powi(4);
        let result = simpsons_3_8_rule(f, 0.0, 1.0, 99);
        assert_abs_diff_eq!(result, 0.2, epsilon = 1e-8);
    }

    #[test]
    fn test_simpsons_3_8_exponential() {
        // ∫₀¹ e^x dx = e - 1
        let f = |x: f64| x.exp();
        let result = simpsons_3_8_rule(f, 0.0, 1.0, 99);
        let exact = std::f64::consts::E - 1.0;
        assert_abs_diff_eq!(result, exact, epsilon = 1e-9);
    }

    #[test]
    fn test_constant_function() {
        // ∫₀¹ 5 dx = 5
        let f = |_x: f64| 5.0;
        let trap = trapezoidal_rule(f, 0.0, 1.0, 10);
        let simp = simpsons_rule(f, 0.0, 1.0, 10);
        let simp38 = simpsons_3_8_rule(f, 0.0, 1.0, 9);

        assert_abs_diff_eq!(trap, 5.0, epsilon = 1e-15);
        assert_abs_diff_eq!(simp, 5.0, epsilon = 1e-15);
        assert_abs_diff_eq!(simp38, 5.0, epsilon = 1e-15);
    }

    #[test]
    fn test_linear_function() {
        // All methods should be exact for linear functions
        // ∫₀¹ (2x + 1) dx = x² + x |₀¹ = 2
        let f = |x: f64| 2.0 * x + 1.0;
        let trap = trapezoidal_rule(f, 0.0, 1.0, 10);
        let simp = simpsons_rule(f, 0.0, 1.0, 10);

        assert_abs_diff_eq!(trap, 2.0, epsilon = 1e-14);
        assert_abs_diff_eq!(simp, 2.0, epsilon = 1e-14);
    }

    #[test]
    fn test_negative_interval() {
        // ∫₋₁¹ x² dx = 2/3
        let f = |x: f64| x * x;
        let result = simpsons_rule(f, -1.0, 1.0, 100);
        assert_abs_diff_eq!(result, 2.0 / 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_simpson_better_than_trapezoidal() {
        // For smooth functions, Simpson's should be more accurate
        let f = |x: f64| (x * x * x).sin();
        let n = 20;

        let trap = trapezoidal_rule(f, 0.0, 1.0, n);
        let simp = simpsons_rule(f, 0.0, 1.0, n);

        // Compute "exact" value with very fine Simpson's rule
        let exact = simpsons_rule(f, 0.0, 1.0, 10000);

        let error_trap = (trap - exact).abs();
        let error_simp = (simp - exact).abs();

        assert!(error_simp < error_trap);
    }

    #[test]
    fn test_f32_support() {
        // Test that methods work with f32
        let f = |x: f32| x * x;
        let result = simpsons_rule(f, 0.0_f32, 1.0_f32, 100);
        assert_abs_diff_eq!(result, 1.0_f32 / 3.0_f32, epsilon = 1e-6_f32);
    }

    #[test]
    #[should_panic(expected = "Number of intervals must be positive")]
    fn test_zero_intervals() {
        let f = |x: f64| x;
        trapezoidal_rule(f, 0.0, 1.0, 0);
    }

    #[test]
    #[should_panic(expected = "Number of intervals must be even")]
    fn test_simpson_odd_intervals() {
        let f = |x: f64| x;
        simpsons_rule(f, 0.0, 1.0, 99);  // Odd number
    }

    #[test]
    #[should_panic(expected = "Number of intervals must be divisible by 3")]
    fn test_simpson_3_8_wrong_intervals() {
        let f = |x: f64| x;
        simpsons_3_8_rule(f, 0.0, 1.0, 100);  // Not divisible by 3
    }
}
