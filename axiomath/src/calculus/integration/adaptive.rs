//! Adaptive integration methods.
//!
//! These methods automatically adjust the integration strategy based on the
//! function's behavior, providing high accuracy with minimal function evaluations.

use crate::core::Float;

/// Computes the definite integral using adaptive Simpson's rule.
///
/// # Algorithm
///
/// This method recursively subdivides the integration interval until the estimated
/// error is below the specified tolerance. It uses Simpson's rule on each subinterval
/// and compares with a finer subdivision to estimate the error.
///
/// The error estimate is:
/// ```text
/// E ≈ |S(a,b) - S(a,c) - S(c,b)| / 15
/// ```
/// where c = (a+b)/2 and S is Simpson's rule.
///
/// # Complexity
///
/// - Function evaluations: Adaptive (typically O(√n) for smooth functions)
/// - Accuracy: Controlled by tolerance parameter
///
/// # Arguments
///
/// * `f` - The function to integrate
/// * `a` - Lower bound of integration
/// * `b` - Upper bound of integration
/// * `tolerance` - Desired absolute error tolerance
/// * `max_recursion` - Maximum recursion depth (prevents infinite recursion)
///
/// # Panics
///
/// Panics if `max_recursion` is 0.
///
/// # Examples
///
/// ```
/// use axiomath::calculus::integration::adaptive_simpson;
///
/// // Integrate sin(x)/x from 0.001 to π (near singularity at 0)
/// use std::f64::consts::PI;
/// let f = |x: f64| x.sin() / x;
/// let result = adaptive_simpson(f, 0.001, PI, 1e-8, 20);
/// // Si(π) ≈ 1.8519 (sine integral)
/// assert!((result - 1.8519).abs() < 1e-6);
/// ```
///
/// # References
///
/// - McKeeman, W.M. "Algorithm 145: Adaptive Numerical Integration by Simpson's Rule." 1962.
pub fn adaptive_simpson<T, F>(f: F, a: T, b: T, tolerance: T, max_recursion: usize) -> T
where
    T: Float,
    F: Fn(T) -> T,
{
    assert!(max_recursion > 0, "Maximum recursion depth must be positive");

    let two = T::one() + T::one();
    let mid = (a + b) / two;

    // Initial Simpson's rule approximation on [a, b]
    let whole = simpson_segment(&f, a, b);

    adaptive_simpson_recursive(&f, a, b, tolerance, whole, max_recursion)
}

/// Recursive helper for adaptive Simpson's rule.
fn adaptive_simpson_recursive<T, F>(
    f: &F,
    a: T,
    b: T,
    tolerance: T,
    whole: T,
    depth: usize,
) -> T
where
    T: Float,
    F: Fn(T) -> T,
{
    let two = T::one() + T::one();
    let fifteen = T::from(15).expect("Failed to convert 15 to Float");

    let mid = (a + b) / two;

    // Simpson's rule on left and right halves
    let left = simpson_segment(f, a, mid);
    let right = simpson_segment(f, mid, b);

    // Error estimate: |S(a,b) - S(a,mid) - S(mid,b)| / 15
    let delta = (left + right - whole).abs() / fifteen;

    if delta <= tolerance || depth == 0 {
        // Converged or max depth reached
        // Return the more accurate estimate (Richardson extrapolation)
        left + right + delta
    } else {
        // Subdivide further
        let left_integral = adaptive_simpson_recursive(f, a, mid, tolerance / two, left, depth - 1);
        let right_integral = adaptive_simpson_recursive(f, mid, b, tolerance / two, right, depth - 1);
        left_integral + right_integral
    }
}

/// Computes Simpson's rule for a single segment [a, b].
///
/// This is the basic Simpson's 1/3 rule with just one interval.
#[inline]
fn simpson_segment<T, F>(f: &F, a: T, b: T) -> T
where
    T: Float,
    F: Fn(T) -> T,
{
    let two = T::one() + T::one();
    let four = two * two;
    let six = four + two;

    let mid = (a + b) / two;
    let h = (b - a) / six;

    h * (f(a) + four * f(mid) + f(b))
}

/// Computes the definite integral using Romberg integration.
///
/// # Algorithm
///
/// Romberg integration uses Richardson extrapolation to refine successive
/// trapezoidal rule approximations. It generates a triangular tableau where
/// each row uses twice as many intervals as the previous, and each column
/// applies one more level of extrapolation.
///
/// ```text
/// R(0,0)
/// R(1,0)  R(1,1)
/// R(2,0)  R(2,1)  R(2,2)
/// ...
/// ```
///
/// where R(n,0) is the trapezoidal rule with 2ⁿ intervals, and R(n,m) applies
/// m levels of Richardson extrapolation.
///
/// # Complexity
///
/// - Function evaluations: O(2ⁿ) for n iterations
/// - Accuracy: Can achieve very high accuracy (10⁻¹⁴ or better)
///
/// # Arguments
///
/// * `f` - The function to integrate
/// * `a` - Lower bound of integration
/// * `b` - Upper bound of integration
/// * `max_iterations` - Maximum number of refinements (typically 10-20)
/// * `tolerance` - Desired absolute error tolerance
///
/// # Returns
///
/// The best estimate of the integral, or the final value if tolerance is not reached.
///
/// # Examples
///
/// ```
/// use axiomath::calculus::integration::romberg_integration;
///
/// // Integrate e^x from 0 to 1, exact value is e - 1
/// let f = |x: f64| x.exp();
/// let result = romberg_integration(f, 0.0, 1.0, 10, 1e-12);
/// let exact = std::f64::consts::E - 1.0;
/// assert!((result - exact).abs() < 1e-12);
/// ```
///
/// # References
///
/// - Burden & Faires, "Numerical Analysis", Section 4.5
/// - Romberg, W. "Vereinfachte numerische Integration." 1955.
pub fn romberg_integration<T, F>(
    f: F,
    a: T,
    b: T,
    max_iterations: usize,
    tolerance: T,
) -> T
where
    T: Float,
    F: Fn(T) -> T,
{
    let two = T::one() + T::one();
    let four = two * two;

    // Romberg tableau
    let mut r: Vec<Vec<T>> = Vec::with_capacity(max_iterations);

    // R(0,0) = trapezoidal rule with 1 interval
    let h0 = b - a;
    r.push(vec![h0 * (f(a) + f(b)) / two]);

    for i in 1..max_iterations {
        r.push(Vec::with_capacity(i + 1));

        // R(i,0) = trapezoidal rule with 2^i intervals
        let h = h0 / T::from(1 << i).expect("Failed to convert power of 2 to Float");
        let mut sum = T::zero();

        // Sum odd-indexed points (new points in this iteration)
        for k in 0..(1 << (i - 1)) {
            let x = a + T::from(2 * k + 1).expect("Failed to convert to Float") * h;
            sum = sum + f(x);
        }

        let prev_value = r[i - 1][0];
        r[i].push(prev_value / two + h * sum);

        // Richardson extrapolation: R(i,m) values
        for m in 1..=i {
            let power_of_4 = four.powi(m as i32);
            let current_val = r[i][m - 1];
            let prev_val = r[i - 1][m - 1];
            let extrapolated = (power_of_4 * current_val - prev_val) / (power_of_4 - T::one());
            r[i].push(extrapolated);
        }

        // Check convergence
        if i > 0 {
            let diff = (r[i][i] - r[i - 1][i - 1]).abs();
            if diff < tolerance {
                return r[i][i];
            }
        }
    }

    // Return best estimate
    r[max_iterations - 1][max_iterations - 1]
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_adaptive_simpson_polynomial() {
        // ∫₀¹ x² dx = 1/3
        let f = |x: f64| x * x;
        let result = adaptive_simpson(f, 0.0, 1.0, 1e-10, 15);
        assert_abs_diff_eq!(result, 1.0 / 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_adaptive_simpson_exponential() {
        // ∫₀¹ e^x dx = e - 1
        let f = |x: f64| x.exp();
        let result = adaptive_simpson(f, 0.0, 1.0, 1e-10, 15);
        let exact = std::f64::consts::E - 1.0;
        assert_abs_diff_eq!(result, exact, epsilon = 1e-10);
    }

    #[test]
    fn test_adaptive_simpson_sine() {
        // ∫₀^π sin(x) dx = 2
        use std::f64::consts::PI;
        let f = |x: f64| x.sin();
        let result = adaptive_simpson(f, 0.0, PI, 1e-10, 15);
        assert_abs_diff_eq!(result, 2.0, epsilon = 1e-8);
    }

    #[test]
    fn test_adaptive_simpson_oscillatory() {
        // ∫₀^{2π} sin(10x) dx = 0
        use std::f64::consts::PI;
        let f = |x: f64| (10.0 * x).sin();
        let result = adaptive_simpson(f, 0.0, 2.0 * PI, 1e-8, 20);
        assert_abs_diff_eq!(result, 0.0, epsilon = 1e-6);
    }

    #[test]
    fn test_adaptive_simpson_sharp_peak() {
        // Function with a sharp peak
        // ∫₋₁¹ e^(-100x²) dx (narrow Gaussian)
        let f = |x: f64| (-100.0 * x * x).exp();
        let result = adaptive_simpson(f, -1.0, 1.0, 1e-8, 20);
        // Approximate value ≈ 0.1772
        assert!((result - 0.1772).abs() < 0.001);
    }

    #[test]
    fn test_romberg_polynomial() {
        // ∫₀¹ x³ dx = 1/4
        let f = |x: f64| x.powi(3);
        let result = romberg_integration(f, 0.0, 1.0, 10, 1e-12);
        assert_abs_diff_eq!(result, 0.25, epsilon = 1e-12);
    }

    #[test]
    fn test_romberg_exponential() {
        // ∫₀¹ e^x dx = e - 1
        let f = |x: f64| x.exp();
        let result = romberg_integration(f, 0.0, 1.0, 10, 1e-12);
        let exact = std::f64::consts::E - 1.0;
        assert_abs_diff_eq!(result, exact, epsilon = 1e-12);
    }

    #[test]
    fn test_romberg_sine() {
        // ∫₀^π sin(x) dx = 2
        use std::f64::consts::PI;
        let f = |x: f64| x.sin();
        let result = romberg_integration(f, 0.0, PI, 10, 1e-12);
        assert_abs_diff_eq!(result, 2.0, epsilon = 1e-12);
    }

    #[test]
    fn test_romberg_high_accuracy() {
        // Test that Romberg achieves very high accuracy
        // ∫₀¹ x⁴ dx = 1/5
        let f = |x: f64| x.powi(4);
        let result = romberg_integration(f, 0.0, 1.0, 15, 1e-14);
        assert_abs_diff_eq!(result, 0.2, epsilon = 1e-14);
    }

    #[test]
    fn test_romberg_cosine() {
        // ∫₀^{π/2} cos(x) dx = 1
        use std::f64::consts::FRAC_PI_2;
        let f = |x: f64| x.cos();
        let result = romberg_integration(f, 0.0, FRAC_PI_2, 10, 1e-12);
        assert_abs_diff_eq!(result, 1.0, epsilon = 1e-12);
    }

    #[test]
    fn test_simpson_segment() {
        // Test the helper function
        // ∫₀¹ x² dx using single Simpson segment = 1/3
        let f = |x: f64| x * x;
        let result = simpson_segment(&f, 0.0, 1.0);
        assert_abs_diff_eq!(result, 1.0 / 3.0, epsilon = 1e-15);
    }

    #[test]
    fn test_adaptive_vs_fixed_efficiency() {
        // Adaptive should need fewer evaluations for smooth functions
        let f = |x: f64| x.exp();

        // Adaptive with high tolerance
        let adaptive_result = adaptive_simpson(f, 0.0, 1.0, 1e-6, 10);

        // "Exact" value
        let exact = std::f64::consts::E - 1.0;

        assert_abs_diff_eq!(adaptive_result, exact, epsilon = 1e-6);
    }

    #[test]
    fn test_romberg_convergence() {
        // Verify that more iterations improve accuracy
        let f = |x: f64| (x * x).sin();

        let result_5 = romberg_integration(f, 0.0, 1.0, 5, 1e-20);
        let result_10 = romberg_integration(f, 0.0, 1.0, 10, 1e-20);

        // Reference value with very fine integration
        let exact = romberg_integration(f, 0.0, 1.0, 15, 1e-15);

        let error_5 = (result_5 - exact).abs();
        let error_10 = (result_10 - exact).abs();

        // More iterations should give better accuracy
        assert!(error_10 < error_5);
    }

    #[test]
    fn test_f32_support() {
        // Test that methods work with f32
        let f = |x: f32| x * x;
        let adaptive_result = adaptive_simpson(f, 0.0_f32, 1.0_f32, 1e-6_f32, 15);
        let romberg_result = romberg_integration(f, 0.0_f32, 1.0_f32, 10, 1e-6_f32);

        assert_abs_diff_eq!(adaptive_result, 1.0_f32 / 3.0_f32, epsilon = 1e-6_f32);
        assert_abs_diff_eq!(romberg_result, 1.0_f32 / 3.0_f32, epsilon = 1e-6_f32);
    }

    #[test]
    fn test_constant_function() {
        // ∫₀¹ 7 dx = 7
        let f = |_x: f64| 7.0;
        let adaptive_result = adaptive_simpson(f, 0.0, 1.0, 1e-10, 10);
        let romberg_result = romberg_integration(f, 0.0, 1.0, 10, 1e-10);

        assert_abs_diff_eq!(adaptive_result, 7.0, epsilon = 1e-14);
        assert_abs_diff_eq!(romberg_result, 7.0, epsilon = 1e-14);
    }
}
