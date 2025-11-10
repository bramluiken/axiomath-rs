//! Gamma function and related functions.
//!
//! This module implements the gamma function Γ(x), its logarithm, and incomplete gamma functions,
//! which are essential for probability distributions and combinatorics.
//!
//! # Mathematical Background
//!
//! The gamma function is defined as:
//! ```text
//! Γ(x) = ∫₀^∞ t^(x-1) e^(-t) dt    for x > 0
//! ```
//!
//! Key properties:
//! - Γ(n+1) = n! for non-negative integers n
//! - Γ(x+1) = x·Γ(x) (recursive property)
//! - Γ(1/2) = √π
//!
//! # Implementation
//!
//! Uses the Lanczos approximation for high accuracy across all domains.
//!
//! # References
//!
//! - Lanczos, C., "A Precision Approximation of the Gamma Function", SIAM JNA 1, 1964
//! - Press et al., "Numerical Recipes", Section 6.1

use crate::core::Float;
use num_traits::Float as NumFloat;

/// Computes the gamma function Γ(x).
///
/// The gamma function extends the factorial function to real and complex numbers:
/// ```text
/// Γ(n) = (n-1)!  for positive integers n
/// Γ(x) = ∫₀^∞ t^(x-1) e^(-t) dt  for x > 0
/// ```
///
/// # Algorithm
///
/// Uses the Lanczos approximation with g=7 and n=9 coefficients,
/// achieving accuracy better than 2e-10 for all x > 0.
///
/// For x < 0, uses the reflection formula:
/// ```text
/// Γ(x) = π / (sin(πx) Γ(1-x))
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::special::gamma_function;
///
/// // Γ(1) = 0! = 1
/// assert!((gamma_function(1.0) - 1.0).abs() < 1e-10);
///
/// // Γ(5) = 4! = 24
/// assert!((gamma_function(5.0) - 24.0).abs() < 1e-10);
///
/// // Γ(1/2) = √π
/// assert!((gamma_function(0.5) - std::f64::consts::PI.sqrt()).abs() < 1e-10);
/// ```
///
/// # Special Values
///
/// - Γ(0) = ∞
/// - Γ(-n) = ∞ for non-negative integers n
/// - Γ(∞) = ∞
///
/// # Complexity
///
/// Time: O(1) - Fixed number of operations
/// Space: O(1)
pub fn gamma_function<T: Float>(x: T) -> T {
    // Handle special cases
    if x.is_nan() {
        return T::nan();
    }
    if x.is_infinite() {
        return if x > T::zero() { T::infinity() } else { T::nan() };
    }

    // For negative integers and zero, gamma is undefined (pole)
    if x <= T::zero() {
        let x_floor = x.floor();
        if x == x_floor {
            return T::infinity();
        }
    }

    // For x < 0.5, use reflection formula: Γ(x) = π / (sin(πx) Γ(1-x))
    let half = T::from(0.5).unwrap();
    if x < half {
        let pi = T::from(std::f64::consts::PI).unwrap();
        let pi_x = pi * x;
        return pi / (pi_x.sin() * gamma_function(T::one() - x));
    }

    // Lanczos approximation
    lanczos_gamma(x)
}

/// Lanczos approximation of the gamma function.
///
/// Uses g=7 with 9 coefficients for high precision.
fn lanczos_gamma<T: Float>(x: T) -> T {
    // Lanczos coefficients for g=7, n=9
    let coefficients = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];

    let g = T::from(7.0).unwrap();
    let two_pi = T::from(2.0 * std::f64::consts::PI).unwrap();

    // Compute Lanczos sum
    let x_minus_one = x - T::one();
    let mut ag = T::from(coefficients[0]).unwrap();

    for (i, &coef) in coefficients.iter().enumerate().skip(1) {
        let i_float = T::from(i).unwrap();
        ag = ag + T::from(coef).unwrap() / (x_minus_one + i_float);
    }

    let term1 = (two_pi).sqrt();
    let term2 = (x_minus_one + g + T::from(0.5).unwrap()).powf(x_minus_one + T::from(0.5).unwrap());
    let term3 = (-(x_minus_one + g + T::from(0.5).unwrap())).exp();

    term1 * term2 * term3 * ag
}

/// Computes the natural logarithm of the gamma function ln(Γ(x)).
///
/// This function is more numerically stable than computing log(gamma(x))
/// when Γ(x) would overflow or underflow.
///
/// # Mathematical Properties
///
/// - ln(Γ(x+1)) = ln(Γ(x)) + ln(x)
/// - ln(Γ(n)) = ln((n-1)!) for positive integers n
///
/// # Examples
///
/// ```
/// use axiomath::special::log_gamma_function;
///
/// // ln(Γ(1)) = ln(1) = 0
/// assert!((log_gamma_function(1.0) - 0.0).abs() < 1e-10);
///
/// // ln(Γ(5)) = ln(24)
/// assert!((log_gamma_function(5.0) - 24.0_f64.ln()).abs() < 1e-10);
/// ```
///
/// # Complexity
///
/// Time: O(1) - Fixed number of operations
/// Space: O(1)
pub fn log_gamma_function<T: Float>(x: T) -> T {
    // Handle special cases
    if x.is_nan() || x <= T::zero() {
        return T::nan();
    }
    if x.is_infinite() {
        return T::infinity();
    }

    // For x < 0.5, use reflection formula
    let half = T::from(0.5).unwrap();
    if x < half {
        let pi = T::from(std::f64::consts::PI).unwrap();
        let pi_x = pi * x;
        return pi.ln() - pi_x.sin().ln() - log_gamma_function(T::one() - x);
    }

    // Lanczos approximation for ln(Γ(x))
    let coefficients = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];

    let g = T::from(7.0).unwrap();
    let two_pi = T::from(2.0 * std::f64::consts::PI).unwrap();

    let x_minus_one = x - T::one();
    let mut ag = T::from(coefficients[0]).unwrap();

    for (i, &coef) in coefficients.iter().enumerate().skip(1) {
        let i_float = T::from(i).unwrap();
        ag = ag + T::from(coef).unwrap() / (x_minus_one + i_float);
    }

    let term1 = (two_pi).ln() / T::from(2.0).unwrap();
    let term2 = (x_minus_one + T::from(0.5).unwrap()) * (x_minus_one + g + T::from(0.5).unwrap()).ln();
    let term3 = -(x_minus_one + g + T::from(0.5).unwrap());
    let term4 = ag.ln();

    term1 + term2 + term3 + term4
}

/// Computes the lower incomplete gamma function γ(a, x).
///
/// The lower incomplete gamma function is defined as:
/// ```text
/// γ(a, x) = ∫₀ˣ t^(a-1) e^(-t) dt
/// ```
///
/// It satisfies: γ(a, x) + Γ(a, x) = Γ(a)
///
/// # Algorithm
///
/// Uses series representation for x < a+1, otherwise uses continued fraction
/// for the upper incomplete gamma and computes γ(a,x) = Γ(a) - Γ(a,x).
///
/// # Examples
///
/// ```
/// use axiomath::special::{incomplete_gamma_lower, gamma_function};
///
/// let a = 2.0;
/// let x = 1.0;
/// let result = incomplete_gamma_lower(a, x);
///
/// // For a=2, x=1: γ(2,1) = 1 - 2e^(-1) ≈ 0.2642411177
/// assert!((result - 0.2642411177).abs() < 1e-8);
/// ```
///
/// # Complexity
///
/// Time: O(1) - Fixed number of series/continued fraction terms
/// Space: O(1)
pub fn incomplete_gamma_lower<T: Float>(a: T, x: T) -> T {
    // Handle special cases
    if a.is_nan() || x.is_nan() || a <= T::zero() {
        return T::nan();
    }
    if x < T::zero() {
        return T::nan();
    }
    if x.is_zero() {
        return T::zero();
    }
    if x.is_infinite() {
        return gamma_function(a);
    }

    // Choose method based on x relative to a
    if x < a + T::one() {
        // Use series representation
        incomplete_gamma_series(a, x)
    } else {
        // Use continued fraction for upper incomplete gamma
        let gamma_a = gamma_function(a);
        gamma_a - incomplete_gamma_upper_cf(a, x)
    }
}

/// Series representation for lower incomplete gamma function.
fn incomplete_gamma_series<T: Float>(a: T, x: T) -> T {
    let max_iterations = 100;
    let epsilon = T::epsilon() * T::from(100.0).unwrap();

    let mut sum = T::one() / a;
    let mut term = T::one() / a;
    let mut n = T::one();

    for _ in 0..max_iterations {
        term = term * x / (a + n);
        sum = sum + term;

        if term.abs() < epsilon * sum.abs() {
            break;
        }

        n = n + T::one();
    }

    let exp_term = (-x + a * x.ln()).exp();
    exp_term * sum
}

/// Continued fraction representation for upper incomplete gamma function.
fn incomplete_gamma_upper_cf<T: Float>(a: T, x: T) -> T {
    let max_iterations = 100;
    let epsilon = T::epsilon() * T::from(100.0).unwrap();

    let mut b = x + T::one() - a;
    let mut c = T::one() / T::epsilon();
    let mut d = T::one() / b;
    let mut h = d;

    for i in 1..max_iterations {
        let i_float = T::from(i).unwrap();
        let an = -i_float * (i_float - a);
        b = b + T::from(2.0).unwrap();

        d = an * d + b;
        c = b + an / c;

        if d.abs() < epsilon {
            d = epsilon;
        }
        if c.abs() < epsilon {
            c = epsilon;
        }

        d = T::one() / d;
        let delta = c * d;
        h = h * delta;

        if (delta - T::one()).abs() < epsilon {
            break;
        }
    }

    let exp_term = (-x + a * x.ln()).exp();
    exp_term * h
}

/// Computes the regularized lower incomplete gamma function P(a, x) = γ(a, x) / Γ(a).
///
/// This is the cumulative distribution function (CDF) of the gamma distribution
/// with shape parameter a and scale parameter 1.
///
/// # Mathematical Properties
///
/// - P(a, 0) = 0
/// - P(a, ∞) = 1
/// - 0 ≤ P(a, x) ≤ 1 for all x ≥ 0
///
/// # Examples
///
/// ```
/// use axiomath::special::regularized_incomplete_gamma;
///
/// let a = 2.0;
/// let x = 2.0;
/// let result = regularized_incomplete_gamma(a, x);
///
/// // Result should be between 0 and 1
/// assert!(result >= 0.0 && result <= 1.0);
/// ```
pub fn regularized_incomplete_gamma<T: Float>(a: T, x: T) -> T {
    let gamma_a = gamma_function(a);
    incomplete_gamma_lower(a, x) / gamma_a
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_gamma_integers() {
        // Γ(n) = (n-1)!
        assert_abs_diff_eq!(gamma_function(1.0), 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(gamma_function(2.0), 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(gamma_function(3.0), 2.0, epsilon = 1e-10);
        assert_abs_diff_eq!(gamma_function(4.0), 6.0, epsilon = 1e-10);
        assert_abs_diff_eq!(gamma_function(5.0), 24.0, epsilon = 1e-10);
        assert_abs_diff_eq!(gamma_function(6.0), 120.0, epsilon = 1e-10);
    }

    #[test]
    fn test_gamma_half_integers() {
        // Γ(1/2) = √π
        assert_abs_diff_eq!(
            gamma_function(0.5),
            std::f64::consts::PI.sqrt(),
            epsilon = 1e-10
        );

        // Γ(3/2) = √π/2
        assert_abs_diff_eq!(
            gamma_function(1.5),
            std::f64::consts::PI.sqrt() / 2.0,
            epsilon = 1e-10
        );

        // Γ(5/2) = 3√π/4
        assert_abs_diff_eq!(
            gamma_function(2.5),
            3.0 * std::f64::consts::PI.sqrt() / 4.0,
            epsilon = 1e-10
        );
    }

    #[test]
    fn test_gamma_recursive_property() {
        // Γ(x+1) = x·Γ(x)
        let values = [0.5, 1.0, 2.5, 3.7, 5.2];

        for x in values {
            let gamma_x = gamma_function(x);
            let gamma_x_plus_1 = gamma_function(x + 1.0);
            assert_abs_diff_eq!(gamma_x_plus_1, x * gamma_x, epsilon = 1e-9);
        }
    }

    #[test]
    fn test_gamma_special_values() {
        assert!(gamma_function(f64::NAN).is_nan());
        assert_eq!(gamma_function(f64::INFINITY), f64::INFINITY);
        assert_eq!(gamma_function(0.0), f64::INFINITY);
        assert_eq!(gamma_function(-1.0), f64::INFINITY);
        assert_eq!(gamma_function(-2.0), f64::INFINITY);
    }

    #[test]
    fn test_gamma_negative_non_integers() {
        // For negative non-integers, gamma should be finite
        let result = gamma_function(-0.5);
        assert!(NumFloat::is_finite(result));

        let result2 = gamma_function(-1.5);
        assert!(NumFloat::is_finite(result2));
    }

    #[test]
    fn test_log_gamma_integers() {
        // ln(Γ(n)) = ln((n-1)!)
        assert_abs_diff_eq!(log_gamma_function(1.0), 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(log_gamma_function(2.0), 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(log_gamma_function(3.0), (2.0_f64).ln(), epsilon = 1e-10);
        assert_abs_diff_eq!(log_gamma_function(4.0), (6.0_f64).ln(), epsilon = 1e-10);
        assert_abs_diff_eq!(log_gamma_function(5.0), (24.0_f64).ln(), epsilon = 1e-10);
    }

    #[test]
    fn test_log_gamma_consistency() {
        let values = [0.5, 1.0, 2.5, 5.0, 10.0];

        for x in values {
            let gamma_x = gamma_function(x);
            let log_gamma_x = log_gamma_function(x);
            assert_abs_diff_eq!(log_gamma_x, NumFloat::ln(gamma_x), epsilon = 1e-9);
        }
    }

    #[test]
    fn test_log_gamma_large_values() {
        // For large x, gamma(x) overflows but log_gamma(x) should work
        let x = 100.0;
        let log_gamma_x = log_gamma_function(x);
        assert!(NumFloat::is_finite(log_gamma_x));
        assert!(log_gamma_x > 0.0);
    }

    #[test]
    fn test_incomplete_gamma_limits() {
        let a = 2.0;

        // γ(a, 0) = 0
        assert_abs_diff_eq!(incomplete_gamma_lower(a, 0.0), 0.0, epsilon = 1e-10);

        // γ(a, ∞) = Γ(a)
        let gamma_a = gamma_function(a);
        assert_abs_diff_eq!(
            incomplete_gamma_lower(a, f64::INFINITY),
            gamma_a,
            epsilon = 1e-10
        );
    }

    #[test]
    fn test_incomplete_gamma_known_values() {
        // For a=1: γ(1,x) = 1 - e^(-x)
        let x = 1.0;
        let result = incomplete_gamma_lower(1.0, x);
        let expected = 1.0 - NumFloat::exp(-x);
        assert_abs_diff_eq!(result, expected, epsilon = 1e-8);
    }

    #[test]
    fn test_regularized_incomplete_gamma() {
        let a = 2.0;
        let x = 1.0;
        let p = regularized_incomplete_gamma(a, x);

        // Should be between 0 and 1
        assert!(p >= 0.0 && p <= 1.0);

        // P(a, 0) = 0
        assert_abs_diff_eq!(regularized_incomplete_gamma(a, 0.0), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_regularized_incomplete_gamma_cdf_properties() {
        let a = 3.0;

        // Should be monotonically increasing
        let x_values = [0.5, 1.0, 2.0, 5.0, 10.0];
        let mut prev = 0.0;

        for x in x_values {
            let p = regularized_incomplete_gamma(a, x);
            assert!(p > prev);
            assert!(p <= 1.0);
            prev = p;
        }
    }

    #[test]
    fn test_gamma_f32() {
        let result: f32 = gamma_function(5.0_f32);
        assert_abs_diff_eq!(result, 24.0_f32, epsilon = 1e-4);
    }
}
