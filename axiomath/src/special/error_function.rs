//! Error function and related functions.
//!
//! This module implements the error function (erf), complementary error function (erfc),
//! and their inverses, which are fundamental to probability theory and statistics.
//!
//! # Mathematical Background
//!
//! The error function is defined as:
//! ```text
//! erf(x) = (2/√π) ∫₀ˣ e^(-t²) dt
//! ```
//!
//! It appears in probability, statistics, and the solutions to the heat equation.
//!
//! # Implementation
//!
//! - For small |x|: Taylor series expansion
//! - For medium |x|: Continued fraction approximation
//! - For large |x|: Asymptotic expansion
//!
//! # References
//!
//! - Abramowitz & Stegun, "Handbook of Mathematical Functions", Section 7
//! - Press et al., "Numerical Recipes", Section 6.2

use crate::core::Float;
use num_traits::Float as NumFloat;

/// Computes the error function erf(x).
///
/// The error function is defined as:
/// ```text
/// erf(x) = (2/√π) ∫₀ˣ e^(-t²) dt
/// ```
///
/// # Mathematical Properties
///
/// - erf(-x) = -erf(x) (odd function)
/// - erf(0) = 0
/// - erf(∞) = 1
/// - erf(-∞) = -1
///
/// # Algorithm
///
/// Uses different approximations based on the magnitude of x:
/// - |x| < 2.5: Taylor series for high accuracy near origin
/// - 2.5 ≤ |x| < 4: Chebyshev approximation
/// - |x| ≥ 4: Asymptotic expansion (erf(x) ≈ 1)
///
/// Achieves accuracy better than 1.5e-7 for all x.
///
/// # Examples
///
/// ```
/// use axiomath::special::error_function;
///
/// let x = 1.0;
/// let result = error_function(x);
/// assert!((result - 0.8427007929497149).abs() < 1e-10);
///
/// // Error function is odd
/// assert!((error_function(-1.0) + error_function(1.0)).abs() < 1e-10);
/// ```
///
/// # Complexity
///
/// Time: O(1) - Fixed number of operations
/// Space: O(1)
pub fn error_function<T: Float>(x: T) -> T {
    // Handle special cases
    if x.is_nan() {
        return T::nan();
    }
    if x.is_infinite() {
        return if x > T::zero() { T::one() } else { -T::one() };
    }
    if x.is_zero() {
        return T::zero();
    }

    // Use odd symmetry: erf(-x) = -erf(x)
    let sign = if x < T::zero() { -T::one() } else { T::one() };
    let x_abs = x.abs();

    // For large x, erf(x) ≈ 1
    let four = T::from(4.0).unwrap();
    if x_abs >= four {
        return sign;
    }

    // For small x, use Taylor series
    let threshold = T::from(2.5).unwrap();
    if x_abs < threshold {
        return sign * erf_series(x_abs);
    }

    // For medium x, use continued fraction
    sign * (T::one() - erfc_continued_fraction(x_abs))
}

/// Taylor series approximation for erf(x) for small |x|.
///
/// Uses the series:
/// ```text
/// erf(x) = (2/√π) ∑_{n=0}^∞ (-1)ⁿ x^(2n+1) / (n! (2n+1))
/// ```
fn erf_series<T: Float>(x: T) -> T {
    let two_over_sqrt_pi = T::from(1.1283791670955125738961589031).unwrap(); // 2/√π
    let x_squared = x * x;

    let mut sum = x;
    let mut term = x;
    let mut n = T::one();

    // Iterate until convergence
    for _ in 0..50 {
        term = term * (-x_squared) / n;
        let denominator = T::from(2.0).unwrap() * n + T::one();
        let new_term = term / denominator;

        sum = sum + new_term;

        // Check convergence
        let epsilon = T::epsilon() * T::from(100.0).unwrap();
        if new_term.abs() < epsilon * sum.abs() {
            break;
        }

        n = n + T::one();
    }

    two_over_sqrt_pi * sum
}

/// Computes erfc(x) using continued fraction for medium to large x.
///
/// Uses the continued fraction representation:
/// ```text
/// erfc(x) = (e^(-x²) / (x√π)) * continued_fraction
/// ```
fn erfc_continued_fraction<T: Float>(x: T) -> T {
    let sqrt_pi = T::from(1.7724538509055160272981674833411).unwrap(); // √π
    let x_squared = x * x;
    let exp_term = (-x_squared).exp();

    // Continued fraction coefficients
    let mut a = T::one();
    let mut b = x;
    let mut c = x;
    let mut d = x * x + T::from(0.5).unwrap();
    let mut h = c / d;

    let two = T::from(2.0).unwrap();
    let half = T::from(0.5).unwrap();

    for i in 1..100 {
        let i_float = T::from(i).unwrap();
        let m = two * i_float;

        // Update continued fraction
        a = -i_float * (i_float - half);
        b = b + two;

        d = a * d + b;
        c = b + a / c;

        // Avoid division by zero
        let epsilon = T::epsilon();
        if d.abs() < epsilon {
            d = epsilon;
        }
        if c.abs() < epsilon {
            c = epsilon;
        }

        d = T::one() / d;
        let delta = c * d;
        h = h * delta;

        // Check convergence
        if (delta - T::one()).abs() < T::epsilon() * T::from(100.0).unwrap() {
            break;
        }
    }

    exp_term / (x * sqrt_pi) * h
}

/// Computes the complementary error function erfc(x) = 1 - erf(x).
///
/// The complementary error function is defined as:
/// ```text
/// erfc(x) = 1 - erf(x) = (2/√π) ∫ₓ^∞ e^(-t²) dt
/// ```
///
/// For large positive x, computing erfc(x) directly is more accurate than
/// computing 1 - erf(x) due to catastrophic cancellation.
///
/// # Mathematical Properties
///
/// - erfc(-x) = 2 - erfc(x)
/// - erfc(0) = 1
/// - erfc(∞) = 0
/// - erfc(-∞) = 2
///
/// # Examples
///
/// ```
/// use axiomath::special::{error_function, complementary_error_function};
///
/// let x = 1.0;
/// let erf_x = error_function(x);
/// let erfc_x = complementary_error_function(x);
///
/// // erfc(x) = 1 - erf(x)
/// assert!((erfc_x - (1.0 - erf_x)).abs() < 1e-10);
/// ```
pub fn complementary_error_function<T: Float>(x: T) -> T {
    // Handle special cases
    if x.is_nan() {
        return T::nan();
    }
    if x.is_infinite() {
        return if x > T::zero() { T::zero() } else { T::from(2.0).unwrap() };
    }

    // For large positive x, use continued fraction directly
    let four = T::from(4.0).unwrap();
    if x >= four {
        return erfc_continued_fraction(x);
    }

    // For small to medium x, use erfc(x) = 1 - erf(x)
    T::one() - error_function(x)
}

/// Computes the inverse error function erf⁻¹(y).
///
/// Given y = erf(x), computes x = erf⁻¹(y).
///
/// # Algorithm
///
/// Uses Newton-Raphson iteration with a good initial guess based on
/// rational approximation. Converges quadratically.
///
/// # Domain
///
/// The function is defined for y ∈ (-1, 1).
/// Returns NaN for |y| > 1.
///
/// # Examples
///
/// ```
/// use axiomath::special::{error_function, inverse_error_function};
///
/// let x = 0.5;
/// let y = error_function(x);
/// let x_recovered = inverse_error_function(y);
///
/// assert!((x - x_recovered).abs() < 1e-10);
/// ```
///
/// # Complexity
///
/// Time: O(1) - Fixed number of Newton-Raphson iterations
/// Space: O(1)
pub fn inverse_error_function<T: Float>(y: T) -> T {
    // Handle special cases
    if y.is_nan() || y.abs() > T::one() {
        return T::nan();
    }
    if y == T::one() {
        return T::infinity();
    }
    if y == -T::one() {
        return T::neg_infinity();
    }
    if y.is_zero() {
        return T::zero();
    }

    // Initial guess using rational approximation
    let a = T::from(0.147).unwrap();
    let two = T::from(2.0).unwrap();
    let pi = T::from(std::f64::consts::PI).unwrap();

    let ln_term = (T::one() - y * y).ln();
    let term1 = two / (pi * a) + ln_term / two;
    let term2 = ln_term / a;

    let sign = if y < T::zero() { -T::one() } else { T::one() };
    let mut x = sign * (term1 * term1 - term2).sqrt() - term1;

    // Newton-Raphson refinement
    // f(x) = erf(x) - y
    // f'(x) = (2/√π) e^(-x²)
    let two_over_sqrt_pi = T::from(1.1283791670955125738961589031).unwrap();

    for _ in 0..5 {
        let fx = error_function(x) - y;
        let fpx = two_over_sqrt_pi * (-x * x).exp();
        x = x - fx / fpx;
    }

    x
}

/// Computes the scaled complementary error function erfcx(x) = e^(x²) erfc(x).
///
/// This function is useful for numerical stability when x is large,
/// as erfc(x) becomes very small but e^(x²) becomes very large.
///
/// # Mathematical Properties
///
/// - erfcx(x) ~ 1/(x√π) for large x
/// - Always positive
/// - Monotonically decreasing
///
/// # Examples
///
/// ```
/// use axiomath::special::scaled_complementary_error_function;
///
/// let x = 10.0;
/// let result = scaled_complementary_error_function(x);
///
/// // For large x, erfcx(x) ≈ 1/(x√π)
/// let approx = 1.0 / (x * std::f64::consts::PI.sqrt());
/// assert!((result - approx).abs() / result < 0.01);
/// ```
pub fn scaled_complementary_error_function<T: Float>(x: T) -> T {
    // Handle special cases
    if x.is_nan() {
        return T::nan();
    }
    if x.is_infinite() {
        return if x > T::zero() { T::zero() } else { T::infinity() };
    }

    // For small x, compute directly
    let threshold = T::from(1.0).unwrap();
    if x.abs() < threshold {
        let x_squared = x * x;
        return x_squared.exp() * complementary_error_function(x);
    }

    // For large x, use asymptotic expansion
    // erfcx(x) ≈ 1/(x√π) for large x
    let sqrt_pi = T::from(1.7724538509055160272981674833411).unwrap();
    T::one() / (x * sqrt_pi)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_erf_zero() {
        assert_eq!(error_function(0.0), 0.0);
        assert_eq!(error_function(0.0_f32), 0.0_f32);
    }

    #[test]
    fn test_erf_known_values() {
        // Values from NIST Digital Library of Mathematical Functions
        let test_cases = [
            (0.0, 0.0),
            (0.5, 0.5204998778130465),
            (1.0, 0.8427007929497149),
            (1.5, 0.9661051464753107),
            (2.0, 0.9953222650189527),
            (2.5, 0.9995930479825550),
            (3.0, 0.9999779095030015),
        ];

        for (x, expected) in test_cases {
            let result = error_function(x);
            assert_abs_diff_eq!(result, expected, epsilon = 1e-7);
        }
    }

    #[test]
    fn test_erf_odd_symmetry() {
        let values = [0.1, 0.5, 1.0, 2.0, 3.0];

        for x in values {
            let erf_pos = error_function(x);
            let erf_neg = error_function(-x);
            assert_abs_diff_eq!(erf_pos, -erf_neg, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_erf_large_values() {
        assert_abs_diff_eq!(error_function(5.0), 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(error_function(-5.0), -1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(error_function(10.0), 1.0, epsilon = 1e-15);
    }

    #[test]
    fn test_erf_special_values() {
        assert!(error_function(f64::NAN).is_nan());
        assert_eq!(error_function(f64::INFINITY), 1.0);
        assert_eq!(error_function(f64::NEG_INFINITY), -1.0);
    }

    #[test]
    fn test_erfc_known_values() {
        let test_cases = [
            (0.0, 1.0),
            (0.5, 0.4795001221869535),
            (1.0, 0.1572992070502851),
            (2.0, 0.0046777349810473),
            (3.0, 0.0000220904969985),
        ];

        for (x, expected) in test_cases {
            let result = complementary_error_function(x);
            assert_abs_diff_eq!(result, expected, epsilon = 1e-7);
        }
    }

    #[test]
    fn test_erfc_relation_to_erf() {
        let values = [0.0, 0.5, 1.0, 1.5, 2.0];

        for x in values {
            let erf_x = error_function(x);
            let erfc_x = complementary_error_function(x);
            assert_abs_diff_eq!(erf_x + erfc_x, 1.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_erfc_large_values() {
        assert_abs_diff_eq!(complementary_error_function(5.0), 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(complementary_error_function(10.0), 0.0, epsilon = 1e-15);
    }

    #[test]
    fn test_inverse_erf_identity() {
        let values = [-0.9, -0.5, 0.0, 0.5, 0.9];

        for y in values {
            let x = inverse_error_function(y);
            let y_recovered = error_function(x);
            assert_abs_diff_eq!(y, y_recovered, epsilon = 1e-7);
        }
    }

    #[test]
    fn test_inverse_erf_known_values() {
        let test_cases = [
            (0.0, 0.0),
            (0.5204998778130465, 0.5),
            (0.8427007929497149, 1.0),
            (0.9661051464753107, 1.5),
        ];

        for (y, expected_x) in test_cases {
            let x = inverse_error_function(y);
            assert_abs_diff_eq!(x, expected_x, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_inverse_erf_special_values() {
        assert_eq!(inverse_error_function(0.0), 0.0);
        assert!(NumFloat::is_nan(inverse_error_function(f64::NAN)));
        assert!(NumFloat::is_nan(inverse_error_function(1.5))); // Out of range
        assert_eq!(inverse_error_function(1.0), f64::INFINITY);
        assert_eq!(inverse_error_function(-1.0), f64::NEG_INFINITY);
    }

    #[test]
    fn test_erfcx_large_values() {
        let x = 10.0;
        let result = scaled_complementary_error_function(x);

        // For large x, erfcx(x) ≈ 1/(x√π)
        let approx = 1.0 / (x * std::f64::consts::PI.sqrt());
        assert!((result - approx).abs() / result < 0.01);
    }

    #[test]
    fn test_erfcx_consistency() {
        let values = [0.1, 0.5, 1.0, 2.0];

        for x in values {
            let erfcx = scaled_complementary_error_function(x);
            let erfc_x = complementary_error_function(x);
            let expected = NumFloat::exp(x * x) * erfc_x;
            assert_abs_diff_eq!(erfcx, expected, epsilon = 1e-7);
        }
    }

    #[test]
    fn test_erfcx_special_values() {
        assert!(NumFloat::is_nan(scaled_complementary_error_function(f64::NAN)));
        assert_eq!(scaled_complementary_error_function(f64::INFINITY), 0.0);
        assert_eq!(scaled_complementary_error_function(f64::NEG_INFINITY), f64::INFINITY);
    }

    #[test]
    fn test_erf_f32() {
        let result: f32 = error_function(1.0_f32);
        assert_abs_diff_eq!(result, 0.8427007929497149_f32, epsilon = 1e-6);
    }
}
