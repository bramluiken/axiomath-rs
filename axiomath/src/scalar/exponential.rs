//! Exponential functions implemented using Taylor series.
//!
//! This module provides exponential functions using mathematical series:
//! - e^x using Taylor series expansion
//! - 2^x and other base exponentials

use crate::core::traits::Float;

/// Computes e^x using Taylor series expansion.
///
/// # Algorithm
///
/// Taylor series for e^x:
/// ```text
/// e^x = 1 + x + x²/2! + x³/3! + x⁴/4! + ...
///     = Σ(x^n / n!) for n = 0 to ∞
/// ```
///
/// For large |x|, uses range reduction:
/// - e^x = (e^(x/2))² for |x| > 1
///
/// This ensures faster convergence of the series.
///
/// # Examples
///
/// ```
/// use axiomath::scalar::exponential::exponential;
///
/// let result = exponential(0.0_f64);
/// assert!((result - 1.0).abs() < 1e-2);
///
/// let result = exponential(1.0_f64);
/// assert!((result - std::f64::consts::E).abs() < 1e-2);
///
/// let result = exponential(2.0_f64);
/// assert!((result - 7.389056098930650).abs() < 1e-2);
/// ```
#[inline]
pub fn exponential<T>(x: T) -> T
where
    T: Float,
{
    // Handle special cases
    if x.is_zero() {
        return T::one();
    }
    if x.is_nan() {
        return x;
    }
    if x.is_infinite() {
        if x > T::zero() {
            return T::infinity();
        } else {
            return T::zero();
        }
    }

    // Range reduction for faster convergence
    // If |x| > 1, use e^x = (e^(x/2))^2
    let two = T::one() + T::one();
    if x.abs() > T::one() {
        let half_exp = exponential(x / two);
        return half_exp * half_exp;
    }

    // Taylor series: e^x = Σ(x^n / n!)
    let mut sum = T::one();
    let mut term = T::one();
    let epsilon = T::epsilon() * (T::one() + T::one()); // 2 * epsilon
    let max_iterations = 100;

    for n in 1..max_iterations {
        term = term * x / T::from(n).expect("Failed to convert n to T");
        sum = sum + term;

        // Check for convergence
        if term.abs() <= epsilon * sum.abs() {
            break;
        }
    }

    sum
}

/// Computes 2^x.
///
/// # Algorithm
///
/// Uses the identity: 2^x = e^(x * ln(2))
///
/// # Examples
///
/// ```
/// use axiomath::scalar::exponential::exp2;
///
/// let result = exp2(0.0_f64);
/// assert!((result - 1.0).abs() < 1e-2);
///
/// let result = exp2(10.0_f64);
/// assert!((result - 1024.0).abs() < 1e-2);
/// ```
#[inline]
pub fn exp2<T>(x: T) -> T
where
    T: Float,
{
    // 2^x = e^(x * ln(2))
    // ln(2) ≈ 0.693147180559945309417232121458
    let ln2 = T::from(0.693147180559945309417232121458)
        .expect("Failed to convert ln(2) to T");
    exponential(x * ln2)
}

/// Computes 10^x.
///
/// # Algorithm
///
/// Uses the identity: 10^x = e^(x * ln(10))
///
/// # Examples
///
/// ```
/// use axiomath::scalar::exponential::exp10;
///
/// let result = exp10(0.0_f64);
/// assert!((result - 1.0).abs() < 1e-2);
///
/// let result = exp10(3.0_f64);
/// assert!((result - 1000.0).abs() < 1e-2);
/// ```
#[inline]
pub fn exp10<T>(x: T) -> T
where
    T: Float,
{
    // 10^x = e^(x * ln(10))
    // ln(10) ≈ 2.302585092994045684017991454684
    let ln10 = T::from(2.302585092994045684017991454684)
        .expect("Failed to convert ln(10) to T");
    exponential(x * ln10)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_exponential_zero() {
        assert_abs_diff_eq!(exponential(0.0_f64), 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(exponential(0.0_f32), 1.0, epsilon = 1e-6);
    }

    #[test]
    fn test_exponential_one() {
        // e^1 = e ≈ 2.718281828459045
        assert_abs_diff_eq!(
            exponential(1.0_f64),
            std::f64::consts::E,
            epsilon = 1e-10
        );
    }

    #[test]
    fn test_exponential_positive() {
        // e^2 ≈ 7.389056098930650
        assert_abs_diff_eq!(exponential(2.0_f64), 7.389056098930650, epsilon = 1e-9);

        // e^3 ≈ 20.085536923187668
        assert_abs_diff_eq!(exponential(3.0_f64), 20.085536923187668, epsilon = 1e-8);

        // e^0.5 ≈ 1.648721270700128
        assert_abs_diff_eq!(exponential(0.5_f64), 1.648721270700128, epsilon = 1e-10);
    }

    #[test]
    fn test_exponential_negative() {
        // e^(-1) ≈ 0.367879441171442
        assert_abs_diff_eq!(exponential(-1.0_f64), 0.367879441171442, epsilon = 1e-10);

        // e^(-2) ≈ 0.135335283236613
        assert_abs_diff_eq!(exponential(-2.0_f64), 0.135335283236613, epsilon = 1e-10);
    }

    #[test]
    fn test_exponential_large() {
        // e^10 ≈ 22026.465794806718
        assert_abs_diff_eq!(exponential(10.0_f64), 22026.465794806718, epsilon = 1e-4);
    }

    #[test]
    fn test_exponential_small() {
        // e^0.1 ≈ 1.105170918075648
        assert_abs_diff_eq!(exponential(0.1_f64), 1.105170918075648, epsilon = 1e-10);

        // e^(-0.1) ≈ 0.904837418035960
        assert_abs_diff_eq!(exponential(-0.1_f64), 0.904837418035960, epsilon = 1e-10);
    }

    #[test]
    fn test_exponential_special_values() {
        assert!(exponential(f64::NAN).is_nan());
        assert!(exponential(f64::INFINITY).is_infinite());
        assert!(exponential(f64::INFINITY) > 0.0);
        assert_abs_diff_eq!(exponential(f64::NEG_INFINITY), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_exp2_basic() {
        assert_abs_diff_eq!(exp2(0.0_f64), 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(exp2(1.0_f64), 2.0, epsilon = 1e-10);
        assert_abs_diff_eq!(exp2(2.0_f64), 4.0, epsilon = 1e-9);
        assert_abs_diff_eq!(exp2(10.0_f64), 1024.0, epsilon = 1e-7);
    }

    #[test]
    fn test_exp2_negative() {
        assert_abs_diff_eq!(exp2(-1.0_f64), 0.5, epsilon = 1e-10);
        assert_abs_diff_eq!(exp2(-2.0_f64), 0.25, epsilon = 1e-10);
    }

    #[test]
    fn test_exp10_basic() {
        assert_abs_diff_eq!(exp10(0.0_f64), 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(exp10(1.0_f64), 10.0, epsilon = 1e-9);
        assert_abs_diff_eq!(exp10(2.0_f64), 100.0, epsilon = 1e-8);
        assert_abs_diff_eq!(exp10(3.0_f64), 1000.0, epsilon = 1e-7);
    }

    #[test]
    fn test_exp10_negative() {
        assert_abs_diff_eq!(exp10(-1.0_f64), 0.1, epsilon = 1e-10);
        assert_abs_diff_eq!(exp10(-2.0_f64), 0.01, epsilon = 1e-10);
    }

    #[test]
    fn test_exponential_identity() {
        // e^(a+b) = e^a * e^b
        let a = 1.5_f64;
        let b = 2.3_f64;
        let lhs = exponential(a + b);
        let rhs = exponential(a) * exponential(b);
        assert_abs_diff_eq!(lhs, rhs, epsilon = 1e-8);
    }

    #[test]
    fn test_exponential_inverse() {
        // e^x * e^(-x) = 1
        for x in [0.5, 1.0, 2.0, 5.0] {
            let product = exponential(x) * exponential(-x);
            assert_abs_diff_eq!(product, 1.0, epsilon = 1e-9);
        }
    }
}
