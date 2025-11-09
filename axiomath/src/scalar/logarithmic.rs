//! Logarithmic functions implemented using series and Newton-Raphson method.
//!
//! This module provides logarithm functions using mathematical algorithms:
//! - Natural logarithm using series expansion
//! - Logarithms with different bases using change of base formula

use crate::core::traits::Float;
use crate::scalar::exponential::exponential;

/// Computes the natural logarithm ln(x) using series expansion and Newton-Raphson.
///
/// # Algorithm
///
/// For x close to 1, uses the Taylor series:
/// ```text
/// ln(1 + y) = y - y²/2 + y³/3 - y⁴/4 + ...
/// ```
///
/// For other values, uses range reduction:
/// - ln(x) = ln(m * 2^e) = ln(m) + e * ln(2)
/// where m ∈ [1/√2, √2] and e is an integer
///
/// Alternatively, uses Newton-Raphson: find y such that e^y = x
/// - y_{n+1} = y_n + (x - e^(y_n)) / e^(y_n)
///
/// # Examples
///
/// ```
/// use axiomath::scalar::logarithmic::natural_logarithm;
///
/// let result = natural_logarithm(1.0_f64);
/// assert!((result - 0.0).abs() < 1e-2);
///
/// let result = natural_logarithm(std::f64::consts::E);
/// assert!((result - 1.0).abs() < 1e-2);
/// ```
///
/// # Panics
///
/// Panics if x ≤ 0 (logarithm of non-positive number).
#[inline]
pub fn natural_logarithm<T>(x: T) -> T
where
    T: Float,
{
    // Handle special cases
    if x <= T::zero() {
        panic!("Cannot compute logarithm of non-positive number");
    }
    if x.is_nan() {
        return x;
    }
    if x.is_infinite() {
        return x;
    }
    if x == T::one() {
        return T::zero();
    }

    // For values close to 1, use Taylor series directly
    let two = T::one() + T::one();
    if x > (T::one() / two) && x < (T::one() + T::one()) {
        return ln_series(x);
    }

    // For other values, use range reduction: ln(x) = 2*ln(√x)
    // Reduce x to range [0.5, 2] by repeatedly halving
    let two = T::one() + T::one();
    let mut x_reduced = x;
    let mut power_of_two = T::zero();

    // Reduce large values
    while x_reduced > two {
        x_reduced = x_reduced / two;
        power_of_two = power_of_two + T::one();
    }

    // Reduce small values
    while x_reduced < (T::one() / two) {
        x_reduced = x_reduced * two;
        power_of_two = power_of_two - T::one();
    }

    // Now x_reduced is in [0.5, 2], use Taylor series
    let ln_reduced = ln_series(x_reduced);

    // ln(2) ≈ 0.693147180559945
    let ln2 = T::from(0.693147180559945309417232121458)
        .expect("Failed to convert ln(2) to T");

    ln_reduced + power_of_two * ln2
}

/// Helper function: computes ln(x) using Taylor series for x near 1.
///
/// Uses the series: ln(1 + u) = u - u²/2 + u³/3 - u⁴/4 + ...
/// where u = x - 1
#[inline]
fn ln_series<T>(x: T) -> T
where
    T: Float,
{
    let u = x - T::one();
    let two = T::one() + T::one();

    // For better convergence, if |u| > 0.5, use ln(x) = -ln(1/x)
    if u.abs() > (T::one() / two) {
        return -ln_series(T::one() / x);
    }

    let mut sum = T::zero();
    let mut term = u;
    let epsilon = T::epsilon() * (T::one() + T::one());
    let max_iterations = 100;

    for n in 1..max_iterations {
        let n_float = T::from(n).expect("Failed to convert n to T");
        let contribution = term / n_float;
        sum = sum + contribution;

        if contribution.abs() <= epsilon * sum.abs() {
            break;
        }

        term = -term * u; // Alternating series
    }

    sum
}

/// Computes logarithm with arbitrary base: log_base(x).
///
/// # Algorithm
///
/// Uses change of base formula:
/// ```text
/// log_b(x) = ln(x) / ln(b)
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::scalar::logarithmic::logarithm_base;
///
/// let result = logarithm_base(8.0_f64, 2.0_f64);
/// assert!((result - 3.0).abs() < 1e-2);
///
/// let result = logarithm_base(1000.0_f64, 10.0_f64);
/// assert!((result - 3.0).abs() < 1e-2);
/// ```
///
/// # Panics
///
/// - Panics if x ≤ 0
/// - Panics if base ≤ 0 or base == 1
#[inline]
pub fn logarithm_base<T>(x: T, base: T) -> T
where
    T: Float,
{
    if base <= T::zero() || base == T::one() {
        panic!("Logarithm base must be positive and not equal to 1");
    }

    natural_logarithm(x) / natural_logarithm(base)
}

/// Computes log₂(x) (logarithm base 2).
///
/// # Examples
///
/// ```
/// use axiomath::scalar::logarithmic::logarithm_2;
///
/// let result = logarithm_2(8.0_f64);
/// assert!((result - 3.0).abs() < 1e-2);
///
/// let result = logarithm_2(1024.0_f64);
/// assert!((result - 10.0).abs() < 1e-2);
/// ```
///
/// # Panics
///
/// Panics if x ≤ 0.
#[inline]
pub fn logarithm_2<T>(x: T) -> T
where
    T: Float,
{
    // ln(2) ≈ 0.693147180559945309417232121458
    let ln2 = T::from(0.693147180559945309417232121458)
        .expect("Failed to convert ln(2) to T");
    natural_logarithm(x) / ln2
}

/// Computes log₁₀(x) (logarithm base 10).
///
/// # Examples
///
/// ```
/// use axiomath::scalar::logarithmic::logarithm_10;
///
/// let result = logarithm_10(100.0_f64);
/// assert!((result - 2.0).abs() < 1e-2);
///
/// let result = logarithm_10(1000.0_f64);
/// assert!((result - 3.0).abs() < 1e-2);
/// ```
///
/// # Panics
///
/// Panics if x ≤ 0.
#[inline]
pub fn logarithm_10<T>(x: T) -> T
where
    T: Float,
{
    // ln(10) ≈ 2.302585092994045684017991454684
    let ln10 = T::from(2.302585092994045684017991454684)
        .expect("Failed to convert ln(10) to T");
    natural_logarithm(x) / ln10
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_natural_logarithm_one() {
        assert_abs_diff_eq!(natural_logarithm(1.0_f64), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_natural_logarithm_e() {
        // ln(e) = 1
        let e = std::f64::consts::E;
        assert_abs_diff_eq!(natural_logarithm(e), 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_natural_logarithm_basic() {
        // ln(2) ≈ 0.693147180559945
        assert_abs_diff_eq!(
            natural_logarithm(2.0_f64),
            0.693147180559945,
            epsilon = 1e-10
        );

        // ln(10) ≈ 2.302585092994046
        assert_abs_diff_eq!(
            natural_logarithm(10.0_f64),
            2.302585092994046,
            epsilon = 1e-9
        );

        // ln(0.5) ≈ -0.693147180559945
        assert_abs_diff_eq!(
            natural_logarithm(0.5_f64),
            -0.693147180559945,
            epsilon = 1e-10
        );
    }

    #[test]
    fn test_natural_logarithm_powers_of_e() {
        // ln(e²) = 2
        let e2 = std::f64::consts::E * std::f64::consts::E;
        assert_abs_diff_eq!(natural_logarithm(e2), 2.0, epsilon = 1e-9);

        // ln(e³) = 3
        let e3 = std::f64::consts::E.powi(3);
        assert_abs_diff_eq!(natural_logarithm(e3), 3.0, epsilon = 1e-9);
    }

    #[test]
    #[should_panic(expected = "Cannot compute logarithm of non-positive number")]
    fn test_natural_logarithm_zero() {
        natural_logarithm(0.0_f64);
    }

    #[test]
    #[should_panic(expected = "Cannot compute logarithm of non-positive number")]
    fn test_natural_logarithm_negative() {
        natural_logarithm(-1.0_f64);
    }

    #[test]
    fn test_logarithm_base_basic() {
        // log₂(8) = 3
        assert_abs_diff_eq!(logarithm_base(8.0_f64, 2.0), 3.0, epsilon = 1e-10);

        // log₁₀(1000) = 3
        assert_abs_diff_eq!(logarithm_base(1000.0_f64, 10.0), 3.0, epsilon = 1e-9);

        // log₅(125) = 3
        assert_abs_diff_eq!(logarithm_base(125.0_f64, 5.0), 3.0, epsilon = 1e-9);
    }

    #[test]
    fn test_logarithm_2_basic() {
        assert_abs_diff_eq!(logarithm_2(2.0_f64), 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(logarithm_2(4.0_f64), 2.0, epsilon = 1e-10);
        assert_abs_diff_eq!(logarithm_2(8.0_f64), 3.0, epsilon = 1e-10);
        assert_abs_diff_eq!(logarithm_2(1024.0_f64), 10.0, epsilon = 1e-9);
    }

    #[test]
    fn test_logarithm_10_basic() {
        assert_abs_diff_eq!(logarithm_10(10.0_f64), 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(logarithm_10(100.0_f64), 2.0, epsilon = 1e-10);
        assert_abs_diff_eq!(logarithm_10(1000.0_f64), 3.0, epsilon = 1e-9);
    }

    #[test]
    fn test_logarithm_identity() {
        // log_b(x * y) = log_b(x) + log_b(y)
        let x = 4.0_f64;
        let y = 8.0_f64;
        let base = 2.0_f64;

        let lhs = logarithm_base(x * y, base);
        let rhs = logarithm_base(x, base) + logarithm_base(y, base);
        assert_abs_diff_eq!(lhs, rhs, epsilon = 1e-9);
    }

    #[test]
    fn test_logarithm_inverse() {
        // exp(ln(x)) = x
        for x in [0.5, 1.0, 2.0, 5.0, 10.0] {
            let ln_x = natural_logarithm(x);
            let result = exponential(ln_x);
            assert_abs_diff_eq!(result, x, epsilon = 1e-9);
        }
    }

    #[test]
    fn test_logarithm_power_rule() {
        // ln(x^n) = n * ln(x)
        let x = 3.0_f64;
        let n = 4.0_f64;

        let lhs = natural_logarithm(x.powf(n));
        let rhs = n * natural_logarithm(x);
        assert_abs_diff_eq!(lhs, rhs, epsilon = 1e-9);
    }
}
