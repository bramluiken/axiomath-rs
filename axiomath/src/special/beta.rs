//! Beta function and related functions.
//!
//! This module implements the beta function B(a,b), its logarithm, and incomplete beta functions,
//! which are essential for probability distributions and Bayesian statistics.
//!
//! # Mathematical Background
//!
//! The beta function is defined as:
//! ```text
//! B(a,b) = ∫₀¹ t^(a-1) (1-t)^(b-1) dt    for a, b > 0
//! ```
//!
//! It is related to the gamma function by:
//! ```text
//! B(a,b) = Γ(a)Γ(b) / Γ(a+b)
//! ```
//!
//! # References
//!
//! - Abramowitz & Stegun, "Handbook of Mathematical Functions", Section 6.2
//! - Press et al., "Numerical Recipes", Section 6.4

use crate::core::Float;
use crate::special::gamma::{gamma_function, log_gamma_function};
use num_traits::Float as NumFloat;

/// Computes the beta function B(a,b).
///
/// The beta function is defined as:
/// ```text
/// B(a,b) = ∫₀¹ t^(a-1) (1-t)^(b-1) dt = Γ(a)Γ(b) / Γ(a+b)
/// ```
///
/// # Mathematical Properties
///
/// - B(a,b) = B(b,a) (symmetric)
/// - B(a,1) = 1/a
/// - B(1,1) = 1
/// - B(a,b) = (a-1)! (b-1)! / (a+b-1)! for positive integers
///
/// # Examples
///
/// ```
/// use axiomath::special::beta_function;
///
/// // B(1,1) = 1
/// assert!((beta_function(1.0, 1.0) - 1.0).abs() < 1e-10);
///
/// // B(2,3) = 1/12
/// assert!((beta_function(2.0, 3.0) - 1.0/12.0).abs() < 1e-10);
///
/// // Symmetry: B(a,b) = B(b,a)
/// assert!((beta_function(2.0, 3.0) - beta_function(3.0, 2.0)).abs() < 1e-10);
/// ```
///
/// # Complexity
///
/// Time: O(1) - Uses gamma function
/// Space: O(1)
pub fn beta_function<T: Float>(a: T, b: T) -> T {
    // Handle special cases
    if a.is_nan() || b.is_nan() || a <= T::zero() || b <= T::zero() {
        return T::nan();
    }
    if a.is_infinite() || b.is_infinite() {
        return T::zero();
    }

    // Use log-gamma for numerical stability
    let log_beta = log_beta_function(a, b);
    log_beta.exp()
}

/// Computes the natural logarithm of the beta function ln(B(a,b)).
///
/// This function is more numerically stable than computing log(beta(a,b))
/// when B(a,b) would overflow or underflow.
///
/// Uses the identity:
/// ```text
/// ln(B(a,b)) = ln(Γ(a)) + ln(Γ(b)) - ln(Γ(a+b))
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::special::{beta_function, log_beta_function};
///
/// let a = 5.0;
/// let b = 3.0;
///
/// // ln(B(a,b)) should equal ln(B(a,b))
/// let log_beta = log_beta_function(a, b);
/// let beta_val = beta_function(a, b);
/// assert!((log_beta - beta_val.ln()).abs() < 1e-10);
/// ```
///
/// # Complexity
///
/// Time: O(1) - Uses log-gamma function
/// Space: O(1)
pub fn log_beta_function<T: Float>(a: T, b: T) -> T {
    // Handle special cases
    if a.is_nan() || b.is_nan() || a <= T::zero() || b <= T::zero() {
        return T::nan();
    }

    // ln(B(a,b)) = ln(Γ(a)) + ln(Γ(b)) - ln(Γ(a+b))
    log_gamma_function(a) + log_gamma_function(b) - log_gamma_function(a + b)
}

/// Computes the incomplete beta function B(x; a, b).
///
/// The incomplete beta function is defined as:
/// ```text
/// B(x; a, b) = ∫₀ˣ t^(a-1) (1-t)^(b-1) dt    for 0 ≤ x ≤ 1
/// ```
///
/// It satisfies: B(x; a, b) + B(1-x; b, a) = B(a, b)
///
/// # Algorithm
///
/// Uses continued fraction representation with modified Lentz's method.
///
/// # Domain
///
/// - a, b > 0
/// - 0 ≤ x ≤ 1
///
/// # Examples
///
/// ```
/// use axiomath::special::incomplete_beta;
///
/// let a = 2.0;
/// let b = 3.0;
/// let x = 0.5;
///
/// let result = incomplete_beta(x, a, b);
/// assert!(result > 0.0 && result.is_finite());
/// ```
///
/// # Complexity
///
/// Time: O(1) - Fixed number of continued fraction terms
/// Space: O(1)
pub fn incomplete_beta<T: Float>(x: T, a: T, b: T) -> T {
    // Handle special cases
    if a.is_nan() || b.is_nan() || x.is_nan() {
        return T::nan();
    }
    if a <= T::zero() || b <= T::zero() {
        return T::nan();
    }
    if x < T::zero() || x > T::one() {
        return T::nan();
    }
    if x.is_zero() {
        return T::zero();
    }
    if x == T::one() {
        return beta_function(a, b);
    }

    // Use symmetry relation for better convergence
    // If x > (a+1)/(a+b+2), use B(x;a,b) = B(a,b) - B(1-x;b,a)
    let symmetry_point = (a + T::one()) / (a + b + T::from(2.0).unwrap());
    if x > symmetry_point {
        let beta_ab = beta_function(a, b);
        let inc_beta_complement = incomplete_beta(T::one() - x, b, a);
        return beta_ab - inc_beta_complement;
    }

    // Compute using continued fraction
    let bt = x.powf(a) * (T::one() - x).powf(b) * beta_function(a, b);
    incomplete_beta_continued_fraction(x, a, b) * bt
}

/// Continued fraction for incomplete beta function.
///
/// Computes I_x(a,b) / B(a,b) using modified Lentz's method.
fn incomplete_beta_continued_fraction<T: Float>(x: T, a: T, b: T) -> T {
    let max_iterations = 200;
    let epsilon = T::epsilon() * T::from(100.0).unwrap();

    let qab = a + b;
    let qap = a + T::one();
    let qam = a - T::one();

    // First step
    let mut c = T::one();
    let mut d = T::one() - qab * x / qap;

    if d.abs() < epsilon {
        d = epsilon;
    }
    d = T::one() / d;
    let mut h = d;

    let two = T::from(2.0).unwrap();

    for m in 1..=max_iterations {
        let m_float = T::from(m).unwrap();
        let m2 = two * m_float;

        // Even step
        let aa = m_float * (b - m_float) * x / ((qam + m2) * (a + m2));
        d = T::one() + aa * d;

        if d.abs() < epsilon {
            d = epsilon;
        }
        c = T::one() + aa / c;
        if c.abs() < epsilon {
            c = epsilon;
        }

        d = T::one() / d;
        h = h * d * c;

        // Odd step
        let aa = -(a + m_float) * (qab + m_float) * x / ((a + m2) * (qap + m2));
        d = T::one() + aa * d;

        if d.abs() < epsilon {
            d = epsilon;
        }
        c = T::one() + aa / c;
        if c.abs() < epsilon {
            c = epsilon;
        }

        d = T::one() / d;
        let delta = d * c;
        h = h * delta;

        // Check convergence
        if (delta - T::one()).abs() < epsilon {
            break;
        }
    }

    T::one() / (a * h)
}

/// Computes the regularized incomplete beta function I_x(a, b) = B(x; a, b) / B(a, b).
///
/// This is the cumulative distribution function (CDF) of the beta distribution.
///
/// # Mathematical Properties
///
/// - I_0(a, b) = 0
/// - I_1(a, b) = 1
/// - 0 ≤ I_x(a, b) ≤ 1 for all 0 ≤ x ≤ 1
/// - I_x(a, b) + I_{1-x}(b, a) = 1 (symmetry relation)
///
/// # Examples
///
/// ```
/// use axiomath::special::regularized_incomplete_beta;
///
/// let a = 2.0;
/// let b = 3.0;
/// let x = 0.5;
///
/// let result = regularized_incomplete_beta(x, a, b);
///
/// // Should be between 0 and 1
/// assert!(result >= 0.0 && result <= 1.0);
///
/// // I_0(a,b) = 0
/// assert!((regularized_incomplete_beta(0.0, a, b) - 0.0).abs() < 1e-10);
///
/// // I_1(a,b) = 1
/// assert!((regularized_incomplete_beta(1.0, a, b) - 1.0).abs() < 1e-10);
/// ```
///
/// # Complexity
///
/// Time: O(1) - Uses continued fraction
/// Space: O(1)
pub fn regularized_incomplete_beta<T: Float>(x: T, a: T, b: T) -> T {
    // Handle special cases
    if a.is_nan() || b.is_nan() || x.is_nan() {
        return T::nan();
    }
    if a <= T::zero() || b <= T::zero() {
        return T::nan();
    }
    if x < T::zero() || x > T::one() {
        return T::nan();
    }
    if x.is_zero() {
        return T::zero();
    }
    if x == T::one() {
        return T::one();
    }

    // Use symmetry for better convergence
    let symmetry_point = (a + T::one()) / (a + b + T::from(2.0).unwrap());
    if x > symmetry_point {
        return T::one() - regularized_incomplete_beta(T::one() - x, b, a);
    }

    // Compute using continued fraction
    let bt = x.powf(a) * (T::one() - x).powf(b) / beta_function(a, b);
    incomplete_beta_continued_fraction(x, a, b) * bt
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_beta_special_values() {
        // B(1,1) = 1
        assert_abs_diff_eq!(beta_function(1.0, 1.0), 1.0, epsilon = 1e-10);

        // B(a,1) = 1/a
        assert_abs_diff_eq!(beta_function(2.0, 1.0), 0.5, epsilon = 1e-10);
        assert_abs_diff_eq!(beta_function(5.0, 1.0), 0.2, epsilon = 1e-10);

        // B(1,b) = 1/b
        assert_abs_diff_eq!(beta_function(1.0, 3.0), 1.0/3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_beta_symmetry() {
        // B(a,b) = B(b,a)
        let test_cases = [
            (1.0, 2.0),
            (2.0, 3.0),
            (3.5, 2.5),
            (0.5, 1.5),
        ];

        for (a, b) in test_cases {
            let b_ab = beta_function(a, b);
            let b_ba = beta_function(b, a);
            assert_abs_diff_eq!(b_ab, b_ba, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_beta_integers() {
        // B(m,n) = (m-1)!(n-1)! / (m+n-1)! for positive integers
        // B(2,3) = 1!2! / 4! = 2/24 = 1/12
        assert_abs_diff_eq!(beta_function(2.0, 3.0), 1.0/12.0, epsilon = 1e-10);

        // B(3,2) = 2!1! / 4! = 2/24 = 1/12
        assert_abs_diff_eq!(beta_function(3.0, 2.0), 1.0/12.0, epsilon = 1e-10);

        // B(3,3) = 2!2! / 5! = 4/120 = 1/30
        assert_abs_diff_eq!(beta_function(3.0, 3.0), 1.0/30.0, epsilon = 1e-10);
    }

    #[test]
    fn test_beta_gamma_relation() {
        // B(a,b) = Γ(a)Γ(b) / Γ(a+b)
        let test_cases = [
            (1.5, 2.5),
            (2.0, 3.0),
            (0.5, 0.5),
        ];

        for (a, b) in test_cases {
            let beta_val = beta_function(a, b);
            let gamma_a = gamma_function(a);
            let gamma_b = gamma_function(b);
            let gamma_ab = gamma_function(a + b);
            let expected = gamma_a * gamma_b / gamma_ab;
            assert_abs_diff_eq!(beta_val, expected, epsilon = 1e-9);
        }
    }

    #[test]
    fn test_log_beta_consistency() {
        let test_cases = [
            (1.0, 1.0),
            (2.0, 3.0),
            (5.0, 7.0),
        ];

        for (a, b) in test_cases {
            let beta_val = beta_function(a, b);
            let log_beta_val = log_beta_function(a, b);
            assert_abs_diff_eq!(log_beta_val, NumFloat::ln(beta_val), epsilon = 1e-9);
        }
    }

    #[test]
    fn test_log_beta_large_values() {
        // For large a, b, beta overflows but log_beta should work
        let a = 50.0;
        let b = 50.0;
        let log_beta = log_beta_function(a, b);
        assert!(NumFloat::is_finite(log_beta));
        assert!(log_beta < 0.0); // B(a,b) < 1 for large a,b
    }

    #[test]
    fn test_incomplete_beta_limits() {
        let a = 2.0;
        let b = 3.0;

        // B(0; a, b) = 0
        assert_abs_diff_eq!(incomplete_beta(0.0, a, b), 0.0, epsilon = 1e-10);

        // B(1; a, b) = B(a, b)
        let beta_ab = beta_function(a, b);
        assert_abs_diff_eq!(incomplete_beta(1.0, a, b), beta_ab, epsilon = 1e-9);
    }

    #[test]
    fn test_incomplete_beta_symmetry() {
        // B(x;a,b) + B(1-x;b,a) = B(a,b)
        let a = 2.0;
        let b = 3.0;
        let x = 0.6;

        let inc_beta_x = incomplete_beta(x, a, b);
        let inc_beta_complement = incomplete_beta(1.0 - x, b, a);
        let beta_ab = beta_function(a, b);

        assert_abs_diff_eq!(inc_beta_x + inc_beta_complement, beta_ab, epsilon = 1e-8);
    }

    #[test]
    fn test_regularized_incomplete_beta_limits() {
        let a = 2.0;
        let b = 3.0;

        // I_0(a,b) = 0
        assert_abs_diff_eq!(regularized_incomplete_beta(0.0, a, b), 0.0, epsilon = 1e-10);

        // I_1(a,b) = 1
        assert_abs_diff_eq!(regularized_incomplete_beta(1.0, a, b), 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_regularized_incomplete_beta_range() {
        let a = 2.0;
        let b = 3.0;
        let x_values = [0.1, 0.25, 0.5, 0.75, 0.9];

        for x in x_values {
            let result = regularized_incomplete_beta(x, a, b);
            assert!(result >= 0.0 && result <= 1.0);
        }
    }

    #[test]
    fn test_regularized_incomplete_beta_monotonic() {
        // I_x(a,b) should be monotonically increasing in x
        let a = 2.0;
        let b = 3.0;
        let x_values = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0];

        let mut prev = -1e-10; // Small negative tolerance for numerical errors
        for x in x_values {
            let result = regularized_incomplete_beta(x, a, b);
            assert!(result >= prev - 1e-10, "result {} should be >= prev {}", result, prev);
            prev = result;
        }
    }

    #[test]
    fn test_regularized_incomplete_beta_symmetry() {
        // I_x(a,b) + I_{1-x}(b,a) = 1
        let a = 2.5;
        let b = 3.5;
        let x = 0.7;

        let i_x = regularized_incomplete_beta(x, a, b);
        let i_complement = regularized_incomplete_beta(1.0 - x, b, a);

        assert_abs_diff_eq!(i_x + i_complement, 1.0, epsilon = 1e-9);
    }

    #[test]
    fn test_beta_special_cases() {
        assert!(NumFloat::is_nan(beta_function(f64::NAN, 1.0)));
        assert!(NumFloat::is_nan(beta_function(1.0, f64::NAN)));
        assert!(NumFloat::is_nan(beta_function(-1.0, 1.0)));
        assert!(NumFloat::is_nan(beta_function(1.0, -1.0)));
        assert_eq!(beta_function(f64::INFINITY, 1.0), 0.0);
    }

    #[test]
    fn test_beta_f32() {
        let result: f32 = beta_function(2.0_f32, 3.0_f32);
        assert_abs_diff_eq!(result, 1.0_f32/12.0_f32, epsilon = 1e-6);
    }
}
