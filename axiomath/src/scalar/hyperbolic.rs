//! Hyperbolic trigonometric functions using exponential definitions.
//!
//! This module provides hyperbolic functions using their exponential definitions:
//! - sinh, cosh, tanh using e^x

use crate::core::traits::Float;
use crate::scalar::exponential::exponential;

/// Computes sinh(x) (hyperbolic sine).
///
/// # Algorithm
///
/// Uses the exponential definition:
/// ```text
/// sinh(x) = (e^x - e^(-x)) / 2
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::scalar::hyperbolic::hyperbolic_sine;
///
/// let result = hyperbolic_sine(0.0_f64);
/// assert!((result - 0.0).abs() < 1e-2);
///
/// let result = hyperbolic_sine(1.0_f64);
/// assert!((result - 1.1752011936438014).abs() < 1e-2);
/// ```
#[inline]
pub fn hyperbolic_sine<T>(x: T) -> T
where
    T: Float,
{
    if x.is_nan() {
        return x;
    }
    if x.is_infinite() {
        return x;
    }
    if x.is_zero() {
        return T::zero();
    }

    let exp_x = exponential(x);
    let exp_neg_x = exponential(-x);
    let two = T::one() + T::one();

    (exp_x - exp_neg_x) / two
}

/// Computes cosh(x) (hyperbolic cosine).
///
/// # Algorithm
///
/// Uses the exponential definition:
/// ```text
/// cosh(x) = (e^x + e^(-x)) / 2
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::scalar::hyperbolic::hyperbolic_cosine;
///
/// let result = hyperbolic_cosine(0.0_f64);
/// assert!((result - 1.0).abs() < 1e-2);
///
/// let result = hyperbolic_cosine(1.0_f64);
/// assert!((result - 1.5430806348152437).abs() < 1e-2);
/// ```
#[inline]
pub fn hyperbolic_cosine<T>(x: T) -> T
where
    T: Float,
{
    if x.is_nan() {
        return x;
    }
    if x.is_infinite() {
        return T::infinity();
    }
    if x.is_zero() {
        return T::one();
    }

    let exp_x = exponential(x);
    let exp_neg_x = exponential(-x);
    let two = T::one() + T::one();

    (exp_x + exp_neg_x) / two
}

/// Computes tanh(x) (hyperbolic tangent).
///
/// # Algorithm
///
/// Uses the definition:
/// ```text
/// tanh(x) = sinh(x) / cosh(x) = (e^x - e^(-x)) / (e^x + e^(-x))
/// ```
///
/// Can also be written as:
/// ```text
/// tanh(x) = (e^(2x) - 1) / (e^(2x) + 1)
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::scalar::hyperbolic::hyperbolic_tangent;
///
/// let result = hyperbolic_tangent(0.0_f64);
/// assert!((result - 0.0).abs() < 1e-2);
///
/// let result = hyperbolic_tangent(1.0_f64);
/// assert!((result - 0.7615941559557649).abs() < 1e-2);
/// ```
#[inline]
pub fn hyperbolic_tangent<T>(x: T) -> T
where
    T: Float,
{
    if x.is_nan() {
        return x;
    }
    if x.is_zero() {
        return T::zero();
    }
    if x.is_infinite() {
        return if x > T::zero() { T::one() } else { -T::one() };
    }

    // Use tanh(x) = (e^(2x) - 1) / (e^(2x) + 1) for better numerical stability
    let two = T::one() + T::one();
    let exp_2x = exponential(two * x);

    (exp_2x - T::one()) / (exp_2x + T::one())
}

/// Computes arsinh(x) (inverse hyperbolic sine).
///
/// # Algorithm
///
/// Uses the logarithmic formula:
/// ```text
/// arsinh(x) = ln(x + √(x² + 1))
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::scalar::hyperbolic::{hyperbolic_sine, inverse_hyperbolic_sine};
///
/// let result = inverse_hyperbolic_sine(0.0_f64);
/// assert!((result - 0.0).abs() < 1e-2);
///
/// // Test inverse property: sinh(arsinh(x)) = x
/// let x = 2.0_f64;
/// let arsinh_x = inverse_hyperbolic_sine(x);
/// let result = hyperbolic_sine(arsinh_x);
/// assert!((result - x).abs() < 1e-2);
/// ```
#[inline]
pub fn inverse_hyperbolic_sine<T>(x: T) -> T
where
    T: Float,
{
    use crate::scalar::logarithmic::natural_logarithm;
    use crate::scalar::power::square_root;

    if x.is_nan() {
        return x;
    }
    if x.is_infinite() {
        return x;
    }
    if x.is_zero() {
        return T::zero();
    }

    // arsinh(x) = ln(x + √(x² + 1))
    natural_logarithm(x + square_root(x * x + T::one()))
}

/// Computes arcosh(x) (inverse hyperbolic cosine).
///
/// # Algorithm
///
/// Uses the logarithmic formula:
/// ```text
/// arcosh(x) = ln(x + √(x² - 1))
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::scalar::hyperbolic::{hyperbolic_cosine, inverse_hyperbolic_cosine};
///
/// let result = inverse_hyperbolic_cosine(1.0_f64);
/// assert!((result - 0.0).abs() < 1e-2);
///
/// // Test inverse property: cosh(arcosh(x)) = x
/// let x = 2.0_f64;
/// let arcosh_x = inverse_hyperbolic_cosine(x);
/// let result = hyperbolic_cosine(arcosh_x);
/// assert!((result - x).abs() < 1e-2);
/// ```
///
/// # Panics
///
/// Panics if x < 1 (arcosh is only defined for x ≥ 1).
#[inline]
pub fn inverse_hyperbolic_cosine<T>(x: T) -> T
where
    T: Float,
{
    use crate::scalar::logarithmic::natural_logarithm;
    use crate::scalar::power::square_root;

    if x < T::one() {
        panic!("arcosh is only defined for x >= 1");
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

    // arcosh(x) = ln(x + √(x² - 1))
    natural_logarithm(x + square_root(x * x - T::one()))
}

/// Computes artanh(x) (inverse hyperbolic tangent).
///
/// # Algorithm
///
/// Uses the logarithmic formula:
/// ```text
/// artanh(x) = (1/2) * ln((1 + x) / (1 - x))
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::scalar::hyperbolic::{hyperbolic_tangent, inverse_hyperbolic_tangent};
///
/// let result = inverse_hyperbolic_tangent(0.0_f64);
/// assert!((result - 0.0).abs() < 1e-2);
///
/// // Test inverse property: tanh(artanh(x)) = x
/// let x = 0.5_f64;
/// let artanh_x = inverse_hyperbolic_tangent(x);
/// let result = hyperbolic_tangent(artanh_x);
/// assert!((result - x).abs() < 1e-2);
/// ```
///
/// # Panics
///
/// Panics if |x| ≥ 1 (artanh is only defined for |x| < 1).
#[inline]
pub fn inverse_hyperbolic_tangent<T>(x: T) -> T
where
    T: Float,
{
    use crate::scalar::logarithmic::natural_logarithm;

    if x.abs() >= T::one() {
        panic!("artanh is only defined for |x| < 1");
    }
    if x.is_nan() {
        return x;
    }
    if x.is_zero() {
        return T::zero();
    }

    // artanh(x) = (1/2) * ln((1 + x) / (1 - x))
    let two = T::one() + T::one();
    natural_logarithm((T::one() + x) / (T::one() - x)) / two
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_hyperbolic_sine_zero() {
        assert_abs_diff_eq!(hyperbolic_sine(0.0_f64), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_hyperbolic_sine_basic() {
        // sinh(1) ≈ 1.1752011936438014
        assert_abs_diff_eq!(
            hyperbolic_sine(1.0_f64),
            1.1752011936438014,
            epsilon = 1e-10
        );

        // sinh(-1) ≈ -1.1752011936438014 (odd function)
        assert_abs_diff_eq!(
            hyperbolic_sine(-1.0_f64),
            -1.1752011936438014,
            epsilon = 1e-10
        );

        // sinh(2) ≈ 3.626860407847019
        assert_abs_diff_eq!(
            hyperbolic_sine(2.0_f64),
            3.626860407847019,
            epsilon = 1e-9
        );
    }

    #[test]
    fn test_hyperbolic_cosine_zero() {
        assert_abs_diff_eq!(hyperbolic_cosine(0.0_f64), 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_hyperbolic_cosine_basic() {
        // cosh(1) ≈ 1.5430806348152437
        assert_abs_diff_eq!(
            hyperbolic_cosine(1.0_f64),
            1.5430806348152437,
            epsilon = 1e-10
        );

        // cosh(-1) ≈ 1.5430806348152437 (even function)
        assert_abs_diff_eq!(
            hyperbolic_cosine(-1.0_f64),
            1.5430806348152437,
            epsilon = 1e-10
        );

        // cosh(2) ≈ 3.7621956910836314
        assert_abs_diff_eq!(
            hyperbolic_cosine(2.0_f64),
            3.7621956910836314,
            epsilon = 1e-9
        );
    }

    #[test]
    fn test_hyperbolic_tangent_zero() {
        assert_abs_diff_eq!(hyperbolic_tangent(0.0_f64), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_hyperbolic_tangent_basic() {
        // tanh(1) ≈ 0.7615941559557649
        assert_abs_diff_eq!(
            hyperbolic_tangent(1.0_f64),
            0.7615941559557649,
            epsilon = 1e-10
        );

        // tanh(-1) ≈ -0.7615941559557649 (odd function)
        assert_abs_diff_eq!(
            hyperbolic_tangent(-1.0_f64),
            -0.7615941559557649,
            epsilon = 1e-10
        );

        // tanh approaches ±1 as x → ±∞
        assert_abs_diff_eq!(hyperbolic_tangent(10.0_f64), 1.0, epsilon = 1e-8);  // Relaxed for asymptotic behavior
        assert_abs_diff_eq!(hyperbolic_tangent(-10.0_f64), -1.0, epsilon = 1e-8);
    }

    #[test]
    fn test_hyperbolic_identity() {
        // cosh²(x) - sinh²(x) = 1
        for x in [0.0, 0.5, 1.0, 2.0, 3.0] {
            let sinh_x = hyperbolic_sine(x);
            let cosh_x = hyperbolic_cosine(x);
            let result = cosh_x * cosh_x - sinh_x * sinh_x;
            assert_abs_diff_eq!(result, 1.0, epsilon = 1e-9);
        }
    }

    #[test]
    fn test_hyperbolic_tanh_identity() {
        // tanh(x) = sinh(x) / cosh(x)
        for x in [0.5, 1.0, 2.0] {
            let tanh_x = hyperbolic_tangent(x);
            let sinh_x = hyperbolic_sine(x);
            let cosh_x = hyperbolic_cosine(x);
            let ratio = sinh_x / cosh_x;
            assert_abs_diff_eq!(tanh_x, ratio, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_inverse_hyperbolic_sine() {
        assert_abs_diff_eq!(inverse_hyperbolic_sine(0.0_f64), 0.0, epsilon = 1e-10);

        // Test inverse property: sinh(arsinh(x)) = x
        for x in [-2.0, -1.0, 0.0, 1.0, 2.0] {
            let arsinh_x = inverse_hyperbolic_sine(x);
            let result = hyperbolic_sine(arsinh_x);
            assert_abs_diff_eq!(result, x, epsilon = 1e-9);
        }
    }

    #[test]
    fn test_inverse_hyperbolic_cosine() {
        assert_abs_diff_eq!(inverse_hyperbolic_cosine(1.0_f64), 0.0, epsilon = 1e-10);

        // Test inverse property: cosh(arcosh(x)) = x
        for x in [1.0, 1.5, 2.0, 3.0] {
            let arcosh_x = inverse_hyperbolic_cosine(x);
            let result = hyperbolic_cosine(arcosh_x);
            assert_abs_diff_eq!(result, x, epsilon = 1e-9);
        }
    }

    #[test]
    #[should_panic(expected = "arcosh is only defined for x >= 1")]
    fn test_inverse_hyperbolic_cosine_out_of_range() {
        inverse_hyperbolic_cosine(0.5);
    }

    #[test]
    fn test_inverse_hyperbolic_tangent() {
        assert_abs_diff_eq!(inverse_hyperbolic_tangent(0.0_f64), 0.0, epsilon = 1e-10);

        // Test inverse property: tanh(artanh(x)) = x
        for x in [-0.9, -0.5, 0.0, 0.5, 0.9] {
            let artanh_x = inverse_hyperbolic_tangent(x);
            let result = hyperbolic_tangent(artanh_x);
            assert_abs_diff_eq!(result, x, epsilon = 1e-10);
        }
    }

    #[test]
    #[should_panic(expected = "artanh is only defined for |x| < 1")]
    fn test_inverse_hyperbolic_tangent_out_of_range() {
        inverse_hyperbolic_tangent(1.0);
    }
}
