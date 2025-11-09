//! Comparison utilities for floating-point numbers.
//!
//! Floating-point comparison requires special care due to rounding errors.
//! This module provides functions for approximate equality testing.

use crate::core::Float;
use num_traits::Zero as NumZero;

/// Compare two floating-point values for approximate equality.
///
/// This function uses both absolute and relative tolerance to determine
/// if two values are "close enough" to be considered equal.
///
/// # Mathematical Definition
///
/// Two values a and b are considered equal if:
/// |a - b| ≤ max(abs_tol, rel_tol × max(|a|, |b|))
///
/// # Arguments
///
/// * `a` - First value
/// * `b` - Second value
/// * `abs_tol` - Absolute tolerance
/// * `rel_tol` - Relative tolerance
///
/// # Examples
///
/// ```
/// use axiomath::utils::approximately_equal;
///
/// let a = 1.0 + 0.0001;
/// let b = 1.0;
///
/// assert!(approximately_equal(a, b, 1e-3, 1e-3));
/// assert!(!approximately_equal(a, b, 1e-5, 1e-5));
/// ```
pub fn approximately_equal<T: Float>(a: T, b: T, abs_tol: T, rel_tol: T) -> bool {
    let diff = (a - b).abs();
    let max_val = a.abs().max(b.abs());
    diff <= abs_tol.max(rel_tol * max_val)
}

/// Compare two floating-point values using absolute tolerance only.
///
/// # Examples
///
/// ```
/// use axiomath::utils::abs_diff_eq;
///
/// assert!(abs_diff_eq(1.0, 1.0001, 0.001));
/// assert!(!abs_diff_eq(1.0, 1.0001, 0.00001));
/// ```
#[inline]
pub fn abs_diff_eq<T: Float>(a: T, b: T, tolerance: T) -> bool {
    (a - b).abs() <= tolerance
}

/// Compare two floating-point values using relative tolerance only.
///
/// # Examples
///
/// ```
/// use axiomath::utils::rel_diff_eq;
///
/// assert!(rel_diff_eq(1000.0, 1001.0, 0.01));
/// assert!(!rel_diff_eq(1000.0, 1001.0, 0.0001));
/// ```
#[inline]
pub fn rel_diff_eq<T: Float>(a: T, b: T, tolerance: T) -> bool {
    let diff = (a - b).abs();
    let max_val = a.abs().max(b.abs());
    diff <= tolerance * max_val
}

/// Calculate the absolute error between an approximate and exact value.
///
/// # Mathematical Definition
///
/// Absolute error = |approximate - exact|
///
/// # Examples
///
/// ```
/// use axiomath::utils::absolute_error;
///
/// let approx = 3.14;
/// let exact = 3.14159265;
/// let error = absolute_error(approx, exact);
/// assert!(error < 0.002);
/// ```
#[inline]
pub fn absolute_error<T: Float>(approximate: T, exact: T) -> T {
    (approximate - exact).abs()
}

/// Calculate the relative error between an approximate and exact value.
///
/// # Mathematical Definition
///
/// Relative error = |approximate - exact| / |exact|
///
/// Returns NaN if exact is zero.
///
/// # Examples
///
/// ```
/// use axiomath::utils::relative_error;
/// use num_traits::Float;
///
/// let approx = 100.0;
/// let exact = 99.0;
/// let error = relative_error(approx, exact);
/// assert!(Float::abs(error - 0.0101) < 0.0001);
/// ```
#[inline]
pub fn relative_error<T: Float>(approximate: T, exact: T) -> T {
    if NumZero::is_zero(&exact) {
        return T::nan();
    }
    (approximate - exact).abs() / exact.abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_approximately_equal() {
        assert!(approximately_equal(1.0, 1.0, 0.0, 0.0));
        assert!(approximately_equal(1.0, 1.000001, 1e-5, 1e-5));
        assert!(!approximately_equal(1.0, 1.001, 1e-5, 1e-5));
    }

    #[test]
    fn test_absolute_error() {
        assert_eq!(absolute_error(3.0, 2.5), 0.5);
        assert_eq!(absolute_error(2.5, 3.0), 0.5);
        assert_eq!(absolute_error(1.0, 1.0), 0.0);
    }

    #[test]
    fn test_relative_error() {
        let error = relative_error(100.0_f64, 99.0_f64);
        assert!((error - (1.0_f64 / 99.0_f64)).abs() < 1e-10);
    }
}
