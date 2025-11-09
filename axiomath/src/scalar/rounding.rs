//! Rounding and comparison functions.
//!
//! This module provides rounding and comparison operations:
//! - floor, ceil, round, truncate
//! - min, max, clamp

use crate::core::traits::Float;

/// Returns the largest integer less than or equal to x.
///
/// # Examples
///
/// ```
/// use axiomath::scalar::rounding::floor;
///
/// assert_eq!(floor(3.7), 3.0);
/// assert_eq!(floor(-3.7), -4.0);
/// assert_eq!(floor(5.0), 5.0);
/// ```
#[inline]
pub fn floor<T>(x: T) -> T
where
    T: Float,
{
    x.floor()
}

/// Returns the smallest integer greater than or equal to x.
///
/// # Examples
///
/// ```
/// use axiomath::scalar::rounding::ceiling;
///
/// assert_eq!(ceiling(3.2), 4.0);
/// assert_eq!(ceiling(-3.2), -3.0);
/// assert_eq!(ceiling(5.0), 5.0);
/// ```
#[inline]
pub fn ceiling<T>(x: T) -> T
where
    T: Float,
{
    x.ceil()
}

/// Rounds to the nearest integer. Halfway cases round away from zero.
///
/// # Examples
///
/// ```
/// use axiomath::scalar::rounding::round;
///
/// assert_eq!(round(3.4), 3.0);
/// assert_eq!(round(3.5), 4.0);
/// assert_eq!(round(-3.5), -4.0);
/// ```
#[inline]
pub fn round<T>(x: T) -> T
where
    T: Float,
{
    x.round()
}

/// Returns the integer part of x, truncating towards zero.
///
/// # Examples
///
/// ```
/// use axiomath::scalar::rounding::truncate;
///
/// assert_eq!(truncate(3.7), 3.0);
/// assert_eq!(truncate(-3.7), -3.0);
/// assert_eq!(truncate(5.0), 5.0);
/// ```
#[inline]
pub fn truncate<T>(x: T) -> T
where
    T: Float,
{
    x.trunc()
}

/// Returns the minimum of two values.
///
/// # Examples
///
/// ```
/// use axiomath::scalar::rounding::minimum;
///
/// assert_eq!(minimum(3.0, 5.0), 3.0);
/// assert_eq!(minimum(-1.0, 1.0), -1.0);
/// assert_eq!(minimum(2.5, 2.5), 2.5);
/// ```
#[inline]
pub fn minimum<T>(a: T, b: T) -> T
where
    T: Float,
{
    a.min(b)
}

/// Returns the maximum of two values.
///
/// # Examples
///
/// ```
/// use axiomath::scalar::rounding::maximum;
///
/// assert_eq!(maximum(3.0, 5.0), 5.0);
/// assert_eq!(maximum(-1.0, 1.0), 1.0);
/// assert_eq!(maximum(2.5, 2.5), 2.5);
/// ```
#[inline]
pub fn maximum<T>(a: T, b: T) -> T
where
    T: Float,
{
    a.max(b)
}

/// Clamps a value between a minimum and maximum.
///
/// Returns:
/// - `min` if `x < min`
/// - `max` if `x > max`
/// - `x` otherwise
///
/// # Examples
///
/// ```
/// use axiomath::scalar::rounding::clamp;
///
/// assert_eq!(clamp(5.0, 0.0, 10.0), 5.0);
/// assert_eq!(clamp(-5.0, 0.0, 10.0), 0.0);
/// assert_eq!(clamp(15.0, 0.0, 10.0), 10.0);
/// ```
///
/// # Panics
///
/// Panics if `min > max`.
#[inline]
pub fn clamp<T>(x: T, min: T, max: T) -> T
where
    T: Float,
{
    assert!(min <= max, "min must be less than or equal to max");
    x.clamp(min, max)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_floor() {
        assert_eq!(floor(3.7), 3.0);
        assert_eq!(floor(3.0), 3.0);
        assert_eq!(floor(-3.2), -4.0);
        assert_eq!(floor(-3.0), -3.0);
        assert_eq!(floor(0.0), 0.0);
    }

    #[test]
    fn test_ceiling() {
        assert_eq!(ceiling(3.2), 4.0);
        assert_eq!(ceiling(3.0), 3.0);
        assert_eq!(ceiling(-3.7), -3.0);
        assert_eq!(ceiling(-3.0), -3.0);
        assert_eq!(ceiling(0.0), 0.0);
    }

    #[test]
    fn test_round() {
        assert_eq!(round(3.4), 3.0);
        assert_eq!(round(3.5), 4.0);
        assert_eq!(round(3.6), 4.0);
        assert_eq!(round(-3.4), -3.0);
        assert_eq!(round(-3.5), -4.0);
        assert_eq!(round(0.0), 0.0);
    }

    #[test]
    fn test_truncate() {
        assert_eq!(truncate(3.7), 3.0);
        assert_eq!(truncate(3.2), 3.0);
        assert_eq!(truncate(-3.7), -3.0);
        assert_eq!(truncate(-3.2), -3.0);
        assert_eq!(truncate(0.0), 0.0);
        assert_eq!(truncate(5.0), 5.0);
    }

    #[test]
    fn test_minimum() {
        assert_eq!(minimum(3.0, 5.0), 3.0);
        assert_eq!(minimum(5.0, 3.0), 3.0);
        assert_eq!(minimum(-1.0, 1.0), -1.0);
        assert_eq!(minimum(2.5, 2.5), 2.5);
        assert_eq!(minimum(0.0, 0.0), 0.0);
    }

    #[test]
    fn test_maximum() {
        assert_eq!(maximum(3.0, 5.0), 5.0);
        assert_eq!(maximum(5.0, 3.0), 5.0);
        assert_eq!(maximum(-1.0, 1.0), 1.0);
        assert_eq!(maximum(2.5, 2.5), 2.5);
        assert_eq!(maximum(0.0, 0.0), 0.0);
    }

    #[test]
    fn test_clamp() {
        assert_eq!(clamp(5.0, 0.0, 10.0), 5.0);
        assert_eq!(clamp(-5.0, 0.0, 10.0), 0.0);
        assert_eq!(clamp(15.0, 0.0, 10.0), 10.0);
        assert_eq!(clamp(0.0, 0.0, 10.0), 0.0);
        assert_eq!(clamp(10.0, 0.0, 10.0), 10.0);
    }

    #[test]
    #[should_panic(expected = "min must be less than or equal to max")]
    fn test_clamp_invalid_range() {
        clamp(5.0, 10.0, 0.0);
    }
}
