//! Vector angle operations.
//!
//! This module provides functions for computing angles between vectors and
//! testing for parallel and perpendicular relationships.
//!
//! # Mathematical Background
//!
//! The angle between two vectors can be computed using the dot product:
//! ```text
//! cos(θ) = (a · b) / (‖a‖ ‖b‖)
//! ```

use crate::core::traits::Float;
use crate::vector::{dot_product, magnitude, Vector};

/// Computes the angle between two vectors in radians.
///
/// # Mathematical Definition
///
/// ```text
/// θ = arccos((a · b) / (‖a‖ ‖b‖))
/// ```
///
/// The result is in the range [0, π] radians.
///
/// # Returns
///
/// `None` if either vector is zero (angle is undefined).
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, angle_between};
/// use std::f64::consts::PI;
///
/// // 90 degree angle
/// let i = vec3(1.0, 0.0, 0.0);
/// let j = vec3(0.0, 1.0, 0.0);
/// let angle = angle_between(&i, &j).unwrap();
/// assert!((angle - PI / 2.0).abs() < 1e-10);
///
/// // 180 degree angle
/// let a = vec3(1.0, 0.0, 0.0);
/// let b = vec3(-1.0, 0.0, 0.0);
/// let angle = angle_between(&a, &b).unwrap();
/// assert!((angle - PI).abs() < 1e-10);
/// ```
#[inline]
pub fn angle_between<T, const N: usize>(a: &Vector<T, N>, b: &Vector<T, N>) -> Option<T>
where
    T: Float,
{
    let mag_a = magnitude(a);
    let mag_b = magnitude(b);

    // Check if either vector is zero
    let epsilon = T::epsilon() * T::from(10).unwrap();
    if mag_a <= epsilon || mag_b <= epsilon {
        return None;
    }

    // Compute cos(θ) = (a · b) / (‖a‖ ‖b‖)
    let dot = dot_product(a, b);
    let cos_theta = dot / (mag_a * mag_b);

    // Clamp to [-1, 1] to handle floating-point errors
    let cos_theta = cos_theta.max(-T::one()).min(T::one());

    // Use Float trait's acos instead of our scalar module to avoid stack overflow
    Some(cos_theta.acos())
}

/// Checks if two vectors are parallel (or antiparallel).
///
/// Two vectors are parallel if the angle between them is 0 or π,
/// which means one is a scalar multiple of the other.
///
/// # Mathematical Definition
///
/// Vectors are parallel if:
/// ```text
/// |cos(θ)| ≈ 1  (within tolerance)
/// ```
///
/// # Arguments
///
/// * `a` - First vector
/// * `b` - Second vector
/// * `tolerance` - Angular tolerance in radians
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, is_parallel};
///
/// let a = vec3(1.0, 2.0, 3.0);
/// let b = vec3(2.0, 4.0, 6.0);  // b = 2a
/// assert!(is_parallel(&a, &b, 1e-10));
///
/// let c = vec3(-1.0, -2.0, -3.0);  // c = -a (antiparallel)
/// assert!(is_parallel(&a, &c, 1e-10));
///
/// let d = vec3(1.0, 0.0, 0.0);
/// assert!(!is_parallel(&a, &d, 1e-10));
/// ```
#[inline]
pub fn is_parallel<T, const N: usize>(a: &Vector<T, N>, b: &Vector<T, N>, tolerance: T) -> bool
where
    T: Float,
{
    match angle_between(a, b) {
        None => false, // Zero vectors are not parallel
        Some(angle) => {
            // Parallel if angle ≈ 0 or angle ≈ π
            let pi = T::from(std::f64::consts::PI).unwrap();
            angle.abs() <= tolerance || (angle - pi).abs() <= tolerance
        }
    }
}

/// Checks if two vectors are perpendicular (orthogonal).
///
/// Two vectors are perpendicular if the angle between them is π/2,
/// which means their dot product is zero.
///
/// # Mathematical Definition
///
/// Vectors are perpendicular if:
/// ```text
/// a · b ≈ 0  (within tolerance)
/// ```
///
/// This is more numerically stable than checking if the angle equals π/2.
///
/// # Arguments
///
/// * `a` - First vector
/// * `b` - Second vector
/// * `tolerance` - Dot product tolerance
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, is_perpendicular};
///
/// let i = vec3(1.0, 0.0, 0.0);
/// let j = vec3(0.0, 1.0, 0.0);
/// assert!(is_perpendicular(&i, &j, 1e-10));
///
/// let a = vec3(1.0, 2.0, 3.0);
/// let b = vec3(2.0, 4.0, 6.0);
/// assert!(!is_perpendicular(&a, &b, 1e-10));
/// ```
#[inline]
pub fn is_perpendicular<T, const N: usize>(
    a: &Vector<T, N>,
    b: &Vector<T, N>,
    tolerance: T,
) -> bool
where
    T: Float,
{
    let dot = dot_product(a, b);
    dot.abs() <= tolerance
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vector::vec3;
    use std::f64::consts::PI;

    #[test]
    fn test_angle_between_perpendicular() {
        let i = vec3(1.0, 0.0, 0.0);
        let j = vec3(0.0, 1.0, 0.0);
        let angle = angle_between(&i, &j).unwrap();

        assert!((angle - PI / 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_angle_between_parallel() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(2.0, 4.0, 6.0);
        let angle: f64 = angle_between(&a, &b).unwrap();

        assert!(angle.abs() < 1e-10); // Should be 0
    }

    #[test]
    fn test_angle_between_antiparallel() {
        let a = vec3(1.0, 0.0, 0.0);
        let b = vec3(-1.0, 0.0, 0.0);
        let angle = angle_between(&a, &b).unwrap();

        assert!((angle - PI).abs() < 1e-10);
    }

    #[test]
    fn test_angle_between_same_vector() {
        let v = vec3(1.0, 2.0, 3.0);
        let angle: f64 = angle_between(&v, &v).unwrap();

        assert!(angle.abs() < 1e-10); // Should be 0
    }

    #[test]
    fn test_angle_between_zero_vector() {
        let v = vec3(1.0, 2.0, 3.0);
        let zero = vec3(0.0, 0.0, 0.0);

        assert!(angle_between(&v, &zero).is_none());
        assert!(angle_between(&zero, &v).is_none());
    }

    #[test]
    fn test_angle_between_45_degrees() {
        let a = vec3(1.0, 0.0, 0.0);
        let b = vec3(1.0, 1.0, 0.0);
        let angle: f64 = angle_between(&a, &b).unwrap();

        assert!((angle - PI / 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_angle_between_120_degrees() {
        let a = vec3(1.0, 0.0, 0.0);
        let b = vec3(-0.5, 0.866025403784, 0.0); // 120 degrees
        let angle: f64 = angle_between(&a, &b).unwrap();

        assert!((angle - 2.0 * PI / 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_is_parallel_same_direction() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(2.0, 4.0, 6.0);

        assert!(is_parallel(&a, &b, 1e-10));
    }

    #[test]
    fn test_is_parallel_opposite_direction() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(-1.0, -2.0, -3.0);

        assert!(is_parallel(&a, &b, 1e-10));
    }

    #[test]
    fn test_is_parallel_same_vector() {
        let v = vec3(1.0, 2.0, 3.0);

        assert!(is_parallel(&v, &v, 1e-10));
    }

    #[test]
    fn test_is_not_parallel() {
        let a = vec3(1.0, 0.0, 0.0);
        let b = vec3(0.0, 1.0, 0.0);

        assert!(!is_parallel(&a, &b, 1e-10));
    }

    #[test]
    fn test_is_parallel_with_tolerance() {
        let a = vec3(1.0, 0.0, 0.0);
        let b = vec3(1.0, 0.01, 0.0); // Nearly parallel

        // Should not be parallel with tight tolerance
        assert!(!is_parallel(&a, &b, 1e-10));

        // Should be parallel with loose tolerance
        assert!(is_parallel(&a, &b, 0.1));
    }

    #[test]
    fn test_is_perpendicular_orthogonal() {
        let i = vec3(1.0, 0.0, 0.0);
        let j = vec3(0.0, 1.0, 0.0);

        assert!(is_perpendicular(&i, &j, 1e-10));
    }

    #[test]
    fn test_is_perpendicular_3d() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(2.0, -1.0, 0.0); // Perpendicular: 1*2 + 2*(-1) + 3*0 = 0

        assert!(is_perpendicular(&a, &b, 1e-10));
    }

    #[test]
    fn test_is_not_perpendicular() {
        let a = vec3(1.0, 0.0, 0.0);
        let b = vec3(1.0, 1.0, 0.0);

        assert!(!is_perpendicular(&a, &b, 1e-10));
    }

    #[test]
    fn test_is_perpendicular_with_tolerance() {
        let a = vec3(1.0, 0.0, 0.0);
        let b = vec3(0.0, 1.0, 0.001); // Nearly perpendicular

        // Dot product = 0.0, so should be perpendicular
        assert!(is_perpendicular(&a, &b, 1e-2));
    }

    #[test]
    fn test_perpendicular_parallel_mutually_exclusive() {
        let a = vec3(1.0, 2.0, 3.0);
        let b_parallel = vec3(2.0, 4.0, 6.0);
        let b_perp = vec3(2.0, -1.0, 0.0);
        let b_neither = vec3(1.0, 0.0, 0.0);

        // Parallel vectors are not perpendicular
        assert!(is_parallel(&a, &b_parallel, 1e-10));
        assert!(!is_perpendicular(&a, &b_parallel, 1e-10));

        // Perpendicular vectors are not parallel
        assert!(!is_parallel(&a, &b_perp, 1e-10));
        assert!(is_perpendicular(&a, &b_perp, 1e-10));

        // Some vectors are neither
        assert!(!is_parallel(&a, &b_neither, 1e-10));
        assert!(!is_perpendicular(&a, &b_neither, 1e-10));
    }

    #[test]
    fn test_angle_symmetry() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(4.0, 5.0, 6.0);

        // angle(a, b) = angle(b, a)
        let angle_ab: f64 = angle_between(&a, &b).unwrap();
        let angle_ba: f64 = angle_between(&b, &a).unwrap();

        assert!((angle_ab - angle_ba).abs() < 1e-10);
    }

    #[test]
    fn test_parallel_symmetry() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(2.0, 4.0, 6.0);

        // is_parallel(a, b) = is_parallel(b, a)
        assert_eq!(is_parallel(&a, &b, 1e-10), is_parallel(&b, &a, 1e-10));
    }

    #[test]
    fn test_perpendicular_symmetry() {
        let a = vec3(1.0, 0.0, 0.0);
        let b = vec3(0.0, 1.0, 0.0);

        // is_perpendicular(a, b) = is_perpendicular(b, a)
        assert_eq!(
            is_perpendicular(&a, &b, 1e-10),
            is_perpendicular(&b, &a, 1e-10)
        );
    }

    #[test]
    fn test_angle_range() {
        let a = vec3(1.0, 0.0, 0.0);

        // Test various angles
        for i in 0..=12 {
            let theta = (i as f64) * PI / 12.0;
            let b = vec3(theta.cos(), theta.sin(), 0.0);
            let angle: f64 = angle_between(&a, &b).unwrap();

            // Angle should be in [0, π]
            assert!(angle >= 0.0);
            assert!(angle <= PI + 1e-10);
        }
    }
}
