//! Vector projection and reflection operations.
//!
//! This module provides operations for projecting vectors onto other vectors,
//! computing rejection (perpendicular component), and reflecting vectors across
//! planes defined by normal vectors.
//!
//! # Mathematical Background
//!
//! - **Projection**: The component of one vector in the direction of another
//! - **Rejection**: The component perpendicular to the projection
//! - **Reflection**: Mirroring a vector across a plane

use crate::core::traits::Scalar;
use crate::vector::{dot_product, magnitude_squared, normalize, Vector};

/// Projects vector `v` onto vector `onto`.
///
/// Returns the component of `v` in the direction of `onto`.
///
/// # Mathematical Definition
///
/// ```text
/// proj_onto(v) = ((v · onto) / (onto · onto)) * onto
///              = ((v · onto) / ‖onto‖²) * onto
/// ```
///
/// # Returns
///
/// `None` if `onto` is the zero vector (undefined projection).
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, project};
///
/// let v = vec3(3.0, 4.0, 0.0);
/// let onto = vec3(1.0, 0.0, 0.0);
/// let proj = project(&v, &onto).unwrap();
///
/// assert_eq!(proj[0], 3.0);
/// assert_eq!(proj[1], 0.0);
/// assert_eq!(proj[2], 0.0);
/// ```
#[inline]
pub fn project<T, const N: usize>(v: &Vector<T, N>, onto: &Vector<T, N>) -> Option<Vector<T, N>>
where
    T: Scalar,
{
    let onto_mag_sq = magnitude_squared(onto);

    // Check if onto is zero vector
    if onto_mag_sq == T::zero() {
        return None;
    }

    let scalar = dot_product(v, onto) / onto_mag_sq;
    Some(onto.map(|x| x * scalar))
}

/// Computes the rejection of vector `v` from vector `from`.
///
/// Returns the component of `v` perpendicular to `from`.
/// This is equivalent to `v - project(v, from)`.
///
/// # Mathematical Definition
///
/// ```text
/// rej_from(v) = v - proj_from(v)
/// ```
///
/// The rejection is perpendicular to `from`: `rej · from = 0`
///
/// # Returns
///
/// `None` if `from` is the zero vector.
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, reject, dot_product};
///
/// let v = vec3(3.0, 4.0, 0.0);
/// let from = vec3(1.0, 0.0, 0.0);
/// let rej = reject(&v, &from).unwrap();
///
/// assert_eq!(rej[0], 0.0);
/// assert_eq!(rej[1], 4.0);
/// assert_eq!(rej[2], 0.0);
///
/// // Rejection is perpendicular to 'from'
/// assert_eq!(dot_product(&rej, &from), 0.0);
/// ```
#[inline]
pub fn reject<T, const N: usize>(v: &Vector<T, N>, from: &Vector<T, N>) -> Option<Vector<T, N>>
where
    T: Scalar,
{
    let proj = project(v, from)?;
    Some(v.zip_with(&proj, |vi, pi| vi - pi))
}

/// Reflects vector `v` across a plane defined by normal vector `normal`.
///
/// The normal vector defines a plane perpendicular to it, and the function
/// returns the mirror image of `v` across this plane.
///
/// # Mathematical Definition
///
/// ```text
/// reflect(v, n) = v - 2 * proj_n(v)
///               = v - 2 * ((v · n̂) * n̂)
/// ```
///
/// where n̂ is the normalized normal vector.
///
/// # Properties
///
/// - Reflecting twice returns the original vector: `reflect(reflect(v, n), n) = v`
/// - The reflected vector has the same magnitude: `‖reflect(v, n)‖ = ‖v‖`
///
/// # Returns
///
/// `None` if `normal` is the zero vector.
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, reflect};
///
/// let v = vec3(1.0, -1.0, 0.0);
/// let normal = vec3(0.0, 1.0, 0.0); // Reflect across xz-plane
/// let reflected = reflect(&v, &normal).unwrap();
///
/// assert_eq!(reflected[0], 1.0);
/// assert_eq!(reflected[1], 1.0);  // y-component flipped
/// assert_eq!(reflected[2], 0.0);
/// ```
#[inline]
pub fn reflect<T, const N: usize>(v: &Vector<T, N>, normal: &Vector<T, N>) -> Option<Vector<T, N>>
where
    T: crate::core::traits::Float,
{
    // Normalize the normal vector
    let n = normalize(normal)?;

    // Compute projection: (v · n̂) * n̂
    let proj_scalar = dot_product(v, &n);
    let two = T::one() + T::one();

    // reflect = v - 2 * proj_n(v)
    Some(v.zip_with(&n, |vi, ni| vi - two * proj_scalar * ni))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vector::{dot_product, magnitude, vec3};

    #[test]
    fn test_project() {
        let v = vec3(3.0, 4.0, 0.0);
        let onto = vec3(1.0, 0.0, 0.0);
        let proj = project(&v, &onto).unwrap();

        assert_eq!(proj[0], 3.0);
        assert_eq!(proj[1], 0.0);
        assert_eq!(proj[2], 0.0);
    }

    #[test]
    fn test_project_parallel() {
        let v = vec3(2.0, 4.0, 6.0);
        let onto = vec3(1.0, 2.0, 3.0);
        let proj = project(&v, &onto).unwrap();

        // Projecting onto parallel vector gives the original vector
        assert!((proj[0] - v[0]).abs() < 1e-10);
        assert!((proj[1] - v[1]).abs() < 1e-10);
        assert!((proj[2] - v[2]).abs() < 1e-10);
    }

    #[test]
    fn test_project_perpendicular() {
        let v = vec3(1.0, 0.0, 0.0);
        let onto = vec3(0.0, 1.0, 0.0);
        let proj = project(&v, &onto).unwrap();

        // Projecting onto perpendicular vector gives zero vector
        assert_eq!(proj[0], 0.0);
        assert_eq!(proj[1], 0.0);
        assert_eq!(proj[2], 0.0);
    }

    #[test]
    fn test_project_zero_onto() {
        let v = vec3(1.0, 2.0, 3.0);
        let zero = vec3(0.0, 0.0, 0.0);

        // Cannot project onto zero vector
        assert!(project(&v, &zero).is_none());
    }

    #[test]
    fn test_reject() {
        let v = vec3(3.0, 4.0, 0.0);
        let from = vec3(1.0, 0.0, 0.0);
        let rej = reject(&v, &from).unwrap();

        assert_eq!(rej[0], 0.0);
        assert_eq!(rej[1], 4.0);
        assert_eq!(rej[2], 0.0);

        // Rejection is perpendicular to 'from'
        assert_eq!(dot_product(&rej, &from), 0.0);
    }

    #[test]
    fn test_project_reject_decomposition() {
        let v = vec3(3.0, 4.0, 5.0);
        let onto = vec3(1.0, 1.0, 0.0);

        let proj = project(&v, &onto).unwrap();
        let rej = reject(&v, &onto).unwrap();

        // v = proj + rej
        let reconstructed = &proj + &rej;
        assert!((reconstructed[0] - v[0]).abs() < 1e-10);
        assert!((reconstructed[1] - v[1]).abs() < 1e-10);
        assert!((reconstructed[2] - v[2]).abs() < 1e-10);

        // proj and rej are perpendicular
        let dot = dot_product(&proj, &rej);
        assert!(dot.abs() < 1e-10);
    }

    #[test]
    fn test_reject_perpendicular() {
        let v = vec3(1.0, 0.0, 0.0);
        let from = vec3(0.0, 1.0, 0.0);
        let rej = reject(&v, &from).unwrap();

        // Rejecting from perpendicular vector gives original vector
        assert!((rej[0] - v[0]).abs() < 1e-10);
        assert!((rej[1] - v[1]).abs() < 1e-10);
        assert!((rej[2] - v[2]).abs() < 1e-10);
    }

    #[test]
    fn test_reflect() {
        let v = vec3(1.0, -1.0, 0.0);
        let normal = vec3(0.0, 1.0, 0.0); // Reflect across xz-plane
        let reflected = reflect(&v, &normal).unwrap();

        assert!((reflected[0] - 1.0).abs() < 1e-10);
        assert!((reflected[1] - 1.0).abs() < 1e-10); // y-component flipped
        assert!((reflected[2] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_reflect_magnitude_preserved() {
        let v = vec3(3.0, 4.0, 5.0);
        let normal = vec3(1.0, 0.0, 0.0);
        let reflected = reflect(&v, &normal).unwrap();

        // Magnitude should be preserved
        let v_mag = magnitude(&v);
        let r_mag = magnitude(&reflected);
        assert!((v_mag - r_mag).abs() < 1e-10);
    }

    #[test]
    fn test_reflect_twice() {
        let v = vec3(3.0, 4.0, 5.0);
        let normal = vec3(1.0, 1.0, 0.0);

        let reflected_once = reflect(&v, &normal).unwrap();
        let reflected_twice = reflect(&reflected_once, &normal).unwrap();

        // Reflecting twice should give original vector
        assert!((reflected_twice[0] - v[0]).abs() < 1e-10);
        assert!((reflected_twice[1] - v[1]).abs() < 1e-10);
        assert!((reflected_twice[2] - v[2]).abs() < 1e-10);
    }

    #[test]
    fn test_reflect_perpendicular_to_normal() {
        let normal = vec3(0.0, 1.0, 0.0);
        let v = vec3(1.0, 0.0, 0.0); // Perpendicular to normal

        let reflected = reflect(&v, &normal).unwrap();

        // Vector perpendicular to normal should be unchanged
        assert!((reflected[0] - v[0]).abs() < 1e-10);
        assert!((reflected[1] - v[1]).abs() < 1e-10);
        assert!((reflected[2] - v[2]).abs() < 1e-10);
    }

    #[test]
    fn test_reflect_parallel_to_normal() {
        let normal = vec3(0.0, 1.0, 0.0);
        let v = vec3(0.0, 2.0, 0.0); // Parallel to normal

        let reflected = reflect(&v, &normal).unwrap();

        // Vector parallel to normal should be reversed
        assert!((reflected[0] - 0.0).abs() < 1e-10);
        assert!((reflected[1] - (-2.0)).abs() < 1e-10);
        assert!((reflected[2] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_reflect_zero_normal() {
        let v = vec3(1.0, 2.0, 3.0);
        let zero = vec3(0.0, 0.0, 0.0);

        // Cannot reflect across zero normal
        assert!(reflect(&v, &zero).is_none());
    }

    #[test]
    fn test_reflect_with_non_unit_normal() {
        let v = vec3(1.0, -1.0, 0.0);
        let normal = vec3(0.0, 2.0, 0.0); // Non-unit normal
        let reflected = reflect(&v, &normal).unwrap();

        // Should still work correctly with non-unit normal
        assert!((reflected[0] - 1.0).abs() < 1e-10);
        assert!((reflected[1] - 1.0).abs() < 1e-10);
        assert!((reflected[2] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_project_magnitude() {
        let v = vec3(3.0, 4.0, 0.0);
        let onto = vec3(1.0, 0.0, 0.0);
        let proj = project(&v, &onto).unwrap();

        // Projection magnitude should not exceed original magnitude
        assert!(magnitude(&proj) <= magnitude(&v) + 1e-10);
    }

    #[test]
    fn test_pythagorean_theorem() {
        let v = vec3(3.0, 4.0, 0.0);
        let onto = vec3(1.0, 1.0, 0.0);

        let proj = project(&v, &onto).unwrap();
        let rej = reject(&v, &onto).unwrap();

        // ‖v‖² = ‖proj‖² + ‖rej‖²
        let v_mag_sq = magnitude(&v) * magnitude(&v);
        let proj_mag_sq = magnitude(&proj) * magnitude(&proj);
        let rej_mag_sq = magnitude(&rej) * magnitude(&rej);

        assert!((v_mag_sq - (proj_mag_sq + rej_mag_sq)).abs() < 1e-10);
    }
}
