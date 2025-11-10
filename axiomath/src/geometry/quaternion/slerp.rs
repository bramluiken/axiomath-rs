//! Spherical linear interpolation (slerp) for quaternions.

use super::operations::{magnitude, normalize, Quaternion};
use num_traits::{Float, One, Zero};

/// Performs spherical linear interpolation between two quaternions.
///
/// # Mathematical Background
///
/// Slerp (Spherical Linear Interpolation) provides smooth interpolation between
/// two rotations represented as unit quaternions. Unlike simple linear interpolation,
/// slerp maintains constant angular velocity and produces the shortest path on
/// the 4D unit hypersphere.
///
/// Given two unit quaternions q₀ and q₁, and parameter t ∈ [0, 1], slerp is defined as:
///
/// ```text
/// slerp(q₀, q₁, t) = (sin((1-t)θ)/sin(θ))q₀ + (sin(tθ)/sin(θ))q₁
/// ```
///
/// where θ is the angle between q₀ and q₁, computed from their dot product:
///
/// ```text
/// cos(θ) = q₀ · q₁ = w₀w₁ + x₀x₁ + y₀y₁ + z₀z₁
/// ```
///
/// Special cases:
/// - When θ ≈ 0 (quaternions very close), fall back to linear interpolation
/// - When cos(θ) < 0 (obtuse angle), negate one quaternion to take shorter path
/// - When t = 0, returns q₀
/// - When t = 1, returns q₁
///
/// # Algorithm
///
/// 1. Normalize both input quaternions
/// 2. Compute dot product (cosine of angle between them)
/// 3. If dot product is negative, negate q₁ to take shorter path
/// 4. If quaternions are very close (|cos θ| ≈ 1), use linear interpolation
/// 5. Otherwise, compute slerp using the formula above
/// 6. Normalize the result to ensure unit quaternion
///
/// # Complexity
///
/// - Time: O(1) - constant number of trigonometric operations
/// - Space: O(1)
///
/// # Parameters
///
/// * `q0` - The starting quaternion (will be normalized)
/// * `q1` - The ending quaternion (will be normalized)
/// * `t` - The interpolation parameter, typically in [0, 1]
///   - t = 0 returns q0
///   - t = 1 returns q1
///   - t ∈ (0, 1) returns an intermediate rotation
///
/// # Returns
///
/// A unit quaternion representing the interpolated rotation
///
/// # References
///
/// - Shoemake, K. (1985). "Animating rotation with quaternion curves", SIGGRAPH
/// - "Visualizing Quaternions" by Andrew J. Hanson
/// - "Understanding Slerp, Then Not Using It" by Jonathan Blow
///
/// # Example
///
/// ```
/// use axiomath::geometry::quaternion::{from_axis_angle, slerp};
/// use axiomath::vector::Vec3;
/// use std::f64::consts::PI;
///
/// // Interpolate between 0° and 90° rotation around Z
/// let axis = Vec3::new(0.0, 0.0, 1.0);
/// let q0 = from_axis_angle(axis, 0.0);
/// let q1 = from_axis_angle(axis, PI / 2.0);
///
/// // At t=0.5, should be 45° rotation
/// let mid = slerp(q0, q1, 0.5);
///
/// // Verify it's a unit quaternion
/// use axiomath::geometry::quaternion::magnitude;
/// assert!((magnitude(mid) - 1.0).abs() < 1e-10);
/// ```
pub fn slerp<T>(q0: Quaternion<T>, q1: Quaternion<T>, t: T) -> Quaternion<T>
where
    T: Float + Zero + One,
{
    // Normalize input quaternions
    let q0 = normalize(q0);
    let mut q1 = normalize(q1);

    // Compute dot product (cosine of angle between quaternions)
    let mut dot = q0.w * q1.w + q0.x * q1.x + q0.y * q1.y + q0.z * q1.z;

    // If dot product is negative, negate one quaternion to take shorter path
    // q and -q represent the same rotation, so we choose the one with positive dot
    if dot < T::zero() {
        q1 = Quaternion::new(-q1.w, -q1.x, -q1.y, -q1.z);
        dot = -dot;
    }

    // Clamp dot to avoid numerical issues with acos
    let dot = if dot > T::one() {
        T::one()
    } else if dot < -T::one() {
        -T::one()
    } else {
        dot
    };

    // Threshold for falling back to linear interpolation
    // When quaternions are very close, slerp becomes numerically unstable
    let threshold = T::one() - T::epsilon() * (T::one() + T::one() + T::one() + T::one()); // 1 - 4ε

    if dot > threshold {
        // Quaternions are very close, use linear interpolation (lerp)
        // This avoids division by small numbers in sin(θ)
        let result = Quaternion::new(
            q0.w + t * (q1.w - q0.w),
            q0.x + t * (q1.x - q0.x),
            q0.y + t * (q1.y - q0.y),
            q0.z + t * (q1.z - q0.z),
        );
        return normalize(result);
    }

    // Compute the angle between the quaternions
    let theta = dot.acos();
    let sin_theta = theta.sin();

    // Compute the slerp coefficients
    let coeff0 = ((T::one() - t) * theta).sin() / sin_theta;
    let coeff1 = (t * theta).sin() / sin_theta;

    // Compute the interpolated quaternion
    let result = Quaternion::new(
        coeff0 * q0.w + coeff1 * q1.w,
        coeff0 * q0.x + coeff1 * q1.x,
        coeff0 * q0.y + coeff1 * q1.y,
        coeff0 * q0.z + coeff1 * q1.z,
    );

    // Normalize to ensure unit quaternion (guards against numerical errors)
    normalize(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::quaternion::from_axis_angle;
    use crate::vector::vec3;
    use num_traits::Float;
    use std::f64::consts::PI;

    #[test]
    fn test_slerp_at_t_zero() {
        let axis = vec3(0.0, 0.0, 1.0);
        let q0 = from_axis_angle(axis, 0.0);
        let q1 = from_axis_angle(axis, PI / 2.0);

        let result = slerp(q0, q1, 0.0);

        // Should return q0
        assert!(Float::abs(result.w - q0.w) < 1e-10);
        assert!(Float::abs(result.x - q0.x) < 1e-10);
        assert!(Float::abs(result.y - q0.y) < 1e-10);
        assert!(Float::abs(result.z - q0.z) < 1e-10);
    }

    #[test]
    fn test_slerp_at_t_one() {
        let axis = vec3(0.0, 0.0, 1.0);
        let q0 = from_axis_angle(axis, 0.0);
        let q1 = from_axis_angle(axis, PI / 2.0);

        let result = slerp(q0, q1, 1.0);

        // Should return q1
        assert!(Float::abs(result.w - q1.w) < 1e-10);
        assert!(Float::abs(result.x - q1.x) < 1e-10);
        assert!(Float::abs(result.y - q1.y) < 1e-10);
        assert!(Float::abs(result.z - q1.z) < 1e-10);
    }

    #[test]
    fn test_slerp_midpoint() {
        let axis = vec3(0.0, 0.0, 1.0);
        let q0 = from_axis_angle(axis, 0.0);
        let q1 = from_axis_angle(axis, PI / 2.0);

        let mid = slerp(q0, q1, 0.5);

        // At t=0.5, should be approximately 45° rotation
        let expected = from_axis_angle(axis, PI / 4.0);

        assert!(Float::abs(mid.w - expected.w) < 1e-10);
        assert!(Float::abs(mid.x - expected.x) < 1e-10);
        assert!(Float::abs(mid.y - expected.y) < 1e-10);
        assert!(Float::abs(mid.z - expected.z) < 1e-10);
    }

    #[test]
    fn test_slerp_produces_unit_quaternion() {
        let axis = vec3(1.0, 1.0, 1.0);
        let q0 = from_axis_angle(axis, 0.0);
        let q1 = from_axis_angle(axis, PI);

        for i in 0..=10 {
            let t = i as f64 / 10.0;
            let result = slerp(q0, q1, t);
            let mag = magnitude(result);
            assert!(
                Float::abs(mag - 1.0) < 1e-10,
                "Magnitude at t={} is {}, expected 1.0",
                t,
                mag
            );
        }
    }

    #[test]
    fn test_slerp_handles_opposite_quaternions() {
        // Test that slerp handles quaternions pointing in opposite directions
        let q0 = Quaternion::new(1.0, 0.0, 0.0, 0.0);
        let q1 = Quaternion::new(-1.0, 0.0, 0.0, 0.0); // Opposite

        let result = slerp(q0, q1, 0.5);

        // Should still produce a valid unit quaternion
        let mag = magnitude(result);
        assert!(Float::abs(mag - 1.0) < 1e-10);
    }

    #[test]
    fn test_slerp_different_axes() {
        // Slerp between rotation around X and rotation around Y
        let axis_x = vec3(1.0, 0.0, 0.0);
        let axis_y = vec3(0.0, 1.0, 0.0);

        let q0 = from_axis_angle(axis_x, PI / 4.0);
        let q1 = from_axis_angle(axis_y, PI / 4.0);

        let result = slerp(q0, q1, 0.5);

        // Should be a unit quaternion
        let mag = magnitude(result);
        assert!(Float::abs(mag - 1.0) < 1e-10);
    }

    #[test]
    fn test_slerp_multiple_steps() {
        let axis = vec3(0.0, 1.0, 0.0);
        let q0 = from_axis_angle(axis, 0.0);
        let q1 = from_axis_angle(axis, PI);

        // Verify smooth interpolation at multiple points
        let mut prev_angle = 0.0;
        for i in 1..=10 {
            let t = i as f64 / 10.0;
            let result = slerp(q0, q1, t);

            // Compute angle from quaternion
            let angle = 2.0 * result.w.acos();

            // Angle should increase monotonically
            assert!(
                angle >= prev_angle,
                "Angle should increase: prev={}, current={}",
                prev_angle,
                angle
            );
            prev_angle = angle;
        }
    }

    #[test]
    fn test_slerp_close_quaternions() {
        // Test the linear interpolation fallback for very close quaternions
        let axis = vec3(1.0, 0.0, 0.0);
        let q0 = from_axis_angle(axis, 0.0);
        let q1 = from_axis_angle(axis, 1e-10); // Very close

        let result = slerp(q0, q1, 0.5);

        // Should still be a valid unit quaternion
        let mag = magnitude(result);
        assert!(Float::abs(mag - 1.0) < 1e-10);
    }
}
