//! Quaternion type and operations.
//!
//! This module provides a comprehensive quaternion implementation for
//! representing and manipulating 3D rotations.

use crate::core::traits::Scalar;
use crate::matrix::{Mat3, Mat4};
use crate::vector::{vec3, Vec3};
use num_traits::{Float, One, Zero};

/// A quaternion representing a rotation in 3D space.
///
/// # Mathematical Background
///
/// A quaternion is a four-dimensional number system that extends complex numbers.
/// Quaternions are particularly useful for representing rotations in 3D space
/// without suffering from gimbal lock, which affects Euler angles.
///
/// A quaternion q is represented as:
///
/// ```text
/// q = w + xi + yj + zk
/// ```
///
/// where w is the scalar (real) part and (x, y, z) is the vector (imaginary) part.
/// The units i, j, k satisfy:
///
/// ```text
/// i² = j² = k² = ijk = -1
/// ```
///
/// For unit quaternions (|q| = 1), there is a correspondence with rotations:
/// - A rotation by angle θ around axis **n** = [nx, ny, nz] is represented by:
///   ```text
///   q = cos(θ/2) + sin(θ/2)(nx*i + ny*j + nz*k)
///   ```
///
/// # Representation
///
/// This implementation uses the (w, x, y, z) convention where:
/// - `w`: scalar/real component
/// - `x`: i component (imaginary)
/// - `y`: j component (imaginary)
/// - `z`: k component (imaginary)
///
/// # References
///
/// - Hamilton, W. R. (1844). "On Quaternions"
/// - Shoemake, K. (1985). "Animating rotation with quaternion curves"
/// - "3D Math Primer for Graphics and Game Development" by Dunn & Parberry
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quaternion<T> {
    /// Scalar (real) component
    pub w: T,
    /// i component (imaginary)
    pub x: T,
    /// j component (imaginary)
    pub y: T,
    /// k component (imaginary)
    pub z: T,
}

impl<T> Quaternion<T>
where
    T: Float + Zero + One,
{
    /// Creates a new quaternion from components.
    ///
    /// # Parameters
    ///
    /// * `w` - Scalar (real) component
    /// * `x` - i component
    /// * `y` - j component
    /// * `z` - k component
    ///
    /// # Example
    ///
    /// ```
    /// use axiomath::geometry::quaternion::Quaternion;
    ///
    /// let q = Quaternion::new(1.0, 0.0, 0.0, 0.0);
    /// assert_eq!(q.w, 1.0);
    /// ```
    pub fn new(w: T, x: T, y: T, z: T) -> Self {
        Quaternion { w, x, y, z }
    }

    /// Creates an identity quaternion representing no rotation.
    ///
    /// The identity quaternion is (1, 0, 0, 0).
    ///
    /// # Example
    ///
    /// ```
    /// use axiomath::geometry::quaternion::Quaternion;
    ///
    /// let q = Quaternion::<f64>::identity();
    /// assert_eq!(q.w, 1.0);
    /// assert_eq!(q.x, 0.0);
    /// ```
    pub fn identity() -> Self {
        Quaternion::new(T::one(), T::zero(), T::zero(), T::zero())
    }
}

/// Multiplies two quaternions using the Hamilton product.
///
/// # Mathematical Background
///
/// The Hamilton product (quaternion multiplication) is defined as:
///
/// ```text
/// q₁ * q₂ = (w₁w₂ - x₁x₂ - y₁y₂ - z₁z₂)
///         + (w₁x₂ + x₁w₂ + y₁z₂ - z₁y₂)i
///         + (w₁y₂ - x₁z₂ + y₁w₂ + z₁x₂)j
///         + (w₁z₂ + x₁y₂ - y₁x₂ + z₁w₂)k
/// ```
///
/// Quaternion multiplication is **not commutative**: q₁ * q₂ ≠ q₂ * q₁ in general.
/// However, it is associative: (q₁ * q₂) * q₃ = q₁ * (q₂ * q₃).
///
/// When representing rotations, multiplying quaternions corresponds to
/// composing rotations: q₁ * q₂ means "first apply rotation q₂, then q₁".
///
/// # Algorithm
///
/// 1. Compute the scalar part using the formula above
/// 2. Compute each of the three imaginary parts
/// 3. Construct the resulting quaternion
///
/// # Complexity
///
/// - Time: O(1) - requires 16 multiplications and 12 additions
/// - Space: O(1)
///
/// # Parameters
///
/// * `q1` - The first quaternion
/// * `q2` - The second quaternion
///
/// # Returns
///
/// The Hamilton product q1 * q2
///
/// # Example
///
/// ```
/// use axiomath::geometry::quaternion::{Quaternion, multiply};
///
/// let q1 = Quaternion::new(1.0, 0.0, 0.0, 0.0);
/// let q2 = Quaternion::new(0.0, 1.0, 0.0, 0.0);
/// let result = multiply(q1, q2);
/// assert_eq!(result.x, 1.0);
/// ```
pub fn multiply<T>(q1: Quaternion<T>, q2: Quaternion<T>) -> Quaternion<T>
where
    T: Float + Zero + One,
{
    let w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;
    let x = q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y;
    let y = q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x;
    let z = q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w;

    Quaternion::new(w, x, y, z)
}

/// Computes the conjugate of a quaternion.
///
/// # Mathematical Background
///
/// The conjugate of a quaternion q = w + xi + yj + zk is defined as:
///
/// ```text
/// q* = w - xi - yj - zk
/// ```
///
/// The conjugate has important properties:
/// - For unit quaternions: q* = q⁻¹ (the conjugate equals the inverse)
/// - q * q* = |q|² (a real number)
/// - (q₁ * q₂)* = q₂* * q₁* (reverses order of multiplication)
///
/// For rotation quaternions, the conjugate represents the inverse rotation.
///
/// # Algorithm
///
/// 1. Keep the scalar component w unchanged
/// 2. Negate all imaginary components (x, y, z)
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Parameters
///
/// * `q` - The quaternion to conjugate
///
/// # Returns
///
/// The conjugate quaternion
///
/// # Example
///
/// ```
/// use axiomath::geometry::quaternion::{Quaternion, conjugate};
///
/// let q = Quaternion::new(1.0, 2.0, 3.0, 4.0);
/// let conj = conjugate(q);
/// assert_eq!(conj.w, 1.0);
/// assert_eq!(conj.x, -2.0);
/// assert_eq!(conj.y, -3.0);
/// assert_eq!(conj.z, -4.0);
/// ```
pub fn conjugate<T>(q: Quaternion<T>) -> Quaternion<T>
where
    T: Float,
{
    Quaternion::new(q.w, -q.x, -q.y, -q.z)
}

/// Computes the magnitude (norm) of a quaternion.
///
/// # Mathematical Background
///
/// The magnitude (or norm) of a quaternion q = w + xi + yj + zk is:
///
/// ```text
/// |q| = sqrt(w² + x² + y² + z²)
/// ```
///
/// This is analogous to the magnitude of a 4D vector. Properties:
/// - |q₁ * q₂| = |q₁| * |q₂|
/// - A quaternion is called "unit" if |q| = 1
/// - Unit quaternions represent pure rotations (no scaling)
///
/// # Algorithm
///
/// 1. Sum the squares of all four components
/// 2. Take the square root
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Parameters
///
/// * `q` - The quaternion
///
/// # Returns
///
/// The magnitude of the quaternion
///
/// # Example
///
/// ```
/// use axiomath::geometry::quaternion::{Quaternion, magnitude};
///
/// let q = Quaternion::new(1.0, 0.0, 0.0, 0.0);
/// let mag = magnitude(q);
/// assert!((mag - 1.0).abs() < 1e-10);
/// ```
pub fn magnitude<T>(q: Quaternion<T>) -> T
where
    T: Float,
{
    (q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z).sqrt()
}

/// Normalizes a quaternion to unit length.
///
/// # Mathematical Background
///
/// A normalized (unit) quaternion has magnitude 1. Normalization is computed as:
///
/// ```text
/// q_normalized = q / |q|
/// ```
///
/// Unit quaternions are essential for representing rotations because:
/// - They preserve distances (no scaling)
/// - Interpolation between unit quaternions (slerp) gives smooth rotations
/// - Concatenating rotations preserves the unit property
///
/// # Algorithm
///
/// 1. Compute the magnitude |q|
/// 2. Divide all components by |q|
/// 3. Handle the edge case where |q| = 0
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Parameters
///
/// * `q` - The quaternion to normalize
///
/// # Returns
///
/// The normalized quaternion, or the identity quaternion if the input has zero magnitude
///
/// # Example
///
/// ```
/// use axiomath::geometry::quaternion::{Quaternion, normalize, magnitude};
///
/// let q = Quaternion::new(2.0, 0.0, 0.0, 0.0);
/// let normalized = normalize(q);
/// let mag = magnitude(normalized);
/// assert!((mag - 1.0).abs() < 1e-10);
/// ```
pub fn normalize<T>(q: Quaternion<T>) -> Quaternion<T>
where
    T: Float + Zero + One,
{
    let mag = magnitude(q);

    // Handle zero-magnitude quaternion
    if mag <= T::epsilon() {
        return Quaternion::identity();
    }

    Quaternion::new(q.w / mag, q.x / mag, q.y / mag, q.z / mag)
}

/// Computes the inverse of a quaternion.
///
/// # Mathematical Background
///
/// The inverse of a quaternion q is defined as the quaternion q⁻¹ such that:
///
/// ```text
/// q * q⁻¹ = q⁻¹ * q = 1 (identity)
/// ```
///
/// The inverse is computed as:
///
/// ```text
/// q⁻¹ = q* / |q|²
/// ```
///
/// where q* is the conjugate and |q| is the magnitude.
///
/// For unit quaternions (|q| = 1), the inverse equals the conjugate: q⁻¹ = q*.
///
/// # Algorithm
///
/// 1. Compute the conjugate q*
/// 2. Compute |q|²
/// 3. Divide each component of q* by |q|²
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Parameters
///
/// * `q` - The quaternion to invert
///
/// # Returns
///
/// The inverse quaternion, or the identity quaternion if the input has zero magnitude
///
/// # Example
///
/// ```
/// use axiomath::geometry::quaternion::{Quaternion, inverse, multiply};
///
/// let q = Quaternion::new(1.0, 2.0, 3.0, 4.0);
/// let q_inv = inverse(q);
/// let result = multiply(q, q_inv);
///
/// // q * q⁻¹ should be approximately identity
/// assert!((result.w - 1.0).abs() < 1e-10);
/// assert!(result.x.abs() < 1e-10);
/// ```
pub fn inverse<T>(q: Quaternion<T>) -> Quaternion<T>
where
    T: Float + Zero + One,
{
    let mag_squared = q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z;

    // Handle zero-magnitude quaternion
    if mag_squared <= T::epsilon() {
        return Quaternion::identity();
    }

    let conj = conjugate(q);
    Quaternion::new(
        conj.w / mag_squared,
        conj.x / mag_squared,
        conj.y / mag_squared,
        conj.z / mag_squared,
    )
}

/// Creates a quaternion from an axis-angle representation.
///
/// # Mathematical Background
///
/// A rotation by angle θ around a unit axis **n** = [nx, ny, nz] can be
/// represented as a quaternion:
///
/// ```text
/// q = cos(θ/2) + sin(θ/2)(nx*i + ny*j + nz*k)
///   = [cos(θ/2), sin(θ/2)*nx, sin(θ/2)*ny, sin(θ/2)*nz]
/// ```
///
/// Note that the angle is halved. This is because quaternions represent
/// rotations in a double-covering of SO(3): q and -q represent the same rotation.
///
/// # Algorithm
///
/// 1. Normalize the input axis to ensure it's a unit vector
/// 2. Compute half_angle = θ/2
/// 3. Compute cos(half_angle) and sin(half_angle)
/// 4. Construct quaternion components
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Parameters
///
/// * `axis` - The rotation axis (will be normalized)
/// * `angle` - The rotation angle in radians
///
/// # Returns
///
/// A unit quaternion representing the rotation
///
/// # Panics
///
/// Panics if the axis is the zero vector.
///
/// # Example
///
/// ```
/// use axiomath::geometry::quaternion::from_axis_angle;
/// use axiomath::vector::vec3;
/// use std::f64::consts::PI;
///
/// // 180 degree rotation around Z axis
/// let axis = vec3(0.0, 0.0, 1.0);
/// let q = from_axis_angle(axis, PI);
/// assert!(q.w.abs() < 1e-10); // cos(90°) = 0
/// assert!((q.z - 1.0).abs() < 1e-10); // sin(90°) = 1
/// ```
pub fn from_axis_angle<T>(axis: Vec3<T>, angle: T) -> Quaternion<T>
where
    T: Scalar + Float + Zero + One,
{
    // Normalize the axis
    let norm = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    assert!(
        norm > T::zero(),
        "Axis must be non-zero for from_axis_angle"
    );

    let nx = axis[0] / norm;
    let ny = axis[1] / norm;
    let nz = axis[2] / norm;

    // Compute half angle
    let half_angle = angle / (T::one() + T::one()); // angle / 2
    let cos_half = half_angle.cos();
    let sin_half = half_angle.sin();

    Quaternion::new(cos_half, sin_half * nx, sin_half * ny, sin_half * nz)
}

/// Converts a quaternion to a 3x3 rotation matrix.
///
/// # Mathematical Background
///
/// A unit quaternion q = w + xi + yj + zk can be converted to a 3x3 rotation matrix:
///
/// ```text
/// R = [ 1-2(y²+z²)   2(xy-wz)     2(xz+wy)   ]
///     [ 2(xy+wz)     1-2(x²+z²)   2(yz-wx)   ]
///     [ 2(xz-wy)     2(yz+wx)     1-2(x²+y²) ]
/// ```
///
/// This derivation comes from the fact that rotating a vector **v** by quaternion q
/// is equivalent to computing: q * [0, v] * q*, where the middle term treats **v**
/// as a pure quaternion (w=0).
///
/// # Algorithm
///
/// 1. Normalize the quaternion to ensure it's unit length
/// 2. Compute commonly used terms (2x², 2y², etc.)
/// 3. Construct the matrix elements using the formula above
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Parameters
///
/// * `q` - The quaternion (will be normalized)
///
/// # Returns
///
/// A 3x3 rotation matrix
///
/// # Example
///
/// ```
/// use axiomath::geometry::quaternion::{from_axis_angle, to_rotation_matrix_3x3};
/// use axiomath::vector::vec3;
/// use std::f64::consts::PI;
///
/// let axis = vec3(0.0, 0.0, 1.0);
/// let q = from_axis_angle(axis, PI / 2.0);
/// let matrix = to_rotation_matrix_3x3(q);
///
/// // Should be a 90-degree rotation around Z
/// assert!(matrix[(0, 0)].abs() < 1e-10); // cos(90°) ≈ 0
/// ```
pub fn to_rotation_matrix_3x3<T>(q: Quaternion<T>) -> Mat3<T>
where
    T: Scalar + Float + Zero + One,
{
    // Normalize to ensure unit quaternion
    let q = normalize(q);

    let two = T::one() + T::one();

    let xx = two * q.x * q.x;
    let yy = two * q.y * q.y;
    let zz = two * q.z * q.z;
    let xy = two * q.x * q.y;
    let xz = two * q.x * q.z;
    let yz = two * q.y * q.z;
    let wx = two * q.w * q.x;
    let wy = two * q.w * q.y;
    let wz = two * q.w * q.z;

    Mat3::from_rows([
        [T::one() - yy - zz, xy - wz, xz + wy],
        [xy + wz, T::one() - xx - zz, yz - wx],
        [xz - wy, yz + wx, T::one() - xx - yy],
    ])
}

/// Converts a quaternion to a 4x4 rotation matrix (homogeneous coordinates).
///
/// # Mathematical Background
///
/// This is the same as [`to_rotation_matrix_3x3`], but returns a 4x4 matrix
/// suitable for use in homogeneous coordinate transformations. The 4x4 matrix
/// has the form:
///
/// ```text
/// [ R₃ₓ₃  0 ]
/// [  0    1 ]
/// ```
///
/// where R₃ₓ₃ is the 3x3 rotation matrix.
///
/// # Parameters
///
/// * `q` - The quaternion (will be normalized)
///
/// # Returns
///
/// A 4x4 rotation matrix in homogeneous coordinates
///
/// # Example
///
/// ```
/// use axiomath::geometry::quaternion::{from_axis_angle, to_rotation_matrix_4x4};
/// use axiomath::vector::vec3;
/// use std::f64::consts::PI;
///
/// let axis = vec3(1.0, 0.0, 0.0);
/// let q = from_axis_angle(axis, PI / 2.0);
/// let matrix = to_rotation_matrix_4x4(q);
///
/// assert_eq!(matrix[(3, 3)], 1.0);
/// ```
pub fn to_rotation_matrix_4x4<T>(q: Quaternion<T>) -> Mat4<T>
where
    T: Scalar + Float + Zero + One,
{
    // Normalize to ensure unit quaternion
    let q = normalize(q);

    let two = T::one() + T::one();

    let xx = two * q.x * q.x;
    let yy = two * q.y * q.y;
    let zz = two * q.z * q.z;
    let xy = two * q.x * q.y;
    let xz = two * q.x * q.z;
    let yz = two * q.y * q.z;
    let wx = two * q.w * q.x;
    let wy = two * q.w * q.y;
    let wz = two * q.w * q.z;

    Mat4::from_rows([
        [T::one() - yy - zz, xy - wz, xz + wy, T::zero()],
        [xy + wz, T::one() - xx - zz, yz - wx, T::zero()],
        [xz - wy, yz + wx, T::one() - xx - yy, T::zero()],
        [T::zero(), T::zero(), T::zero(), T::one()],
    ])
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::Float;
    use std::f64::consts::PI;

    #[test]
    fn test_quaternion_identity() {
        let q = Quaternion::<f64>::identity();
        assert_eq!(q.w, 1.0);
        assert_eq!(q.x, 0.0);
        assert_eq!(q.y, 0.0);
        assert_eq!(q.z, 0.0);
    }

    #[test]
    fn test_multiply_identity() {
        let q = Quaternion::new(1.0, 2.0, 3.0, 4.0);
        let identity = Quaternion::identity();
        let result = multiply(q, identity);
        assert_eq!(result, q);
    }

    #[test]
    fn test_multiply_associative() {
        let q1 = Quaternion::new(1.0, 0.5, 0.0, 0.0);
        let q2 = Quaternion::new(0.0, 1.0, 0.0, 0.0);
        let q3 = Quaternion::new(0.0, 0.0, 1.0, 0.0);

        let result1 = multiply(multiply(q1, q2), q3);
        let result2 = multiply(q1, multiply(q2, q3));

        assert!(Float::abs(result1.w - result2.w) < 1e-10);
        assert!(Float::abs(result1.x - result2.x) < 1e-10);
        assert!(Float::abs(result1.y - result2.y) < 1e-10);
        assert!(Float::abs(result1.z - result2.z) < 1e-10);
    }

    #[test]
    fn test_conjugate() {
        let q = Quaternion::new(1.0, 2.0, 3.0, 4.0);
        let conj = conjugate(q);
        assert_eq!(conj.w, 1.0);
        assert_eq!(conj.x, -2.0);
        assert_eq!(conj.y, -3.0);
        assert_eq!(conj.z, -4.0);
    }

    #[test]
    fn test_magnitude() {
        let q = Quaternion::new(1.0, 0.0, 0.0, 0.0);
        assert!(Float::abs(magnitude(q) - 1.0) < 1e-10);

        let q2 = Quaternion::new(2.0, 2.0, 2.0, 2.0);
        assert!(Float::abs(magnitude(q2) - 4.0) < 1e-10);
    }

    #[test]
    fn test_normalize() {
        let q = Quaternion::new(2.0, 0.0, 0.0, 0.0);
        let normalized = normalize(q);
        assert!(Float::abs(magnitude(normalized) - 1.0) < 1e-10);
        assert!(Float::abs(normalized.w - 1.0) < 1e-10);
    }

    #[test]
    fn test_inverse() {
        let q = Quaternion::new(1.0, 2.0, 3.0, 4.0);
        let q_inv = inverse(q);
        let result = multiply(q, q_inv);

        // Should be approximately identity
        assert!(Float::abs(result.w - 1.0) < 1e-10);
        assert!(Float::abs(result.x) < 1e-10);
        assert!(Float::abs(result.y) < 1e-10);
        assert!(Float::abs(result.z) < 1e-10);
    }

    #[test]
    fn test_from_axis_angle_z() {
        let axis = vec3(0.0, 0.0, 1.0);
        let q = from_axis_angle(axis, PI);
        assert!(Float::abs(q.w) < 1e-10); // cos(π/2) ≈ 0
        assert!(Float::abs(magnitude(q) - 1.0) < 1e-10);
    }

    #[test]
    fn test_from_axis_angle_identity() {
        let axis = vec3(1.0, 0.0, 0.0);
        let q = from_axis_angle(axis, 0.0);
        assert!(Float::abs(q.w - 1.0) < 1e-10);
        assert!(Float::abs(q.x) < 1e-10);
    }

    #[test]
    fn test_to_rotation_matrix_identity() {
        let q = Quaternion::<f64>::identity();
        let mat = to_rotation_matrix_3x3(q);
        let identity = Mat3::<f64>::identity();

        for i in 0..3 {
            for j in 0..3 {
                assert!(Float::abs(mat[(i, j)] - identity[(i, j)]) < 1e-10);
            }
        }
    }

    #[test]
    fn test_quaternion_rotation_consistency() {
        // A 90-degree rotation around Z should match expected matrix
        let axis = vec3(0.0, 0.0, 1.0);
        let q = from_axis_angle(axis, PI / 2.0);
        let mat = to_rotation_matrix_3x3(q);

        // cos(90°) ≈ 0, sin(90°) ≈ 1
        assert!(Float::abs(mat[(0, 0)]) < 1e-10);
        assert!(Float::abs(mat[(0, 1)] + 1.0) < 1e-10);
        assert!(Float::abs(mat[(1, 0)] - 1.0) < 1e-10);
        assert!(Float::abs(mat[(1, 1)]) < 1e-10);
    }

    #[test]
    fn test_to_rotation_matrix_4x4() {
        let q = Quaternion::<f64>::identity();
        let mat = to_rotation_matrix_4x4(q);
        assert_eq!(mat[(3, 3)], 1.0);
        assert_eq!(mat[(0, 3)], 0.0);
        assert_eq!(mat[(3, 0)], 0.0);
    }
}
