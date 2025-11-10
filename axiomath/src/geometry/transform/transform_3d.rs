//! 3D transformation matrices using homogeneous coordinates.
//!
//! This module provides functions to create 4x4 transformation matrices
//! for 3D geometric transformations including translation, rotation around
//! coordinate axes, rotation around arbitrary axes, and scaling.

use crate::core::traits::Scalar;
use crate::matrix::{transpose, Mat4};
use crate::vector::{vec3, Vec3};
use num_traits::{Float, One, Zero};
use std::ops::{Add, Mul, Neg};

/// Creates a 3D translation matrix.
///
/// # Mathematical Background
///
/// Translation in 3D homogeneous coordinates is represented by a 4x4 matrix:
///
/// ```text
/// [ 1  0  0  x ]
/// [ 0  1  0  y ]
/// [ 0  0  1  z ]
/// [ 0  0  0  1 ]
/// ```
///
/// When multiplied with a homogeneous point [px, py, pz, 1]^T, this produces
/// [px + x, py + y, pz + z, 1]^T, translating the point by (x, y, z).
///
/// # Algorithm
///
/// 1. Create an identity 4x4 matrix
/// 2. Set the translation components in the fourth column
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Parameters
///
/// * `x` - Translation amount along the x-axis
/// * `y` - Translation amount along the y-axis
/// * `z` - Translation amount along the z-axis
///
/// # Returns
///
/// A 4x4 translation matrix in homogeneous coordinates
///
/// # Example
///
/// ```
/// use axiomath::geometry::transform::translation_matrix_3d;
/// use axiomath::matrix::Mat4;
///
/// let t = translation_matrix_3d(5.0, 3.0, 2.0);
/// assert_eq!(t[(0, 3)], 5.0);
/// assert_eq!(t[(1, 3)], 3.0);
/// assert_eq!(t[(2, 3)], 2.0);
/// ```
pub fn translation_matrix_3d<T>(x: T, y: T, z: T) -> Mat4<T>
where
    T: Scalar + Float + Zero + One,
{
    Mat4::from_rows([
        [T::one(), T::zero(), T::zero(), x],
        [T::zero(), T::one(), T::zero(), y],
        [T::zero(), T::zero(), T::one(), z],
        [T::zero(), T::zero(), T::zero(), T::one()],
    ])
}

/// Creates a 3D rotation matrix around the X-axis.
///
/// # Mathematical Background
///
/// Rotation around the X-axis in 3D homogeneous coordinates:
///
/// ```text
/// [ 1    0        0      0 ]
/// [ 0  cos(θ)  -sin(θ)  0 ]
/// [ 0  sin(θ)   cos(θ)  0 ]
/// [ 0    0        0      1 ]
/// ```
///
/// This rotates points counterclockwise around the X-axis when looking
/// from positive X towards the origin (right-hand rule).
///
/// # Algorithm
///
/// 1. Compute cos(angle) and sin(angle)
/// 2. Construct the rotation matrix with X-axis unchanged
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Parameters
///
/// * `angle` - Rotation angle in radians (counterclockwise)
///
/// # Returns
///
/// A 4x4 rotation matrix around the X-axis
///
/// # Example
///
/// ```
/// use axiomath::geometry::transform::rotation_matrix_x;
/// use std::f64::consts::PI;
///
/// let r = rotation_matrix_x(PI / 2.0); // 90 degree rotation
/// assert!(r.get(1, 1).abs() < 1e-10); // cos(90°) ≈ 0
/// ```
pub fn rotation_matrix_x<T>(angle: T) -> Mat4<T>
where
    T: Scalar + Float + Zero + One,
{
    let cos_angle = angle.cos();
    let sin_angle = angle.sin();

    Mat4::from_rows([
        [T::one(), T::zero(), T::zero(), T::zero()],
        [T::zero(), cos_angle, -sin_angle, T::zero()],
        [T::zero(), sin_angle, cos_angle, T::zero()],
        [T::zero(), T::zero(), T::zero(), T::one()],
    ])
}

/// Creates a 3D rotation matrix around the Y-axis.
///
/// # Mathematical Background
///
/// Rotation around the Y-axis in 3D homogeneous coordinates:
///
/// ```text
/// [  cos(θ)  0  sin(θ)  0 ]
/// [    0     1    0     0 ]
/// [ -sin(θ)  0  cos(θ)  0 ]
/// [    0     0    0     1 ]
/// ```
///
/// This rotates points counterclockwise around the Y-axis when looking
/// from positive Y towards the origin (right-hand rule).
///
/// # Algorithm
///
/// 1. Compute cos(angle) and sin(angle)
/// 2. Construct the rotation matrix with Y-axis unchanged
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Parameters
///
/// * `angle` - Rotation angle in radians (counterclockwise)
///
/// # Returns
///
/// A 4x4 rotation matrix around the Y-axis
///
/// # Example
///
/// ```
/// use axiomath::geometry::transform::rotation_matrix_y;
/// use std::f64::consts::PI;
///
/// let r = rotation_matrix_y(PI / 2.0); // 90 degree rotation
/// assert!(r.get(0, 0).abs() < 1e-10); // cos(90°) ≈ 0
/// ```
pub fn rotation_matrix_y<T>(angle: T) -> Mat4<T>
where
    T: Scalar + Float + Zero + One,
{
    let cos_angle = angle.cos();
    let sin_angle = angle.sin();

    Mat4::from_rows([
        [cos_angle, T::zero(), sin_angle, T::zero()],
        [T::zero(), T::one(), T::zero(), T::zero()],
        [-sin_angle, T::zero(), cos_angle, T::zero()],
        [T::zero(), T::zero(), T::zero(), T::one()],
    ])
}

/// Creates a 3D rotation matrix around the Z-axis.
///
/// # Mathematical Background
///
/// Rotation around the Z-axis in 3D homogeneous coordinates:
///
/// ```text
/// [ cos(θ)  -sin(θ)  0  0 ]
/// [ sin(θ)   cos(θ)  0  0 ]
/// [   0        0     1  0 ]
/// [   0        0     0  1 ]
/// ```
///
/// This rotates points counterclockwise around the Z-axis when looking
/// from positive Z towards the origin (right-hand rule).
///
/// # Algorithm
///
/// 1. Compute cos(angle) and sin(angle)
/// 2. Construct the rotation matrix with Z-axis unchanged
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Parameters
///
/// * `angle` - Rotation angle in radians (counterclockwise)
///
/// # Returns
///
/// A 4x4 rotation matrix around the Z-axis
///
/// # Example
///
/// ```
/// use axiomath::geometry::transform::rotation_matrix_z;
/// use std::f64::consts::PI;
///
/// let r = rotation_matrix_z(PI / 2.0); // 90 degree rotation
/// assert!(r.get(0, 0).abs() < 1e-10); // cos(90°) ≈ 0
/// ```
pub fn rotation_matrix_z<T>(angle: T) -> Mat4<T>
where
    T: Scalar + Float + Zero + One,
{
    let cos_angle = angle.cos();
    let sin_angle = angle.sin();

    Mat4::from_rows([
        [cos_angle, -sin_angle, T::zero(), T::zero()],
        [sin_angle, cos_angle, T::zero(), T::zero()],
        [T::zero(), T::zero(), T::one(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::one()],
    ])
}

/// Creates a 3D rotation matrix around an arbitrary axis using Rodrigues' formula.
///
/// # Mathematical Background
///
/// Rodrigues' rotation formula provides a way to compute the rotation matrix
/// for a rotation by angle θ around an arbitrary unit axis vector **n** = [nx, ny, nz].
///
/// The rotation matrix R is given by:
///
/// ```text
/// R = I + sin(θ)K + (1 - cos(θ))K²
/// ```
///
/// where I is the identity matrix and K is the skew-symmetric cross-product matrix:
///
/// ```text
/// K = [  0   -nz   ny ]
///     [  nz   0   -nx ]
///     [ -ny   nx   0  ]
/// ```
///
/// Expanding this formula gives the explicit form:
///
/// ```text
/// R = [ nx²(1-c)+c    nxny(1-c)-nzs  nxnz(1-c)+nys ]
///     [ nxny(1-c)+nzs  ny²(1-c)+c    nynz(1-c)-nxs ]
///     [ nxnz(1-c)-nys  nynz(1-c)+nxs  nz²(1-c)+c   ]
/// ```
///
/// where c = cos(θ) and s = sin(θ).
///
/// # Algorithm
///
/// 1. Normalize the input axis to ensure it's a unit vector
/// 2. Compute cos(θ) and sin(θ)
/// 3. Apply Rodrigues' formula element by element
/// 4. Construct the 4x4 homogeneous matrix
///
/// # Complexity
///
/// - Time: O(1) - constant number of operations
/// - Space: O(1)
///
/// # Parameters
///
/// * `axis` - The axis of rotation (will be normalized)
/// * `angle` - Rotation angle in radians (counterclockwise)
///
/// # Returns
///
/// A 4x4 rotation matrix around the specified axis
///
/// # Panics
///
/// Panics if the axis is the zero vector (has zero length).
///
/// # References
///
/// - Rodrigues, O. (1840). "Des lois géométriques qui régissent les déplacements d'un système solide"
/// - "Fundamentals of Computer Graphics" by Shirley & Marschner, 4th edition, Chapter 6
///
/// # Example
///
/// ```
/// use axiomath::geometry::transform::rotation_matrix_axis_angle;
/// use axiomath::vector::Vec3;
/// use std::f64::consts::PI;
///
/// // Rotation around [1, 0, 0] should match rotation_matrix_x
/// let axis = vec3(1.0, 0.0, 0.0);
/// let r = rotation_matrix_axis_angle(axis, PI / 2.0);
/// assert!(r[(1, 1)].abs() < 1e-10); // cos(90°) ≈ 0
/// ```
pub fn rotation_matrix_axis_angle<T>(axis: Vec3<T>, angle: T) -> Mat4<T>
where
    T: Scalar + Float + Zero + One + Add<Output = T> + Mul<Output = T> + Neg<Output = T>,
{
    // Normalize the axis
    let norm = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    assert!(
        norm > T::zero(),
        "Axis must be non-zero for rotation_matrix_axis_angle"
    );

    let nx = axis[0] / norm;
    let ny = axis[1] / norm;
    let nz = axis[2] / norm;

    // Compute trigonometric values
    let c = angle.cos();
    let s = angle.sin();
    let one_minus_c = T::one() - c;

    // Apply Rodrigues' formula
    // R[i][j] = n[i]*n[j]*(1-cos) + δ[i][j]*cos + ε[i][j][k]*n[k]*sin
    // where δ is Kronecker delta and ε is Levi-Civita symbol

    let m00 = nx * nx * one_minus_c + c;
    let m01 = nx * ny * one_minus_c - nz * s;
    let m02 = nx * nz * one_minus_c + ny * s;

    let m10 = nx * ny * one_minus_c + nz * s;
    let m11 = ny * ny * one_minus_c + c;
    let m12 = ny * nz * one_minus_c - nx * s;

    let m20 = nx * nz * one_minus_c - ny * s;
    let m21 = ny * nz * one_minus_c + nx * s;
    let m22 = nz * nz * one_minus_c + c;

    Mat4::from_rows([
        [m00, m01, m02, T::zero()],
        [m10, m11, m12, T::zero()],
        [m20, m21, m22, T::zero()],
        [T::zero(), T::zero(), T::zero(), T::one()],
    ])
}

/// Creates a 3D scaling matrix with non-uniform scale factors.
///
/// # Mathematical Background
///
/// Scaling in 3D homogeneous coordinates is represented by:
///
/// ```text
/// [ sx  0   0   0 ]
/// [ 0   sy  0   0 ]
/// [ 0   0   sz  0 ]
/// [ 0   0   0   1 ]
/// ```
///
/// This scales points by sx, sy, sz along the x, y, z axes respectively.
/// When sx = sy = sz, this is uniform scaling. Otherwise, it's non-uniform.
///
/// # Algorithm
///
/// 1. Create a diagonal matrix with scale factors
/// 2. Set the homogeneous component to 1
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Parameters
///
/// * `sx` - Scale factor along the x-axis
/// * `sy` - Scale factor along the y-axis
/// * `sz` - Scale factor along the z-axis
///
/// # Returns
///
/// A 4x4 scaling matrix in homogeneous coordinates
///
/// # Example
///
/// ```
/// use axiomath::geometry::transform::scaling_matrix_3d;
/// use axiomath::matrix::Mat4;
///
/// let s = scaling_matrix_3d(2.0, 3.0, 4.0);
/// assert_eq!(s[(0, 0)], 2.0);
/// assert_eq!(s[(1, 1)], 3.0);
/// assert_eq!(s[(2, 2)], 4.0);
/// ```
pub fn scaling_matrix_3d<T>(sx: T, sy: T, sz: T) -> Mat4<T>
where
    T: Scalar + Float + Zero + One,
{
    Mat4::from_rows([
        [sx, T::zero(), T::zero(), T::zero()],
        [T::zero(), sy, T::zero(), T::zero()],
        [T::zero(), T::zero(), sz, T::zero()],
        [T::zero(), T::zero(), T::zero(), T::one()],
    ])
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_translation_matrix_3d() {
        let t = translation_matrix_3d(5.0, 3.0, 2.0);
        assert_eq!(t[(0, 3)], 5.0);
        assert_eq!(t[(1, 3)], 3.0);
        assert_eq!(t[(2, 3)], 2.0);
        assert_eq!(t[(3, 3)], 1.0);
    }

    #[test]
    fn test_translation_identity() {
        let t = translation_matrix_3d(0.0, 0.0, 0.0);
        let identity = Mat4::<f64>::identity();
        assert_eq!(t, identity);
    }

    #[test]
    fn test_rotation_matrix_x_zero() {
        let r = rotation_matrix_x(0.0);
        let identity = Mat4::<f64>::identity();
        for i in 0..4 {
            for j in 0..4 {
                assert!(Float::abs(r[(i, j)] - identity[(i, j)]) < 1e-10);
            }
        }
    }

    #[test]
    fn test_rotation_matrix_x_90_degrees() {
        let r = rotation_matrix_x(PI / 2.0);
        assert!(Float::abs(r[(0, 0)] - 1.0) < 1e-10);
        assert!(Float::abs(r[(1, 1)]) < 1e-10); // cos(90°) ≈ 0
        assert!(Float::abs(r[(1, 2)] + 1.0) < 1e-10); // -sin(90°) ≈ -1
        assert!(Float::abs(r[(2, 1)] - 1.0) < 1e-10); // sin(90°) ≈ 1
    }

    #[test]
    fn test_rotation_matrix_y_90_degrees() {
        let r = rotation_matrix_y(PI / 2.0);
        assert!(Float::abs(r[(0, 0)]) < 1e-10); // cos(90°) ≈ 0
        assert!(Float::abs(r[(1, 1)] - 1.0) < 1e-10);
        assert!(Float::abs(r[(0, 2)] - 1.0) < 1e-10); // sin(90°) ≈ 1
        assert!(Float::abs(r[(2, 0)] + 1.0) < 1e-10); // -sin(90°) ≈ -1
    }

    #[test]
    fn test_rotation_matrix_z_90_degrees() {
        let r = rotation_matrix_z(PI / 2.0);
        assert!(Float::abs(r[(0, 0)]) < 1e-10); // cos(90°) ≈ 0
        assert!(Float::abs(r[(0, 1)] + 1.0) < 1e-10); // -sin(90°) ≈ -1
        assert!(Float::abs(r[(1, 0)] - 1.0) < 1e-10); // sin(90°) ≈ 1
        assert!(Float::abs(r[(1, 1)]) < 1e-10); // cos(90°) ≈ 0
        assert!(Float::abs(r[(2, 2)] - 1.0) < 1e-10);
    }

    #[test]
    fn test_rotation_matrix_axis_angle_x_axis() {
        let axis = vec3(1.0, 0.0, 0.0);
        let angle = PI / 2.0;
        let r_axis = rotation_matrix_axis_angle(axis, angle);
        let r_x = rotation_matrix_x(angle);

        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    Float::abs(r_axis[(i, j)] - r_x[(i, j)]) < 1e-10,
                    "Mismatch at ({}, {}): {} vs {}",
                    i,
                    j,
                    r_axis[(i, j)],
                    r_x[(i, j)]
                );
            }
        }
    }

    #[test]
    fn test_rotation_matrix_axis_angle_y_axis() {
        let axis = vec3(0.0, 1.0, 0.0);
        let angle = PI / 2.0;
        let r_axis = rotation_matrix_axis_angle(axis, angle);
        let r_y = rotation_matrix_y(angle);

        for i in 0..4 {
            for j in 0..4 {
                assert!(Float::abs(r_axis[(i, j)] - r_y[(i, j)]) < 1e-10);
            }
        }
    }

    #[test]
    fn test_rotation_matrix_axis_angle_z_axis() {
        let axis = vec3(0.0, 0.0, 1.0);
        let angle = PI / 2.0;
        let r_axis = rotation_matrix_axis_angle(axis, angle);
        let r_z = rotation_matrix_z(angle);

        for i in 0..4 {
            for j in 0..4 {
                assert!(Float::abs(r_axis[(i, j)] - r_z[(i, j)]) < 1e-10);
            }
        }
    }

    #[test]
    fn test_rotation_matrix_axis_angle_arbitrary_axis() {
        use crate::matrix::multiply;
        // Test rotation around axis [1, 1, 1] normalized
        let axis = vec3(1.0, 1.0, 1.0);
        let angle = PI / 3.0; // 60 degrees
        let r = rotation_matrix_axis_angle(axis, angle);

        // The matrix should be orthogonal (R * R^T = I)
        let rt = transpose(&r);
        let identity = Mat4::<f64>::identity();
        let product = multiply(&r, &rt);

        for i in 0..4 {
            for j in 0..4 {
                assert!(Float::abs(product[(i, j)] - identity[(i, j)]) < 1e-10);
            }
        }
    }

    #[test]
    fn test_rotation_matrix_axis_angle_non_normalized() {
        // Should work with non-normalized axes
        let axis = vec3(2.0, 0.0, 0.0); // Non-normalized
        let r = rotation_matrix_axis_angle(axis, PI / 2.0);
        let r_x = rotation_matrix_x(PI / 2.0);

        for i in 0..4 {
            for j in 0..4 {
                assert!(Float::abs(r[(i, j)] - r_x[(i, j)]) < 1e-10);
            }
        }
    }

    #[test]
    #[should_panic(expected = "Axis must be non-zero")]
    fn test_rotation_matrix_axis_angle_zero_axis() {
        let axis = vec3(0.0, 0.0, 0.0);
        rotation_matrix_axis_angle(axis, PI / 2.0);
    }

    #[test]
    fn test_scaling_matrix_3d() {
        let s = scaling_matrix_3d(2.0, 3.0, 4.0);
        assert_eq!(s[(0, 0)], 2.0);
        assert_eq!(s[(1, 1)], 3.0);
        assert_eq!(s[(2, 2)], 4.0);
        assert_eq!(s[(3, 3)], 1.0);
    }

    #[test]
    fn test_scaling_identity() {
        let s = scaling_matrix_3d(1.0, 1.0, 1.0);
        let identity = Mat4::<f64>::identity();
        assert_eq!(s, identity);
    }

    #[test]
    fn test_scaling_uniform() {
        let s = scaling_matrix_3d(2.0, 2.0, 2.0);
        assert_eq!(s[(0, 0)], s[(1, 1)]);
        assert_eq!(s[(1, 1)], s[(2, 2)]);
    }
}
