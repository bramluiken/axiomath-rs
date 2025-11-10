//! 2D transformation matrices using homogeneous coordinates.
//!
//! This module provides functions to create 3x3 transformation matrices
//! for 2D geometric transformations including translation, rotation,
//! scaling, and shearing.

use crate::core::traits::Scalar;
use crate::matrix::Mat3;
use num_traits::{Float, One, Zero};

/// Creates a 2D translation matrix.
///
/// # Mathematical Background
///
/// Translation in 2D homogeneous coordinates is represented by a 3x3 matrix:
///
/// ```text
/// [ 1  0  x ]
/// [ 0  1  y ]
/// [ 0  0  1 ]
/// ```
///
/// When multiplied with a homogeneous point [px, py, 1]^T, this produces
/// [px + x, py + y, 1]^T, effectively translating the point by (x, y).
///
/// # Algorithm
///
/// 1. Create an identity 3x3 matrix
/// 2. Set the translation components in the third column
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
///
/// # Returns
///
/// A 3x3 translation matrix in homogeneous coordinates
///
/// # Example
///
/// ```
/// use axiomath::geometry::transform::translation_matrix_2d;
/// use axiomath::matrix::Mat3;
///
/// let t = translation_matrix_2d(5.0, 3.0);
/// let expected = Mat3::from_array([
///     [1.0, 0.0, 5.0],
///     [0.0, 1.0, 3.0],
///     [0.0, 0.0, 1.0],
/// ]);
/// assert_eq!(t, expected);
/// ```
pub fn translation_matrix_2d<T>(x: T, y: T) -> Mat3<T>
where
    T: Scalar + Float + Zero + One,
{
    Mat3::from_rows([
        [T::one(), T::zero(), x],
        [T::zero(), T::one(), y],
        [T::zero(), T::zero(), T::one()],
    ])
}

/// Creates a 2D rotation matrix.
///
/// # Mathematical Background
///
/// Rotation in 2D homogeneous coordinates is represented by:
///
/// ```text
/// [ cos(θ)  -sin(θ)  0 ]
/// [ sin(θ)   cos(θ)  0 ]
/// [   0        0     1 ]
/// ```
///
/// This rotates points counterclockwise around the origin by angle θ.
/// The rotation preserves distances and angles, making it an orthogonal
/// transformation.
///
/// # Algorithm
///
/// 1. Compute cos(angle) and sin(angle)
/// 2. Construct the rotation matrix using the standard rotation formula
///
/// # Complexity
///
/// - Time: O(1) - involves two trigonometric function calls
/// - Space: O(1)
///
/// # Parameters
///
/// * `angle` - Rotation angle in radians (counterclockwise)
///
/// # Returns
///
/// A 3x3 rotation matrix in homogeneous coordinates
///
/// # Example
///
/// ```
/// use axiomath::geometry::transform::rotation_matrix_2d;
/// use std::f64::consts::PI;
///
/// // 90 degree rotation
/// let r = rotation_matrix_2d(PI / 2.0);
///
/// // cos(90°) ≈ 0, sin(90°) ≈ 1
/// let cos_90 = r[(0, 0)];
/// let sin_90 = r[(1, 0)];
/// assert!((cos_90).abs() < 1e-10);
/// assert!((sin_90 - 1.0).abs() < 1e-10);
/// ```
pub fn rotation_matrix_2d<T>(angle: T) -> Mat3<T>
where
    T: Scalar + Float + Zero + One,
{
    let cos_angle = angle.cos();
    let sin_angle = angle.sin();

    Mat3::from_rows([
        [cos_angle, -sin_angle, T::zero()],
        [sin_angle, cos_angle, T::zero()],
        [T::zero(), T::zero(), T::one()],
    ])
}

/// Creates a 2D scaling matrix with non-uniform scale factors.
///
/// # Mathematical Background
///
/// Scaling in 2D homogeneous coordinates is represented by:
///
/// ```text
/// [ sx  0   0 ]
/// [ 0   sy  0 ]
/// [ 0   0   1 ]
/// ```
///
/// This scales points by sx along the x-axis and sy along the y-axis.
/// When sx = sy, this is uniform scaling. When sx ≠ sy, it's non-uniform.
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
///
/// # Returns
///
/// A 3x3 scaling matrix in homogeneous coordinates
///
/// # Example
///
/// ```
/// use axiomath::geometry::transform::scaling_matrix_2d;
/// use axiomath::matrix::Mat3;
///
/// let s = scaling_matrix_2d(2.0, 3.0);
/// let expected = Mat3::from_array([
///     [2.0, 0.0, 0.0],
///     [0.0, 3.0, 0.0],
///     [0.0, 0.0, 1.0],
/// ]);
/// assert_eq!(s, expected);
/// ```
pub fn scaling_matrix_2d<T>(sx: T, sy: T) -> Mat3<T>
where
    T: Scalar + Float + Zero + One,
{
    Mat3::from_rows([
        [sx, T::zero(), T::zero()],
        [T::zero(), sy, T::zero()],
        [T::zero(), T::zero(), T::one()],
    ])
}

/// Creates a 2D shearing matrix.
///
/// # Mathematical Background
///
/// Shearing (or shear mapping) is a linear transformation that displaces
/// each point in a fixed direction by an amount proportional to its
/// signed distance from a line parallel to that direction.
///
/// The shearing matrix is:
///
/// ```text
/// [ 1   shx  0 ]
/// [ shy  1   0 ]
/// [ 0    0   1 ]
/// ```
///
/// - shx: shear factor along x (y-coordinate affects x displacement)
/// - shy: shear factor along y (x-coordinate affects y displacement)
///
/// # Algorithm
///
/// 1. Create an identity matrix
/// 2. Set the shear components at positions (0,1) and (1,0)
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Parameters
///
/// * `shx` - Shear factor along x-axis
/// * `shy` - Shear factor along y-axis
///
/// # Returns
///
/// A 3x3 shearing matrix in homogeneous coordinates
///
/// # Example
///
/// ```
/// use axiomath::geometry::transform::shearing_matrix_2d;
/// use axiomath::matrix::Mat3;
///
/// let shear = shearing_matrix_2d(0.5, 0.0);
/// let expected = Mat3::from_array([
///     [1.0, 0.5, 0.0],
///     [0.0, 1.0, 0.0],
///     [0.0, 0.0, 1.0],
/// ]);
/// assert_eq!(shear, expected);
/// ```
pub fn shearing_matrix_2d<T>(shx: T, shy: T) -> Mat3<T>
where
    T: Scalar + Float + Zero + One,
{
    Mat3::from_rows([
        [T::one(), shx, T::zero()],
        [shy, T::one(), T::zero()],
        [T::zero(), T::zero(), T::one()],
    ])
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_translation_matrix_2d() {
        let t = translation_matrix_2d(5.0, 3.0);
        let expected = Mat3::from_rows([[1.0, 0.0, 5.0], [0.0, 1.0, 3.0], [0.0, 0.0, 1.0]]);
        assert_eq!(t, expected);
    }

    #[test]
    fn test_translation_identity() {
        let t = translation_matrix_2d(0.0, 0.0);
        let identity = Mat3::<f64>::identity();
        assert_eq!(t, identity);
    }

    #[test]
    fn test_rotation_matrix_2d_zero() {
        let r = rotation_matrix_2d(0.0);
        let identity = Mat3::<f64>::identity();
        for i in 0..3 {
            for j in 0..3 {
                assert!(Float::abs(r[(i, j)] - identity[(i, j)]) < 1e-10);
            }
        }
    }

    #[test]
    fn test_rotation_matrix_2d_90_degrees() {
        let r = rotation_matrix_2d(PI / 2.0);
        assert!(Float::abs(r[(0, 0)]) < 1e-10); // cos(90°) ≈ 0
        assert!(Float::abs(r[(0, 1)] + 1.0) < 1e-10); // -sin(90°) ≈ -1
        assert!(Float::abs(r[(1, 0)] - 1.0) < 1e-10); // sin(90°) ≈ 1
        assert!(Float::abs(r[(1, 1)]) < 1e-10); // cos(90°) ≈ 0
    }

    #[test]
    fn test_rotation_matrix_2d_180_degrees() {
        let r = rotation_matrix_2d(PI);
        assert!(Float::abs(r[(0, 0)] + 1.0) < 1e-10); // cos(180°) ≈ -1
        assert!(Float::abs(r[(0, 1)]) < 1e-10); // -sin(180°) ≈ 0
        assert!(Float::abs(r[(1, 0)]) < 1e-10); // sin(180°) ≈ 0
        assert!(Float::abs(r[(1, 1)] + 1.0) < 1e-10); // cos(180°) ≈ -1
    }

    #[test]
    fn test_scaling_matrix_2d() {
        let s = scaling_matrix_2d(2.0, 3.0);
        let expected = Mat3::from_rows([[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 1.0]]);
        assert_eq!(s, expected);
    }

    #[test]
    fn test_scaling_identity() {
        let s = scaling_matrix_2d(1.0, 1.0);
        let identity = Mat3::<f64>::identity();
        assert_eq!(s, identity);
    }

    #[test]
    fn test_shearing_matrix_2d() {
        let shear = shearing_matrix_2d(0.5, 0.0);
        let expected = Mat3::from_rows([[1.0, 0.5, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        assert_eq!(shear, expected);
    }

    #[test]
    fn test_shearing_identity() {
        let shear = shearing_matrix_2d(0.0, 0.0);
        let identity = Mat3::<f64>::identity();
        assert_eq!(shear, identity);
    }

    #[test]
    fn test_shearing_both_directions() {
        let shear = shearing_matrix_2d(0.5, 0.3);
        assert_eq!(shear[(0, 1)], 0.5);
        assert_eq!(shear[(1, 0)], 0.3);
        assert_eq!(shear[(0, 0)], 1.0);
        assert_eq!(shear[(1, 1)], 1.0);
    }
}
