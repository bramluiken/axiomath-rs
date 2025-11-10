//! Homogeneous coordinate utilities for applying transformations.
//!
//! This module provides helper functions to apply transformation matrices
//! to points and vectors using homogeneous coordinates.

use crate::core::traits::Scalar;
use crate::matrix::{Mat3, Mat4};
use crate::vector::{vec2, vec3, Vec2, Vec3};
use num_traits::{Float, One, Zero};

/// Applies a 2D transformation matrix to a 2D point.
///
/// # Mathematical Background
///
/// In 2D homogeneous coordinates, a point [x, y] is represented as [x, y, 1].
/// A transformation matrix M (3x3) is applied via matrix-vector multiplication:
///
/// ```text
/// [x']   [M00 M01 M02]   [x]
/// [y'] = [M10 M11 M12] * [y]
/// [1 ]   [M20 M21 M22]   [1]
/// ```
///
/// After multiplication, we have [x', y', w]. The result is then converted
/// back to Cartesian coordinates by dividing by w: (x'/w, y'/w).
///
/// For affine transformations (translation, rotation, scaling, shearing),
/// w will always be 1, so the division is unnecessary. However, this
/// implementation supports projective transformations where w ≠ 1.
///
/// # Algorithm
///
/// 1. Extend the 2D point to homogeneous coordinates [x, y, 1]
/// 2. Multiply the transformation matrix by the homogeneous point
/// 3. Divide by the w component to get Cartesian coordinates
///
/// # Complexity
///
/// - Time: O(1) - constant number of operations
/// - Space: O(1)
///
/// # Parameters
///
/// * `point` - The 2D point to transform
/// * `matrix` - The 3x3 transformation matrix
///
/// # Returns
///
/// The transformed 2D point
///
/// # Example
///
/// ```
/// use axiomath::geometry::transform::{apply_transform_2d, translation_matrix_2d};
/// use axiomath::vector::vec2;
///
/// let point = vec2(1.0, 2.0);
/// let transform = translation_matrix_2d(5.0, 3.0);
/// let result = apply_transform_2d(point, &transform);
///
/// assert!((result[0] - 6.0).abs() < 1e-10);
/// assert!((result[1] - 5.0).abs() < 1e-10);
/// ```
pub fn apply_transform_2d<T>(point: Vec2<T>, matrix: &Mat3<T>) -> Vec2<T>
where
    T: Scalar + Float + Zero + One,
{
    // Convert to homogeneous coordinates [x, y, 1]
    let x = point[0];
    let y = point[1];
    let w = T::one();

    // Apply transformation: result = M * [x, y, 1]^T
    let x_prime = matrix[(0, 0)] * x + matrix[(0, 1)] * y + matrix[(0, 2)] * w;
    let y_prime = matrix[(1, 0)] * x + matrix[(1, 1)] * y + matrix[(1, 2)] * w;
    let w_prime = matrix[(2, 0)] * x + matrix[(2, 1)] * y + matrix[(2, 2)] * w;

    // Convert back to Cartesian coordinates
    vec2(x_prime / w_prime, y_prime / w_prime)
}

/// Applies a 3D transformation matrix to a 3D point.
///
/// # Mathematical Background
///
/// In 3D homogeneous coordinates, a point [x, y, z] is represented as [x, y, z, 1].
/// A transformation matrix M (4x4) is applied via matrix-vector multiplication:
///
/// ```text
/// [x']   [M00 M01 M02 M03]   [x]
/// [y']   [M10 M11 M12 M13]   [y]
/// [z'] = [M20 M21 M22 M23] * [z]
/// [w']   [M30 M31 M32 M33]   [1]
/// ```
///
/// After multiplication, we have [x', y', z', w']. The result is then converted
/// back to Cartesian coordinates by dividing by w': (x'/w', y'/w', z'/w').
///
/// For affine transformations, w' will always be 1. For projective transformations
/// (like perspective projection), w' ≠ 1 and the division is essential.
///
/// # Algorithm
///
/// 1. Extend the 3D point to homogeneous coordinates [x, y, z, 1]
/// 2. Multiply the transformation matrix by the homogeneous point
/// 3. Divide by the w component to get Cartesian coordinates
///
/// # Complexity
///
/// - Time: O(1) - constant number of operations
/// - Space: O(1)
///
/// # Parameters
///
/// * `point` - The 3D point to transform
/// * `matrix` - The 4x4 transformation matrix
///
/// # Returns
///
/// The transformed 3D point
///
/// # Example
///
/// ```
/// use axiomath::geometry::transform::{apply_transform_3d, translation_matrix_3d};
/// use axiomath::vector::vec3;
///
/// let point = vec3(1.0, 2.0, 3.0);
/// let transform = translation_matrix_3d(5.0, 3.0, 2.0);
/// let result = apply_transform_3d(point, &transform);
///
/// assert!((result[0] - 6.0).abs() < 1e-10);
/// assert!((result[1] - 5.0).abs() < 1e-10);
/// assert!((result[2] - 5.0).abs() < 1e-10);
/// ```
pub fn apply_transform_3d<T>(point: Vec3<T>, matrix: &Mat4<T>) -> Vec3<T>
where
    T: Scalar + Float + Zero + One,
{
    // Convert to homogeneous coordinates [x, y, z, 1]
    let x = point[0];
    let y = point[1];
    let z = point[2];
    let w = T::one();

    // Apply transformation: result = M * [x, y, z, 1]^T
    let x_prime = matrix[(0, 0)] * x
        + matrix[(0, 1)] * y
        + matrix[(0, 2)] * z
        + matrix[(0, 3)] * w;
    let y_prime = matrix[(1, 0)] * x
        + matrix[(1, 1)] * y
        + matrix[(1, 2)] * z
        + matrix[(1, 3)] * w;
    let z_prime = matrix[(2, 0)] * x
        + matrix[(2, 1)] * y
        + matrix[(2, 2)] * z
        + matrix[(2, 3)] * w;
    let w_prime = matrix[(3, 0)] * x
        + matrix[(3, 1)] * y
        + matrix[(3, 2)] * z
        + matrix[(3, 3)] * w;

    // Convert back to Cartesian coordinates
    vec3(
        x_prime / w_prime,
        y_prime / w_prime,
        z_prime / w_prime,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::transform::{
        rotation_matrix_2d, rotation_matrix_z, scaling_matrix_2d, scaling_matrix_3d,
        translation_matrix_2d, translation_matrix_3d,
    };
    use num_traits::Float;
    use std::f64::consts::PI;

    #[test]
    fn test_apply_transform_2d_translation() {
        let point = vec2(1.0, 2.0);
        let transform = translation_matrix_2d(5.0, 3.0);
        let result = apply_transform_2d(point, &transform);

        assert!(Float::abs(result[0] - 6.0) < 1e-10);
        assert!(Float::abs(result[1] - 5.0) < 1e-10);
    }

    #[test]
    fn test_apply_transform_2d_rotation() {
        let point = vec2(1.0, 0.0);
        let transform = rotation_matrix_2d(PI / 2.0); // 90 degrees
        let result = apply_transform_2d(point, &transform);

        // After 90° rotation, (1, 0) should become approximately (0, 1)
        assert!(Float::abs(result[0]) < 1e-10);
        assert!(Float::abs(result[1] - 1.0) < 1e-10);
    }

    #[test]
    fn test_apply_transform_2d_scaling() {
        let point = vec2(2.0, 3.0);
        let transform = scaling_matrix_2d(2.0, 3.0);
        let result = apply_transform_2d(point, &transform);

        assert!(Float::abs(result[0] - 4.0) < 1e-10);
        assert!(Float::abs(result[1] - 9.0) < 1e-10);
    }

    #[test]
    fn test_apply_transform_2d_identity() {
        let point = vec2(5.0, 7.0);
        let identity = Mat3::<f64>::identity();
        let result = apply_transform_2d(point, &identity);

        assert!(Float::abs(result[0] - point[0]) < 1e-10);
        assert!(Float::abs(result[1] - point[1]) < 1e-10);
    }

    #[test]
    fn test_apply_transform_2d_composition() {
        use crate::matrix::multiply;
        // Test that applying T(5,3) * R(90°) is equivalent to
        // applying R(90°) then T(5,3)
        let point = vec2(1.0, 0.0);

        let rotation = rotation_matrix_2d(PI / 2.0);
        let translation = translation_matrix_2d(5.0, 3.0);

        // First rotate, then translate
        let rotated = apply_transform_2d(point, &rotation);
        let result1 = apply_transform_2d(rotated, &translation);

        // Apply composed transformation
        let composed = multiply(&translation, &rotation);
        let result2 = apply_transform_2d(point, &composed);

        assert!(Float::abs(result1[0] - result2[0]) < 1e-10);
        assert!(Float::abs(result1[1] - result2[1]) < 1e-10);
    }

    #[test]
    fn test_apply_transform_3d_translation() {
        let point = vec3(1.0, 2.0, 3.0);
        let transform = translation_matrix_3d(5.0, 3.0, 2.0);
        let result = apply_transform_3d(point, &transform);

        assert!(Float::abs(result[0] - 6.0) < 1e-10);
        assert!(Float::abs(result[1] - 5.0) < 1e-10);
        assert!(Float::abs(result[2] - 5.0) < 1e-10);
    }

    #[test]
    fn test_apply_transform_3d_rotation() {
        let point = vec3(1.0, 0.0, 0.0);
        let transform = rotation_matrix_z(PI / 2.0); // 90 degrees around Z
        let result = apply_transform_3d(point, &transform);

        // After 90° rotation around Z, (1, 0, 0) should become approximately (0, 1, 0)
        assert!(Float::abs(result[0]) < 1e-10);
        assert!(Float::abs(result[1] - 1.0) < 1e-10);
        assert!(Float::abs(result[2]) < 1e-10);
    }

    #[test]
    fn test_apply_transform_3d_scaling() {
        let point = vec3(2.0, 3.0, 4.0);
        let transform = scaling_matrix_3d(2.0, 3.0, 4.0);
        let result = apply_transform_3d(point, &transform);

        assert!(Float::abs(result[0] - 4.0) < 1e-10);
        assert!(Float::abs(result[1] - 9.0) < 1e-10);
        assert!(Float::abs(result[2] - 16.0) < 1e-10);
    }

    #[test]
    fn test_apply_transform_3d_identity() {
        let point = vec3(5.0, 7.0, 9.0);
        let identity = Mat4::<f64>::identity();
        let result = apply_transform_3d(point, &identity);

        assert!(Float::abs(result[0] - point[0]) < 1e-10);
        assert!(Float::abs(result[1] - point[1]) < 1e-10);
        assert!(Float::abs(result[2] - point[2]) < 1e-10);
    }

    #[test]
    fn test_apply_transform_3d_composition() {
        use crate::matrix::multiply;
        // Test that applying transformations in sequence gives same result
        // as applying the composed matrix
        let point = vec3(1.0, 0.0, 0.0);

        let rotation = rotation_matrix_z(PI / 2.0);
        let translation = translation_matrix_3d(5.0, 3.0, 2.0);

        // First rotate, then translate
        let rotated = apply_transform_3d(point, &rotation);
        let result1 = apply_transform_3d(rotated, &translation);

        // Apply composed transformation
        let composed = multiply(&translation, &rotation);
        let result2 = apply_transform_3d(point, &composed);

        assert!(Float::abs(result1[0] - result2[0]) < 1e-10);
        assert!(Float::abs(result1[1] - result2[1]) < 1e-10);
        assert!(Float::abs(result1[2] - result2[2]) < 1e-10);
    }
}
