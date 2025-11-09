//! Vector types and operations.
//!
//! This module provides generic vector types using const generics,
//! along with all fundamental vector operations.
//!
//! # Overview
//!
//! The core type is [`Vector<T, N>`](Vector), a generic N-dimensional vector where:
//! - `T` is the scalar type (must implement [`Scalar`](crate::core::traits::Scalar))
//! - `N` is the dimension (compile-time constant)
//!
//! # Type Aliases
//!
//! For convenience, type aliases are provided for common dimensions:
//! - [`Vec2<T>`] - 2D vector
//! - [`Vec3<T>`] - 3D vector
//! - [`Vec4<T>`] - 4D vector
//!
//! # Examples
//!
//! ```
//! use axiomath::vector::{Vector, Vec3, vec3};
//!
//! // Create vectors using different methods
//! let v1 = Vector::new([1.0, 2.0, 3.0]);
//! let v2: Vec3<f64> = Vector::from([4.0, 5.0, 6.0]);
//! let v3 = vec3(7.0, 8.0, 9.0);
//!
//! // Vectors are Copy
//! let v4 = v1;
//! assert_eq!(v1[0], v4[0]);
//! ```
//!
//! # Operations
//!
//! The module provides:
//! - **Arithmetic**: Component-wise add, subtract, multiply, divide
//! - **Scalar operations**: Scale by scalar, divide by scalar
//! - **Norms**: Manhattan (L1), Euclidean (L2), infinity (L∞), p-norm
//! - **Products**: Dot product, cross product (3D only), triple product
//! - **Projections**: Project onto vector, reject from vector, reflect across plane
//! - **Angles**: Angle between vectors, parallel/perpendicular tests
//!
//! # Mathematical Background
//!
//! Vectors are fundamental mathematical objects representing quantities with
//! both magnitude and direction. This implementation follows standard linear
//! algebra conventions and provides both free functions and operator overloads.

// Module declarations
pub mod vector_type;
pub mod arithmetic;
pub mod scalar_ops;
pub mod norms;
pub mod products;
pub mod projections;
pub mod angles;

// Re-export the main type
pub use vector_type::Vector;

// Re-export arithmetic operations
pub use arithmetic::{add, divide, multiply, negate, reciprocal, subtract};

// Re-export scalar operations
pub use scalar_ops::{divide_by_scalar, scale};

// Re-export norm operations
pub use norms::{
    distance, distance_squared, euclidean_norm, infinity_norm, is_normalized, manhattan_norm,
    magnitude, magnitude_squared, normalize, normalize_or, p_norm,
};

// Re-export product operations
pub use products::{cross_product, dot_product, triple_product};

// Re-export projection operations
pub use projections::{project, reflect, reject};

// Re-export angle operations
pub use angles::{angle_between, is_parallel, is_perpendicular};

// Type aliases for common dimensions
/// A 2-dimensional vector.
///
/// # Examples
///
/// ```
/// use axiomath::vector::Vec2;
///
/// let v: Vec2<f64> = Vec2::new([1.0, 2.0]);
/// assert_eq!(v[0], 1.0);
/// assert_eq!(v[1], 2.0);
/// ```
pub type Vec2<T> = Vector<T, 2>;

/// A 3-dimensional vector.
///
/// # Examples
///
/// ```
/// use axiomath::vector::Vec3;
///
/// let v: Vec3<f64> = Vec3::new([1.0, 2.0, 3.0]);
/// assert_eq!(v[0], 1.0);
/// assert_eq!(v[1], 2.0);
/// assert_eq!(v[2], 3.0);
/// ```
pub type Vec3<T> = Vector<T, 3>;

/// A 4-dimensional vector.
///
/// # Examples
///
/// ```
/// use axiomath::vector::Vec4;
///
/// let v: Vec4<f64> = Vec4::new([1.0, 2.0, 3.0, 4.0]);
/// assert_eq!(v[0], 1.0);
/// assert_eq!(v[3], 4.0);
/// ```
pub type Vec4<T> = Vector<T, 4>;

// Convenience constructors

/// Creates a 2D vector from components.
///
/// # Examples
///
/// ```
/// use axiomath::vector::vec2;
///
/// let v = vec2(1.0, 2.0);
/// assert_eq!(v[0], 1.0);
/// assert_eq!(v[1], 2.0);
/// ```
#[inline]
pub fn vec2<T>(x: T, y: T) -> Vec2<T>
where
    T: crate::core::traits::Scalar,
{
    Vector::new([x, y])
}

/// Creates a 3D vector from components.
///
/// # Examples
///
/// ```
/// use axiomath::vector::vec3;
///
/// let v = vec3(1.0, 2.0, 3.0);
/// assert_eq!(v[0], 1.0);
/// assert_eq!(v[1], 2.0);
/// assert_eq!(v[2], 3.0);
/// ```
#[inline]
pub fn vec3<T>(x: T, y: T, z: T) -> Vec3<T>
where
    T: crate::core::traits::Scalar,
{
    Vector::new([x, y, z])
}

/// Creates a 4D vector from components.
///
/// # Examples
///
/// ```
/// use axiomath::vector::vec4;
///
/// let v = vec4(1.0, 2.0, 3.0, 4.0);
/// assert_eq!(v[0], 1.0);
/// assert_eq!(v[1], 2.0);
/// assert_eq!(v[2], 3.0);
/// assert_eq!(v[3], 4.0);
/// ```
#[inline]
pub fn vec4<T>(x: T, y: T, z: T, w: T) -> Vec4<T>
where
    T: crate::core::traits::Scalar,
{
    Vector::new([x, y, z, w])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_type_aliases() {
        let v2: Vec2<f64> = Vec2::new([1.0, 2.0]);
        assert_eq!(v2.dimension(), 2);

        let v3: Vec3<f64> = Vec3::new([1.0, 2.0, 3.0]);
        assert_eq!(v3.dimension(), 3);

        let v4: Vec4<f64> = Vec4::new([1.0, 2.0, 3.0, 4.0]);
        assert_eq!(v4.dimension(), 4);
    }

    #[test]
    fn test_convenience_constructors() {
        let v2 = vec2(1.0, 2.0);
        assert_eq!(v2[0], 1.0);
        assert_eq!(v2[1], 2.0);

        let v3 = vec3(1.0, 2.0, 3.0);
        assert_eq!(v3[0], 1.0);
        assert_eq!(v3[1], 2.0);
        assert_eq!(v3[2], 3.0);

        let v4 = vec4(1.0, 2.0, 3.0, 4.0);
        assert_eq!(v4[0], 1.0);
        assert_eq!(v4[1], 2.0);
        assert_eq!(v4[2], 3.0);
        assert_eq!(v4[3], 4.0);
    }
}
