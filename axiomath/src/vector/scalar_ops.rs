//! Scalar-vector operations.
//!
//! This module provides operations between vectors and scalars, such as
//! multiplication and division by a scalar value.
//!
//! # Mathematical Background
//!
//! Scalar multiplication scales each component of a vector by a scalar value:
//! ```text
//! (k * v)[i] = k * v[i]  for all i
//! ```
//!
//! This operation preserves direction but changes magnitude.

use crate::core::traits::Scalar;
use crate::vector::Vector;
use core::ops::{Div, Mul};

/// Multiplies a vector by a scalar (scalar multiplication).
///
/// # Mathematical Definition
///
/// ```text
/// b = k * a  where  b[i] = k * a[i]  for all i
/// ```
///
/// # Properties
///
/// - **Associativity**: `k₁ * (k₂ * v) = (k₁ * k₂) * v`
/// - **Distributivity over vector addition**: `k * (v + w) = k*v + k*w`
/// - **Distributivity over scalar addition**: `(k₁ + k₂) * v = k₁*v + k₂*v`
/// - **Identity**: `1 * v = v`
/// - **Zero annihilation**: `0 * v = 0`
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, scale};
///
/// let v = vec3(1.0, 2.0, 3.0);
/// let scaled = scale(&v, 2.5);
///
/// assert_eq!(scaled[0], 2.5);
/// assert_eq!(scaled[1], 5.0);
/// assert_eq!(scaled[2], 7.5);
/// ```
#[inline]
pub fn scale<T, const N: usize>(v: &Vector<T, N>, scalar: T) -> Vector<T, N>
where
    T: Scalar,
{
    v.map(|x| x * scalar)
}

/// Divides a vector by a scalar.
///
/// # Mathematical Definition
///
/// ```text
/// b = a / k  where  b[i] = a[i] / k  for all i
/// ```
///
/// Equivalent to multiplication by the reciprocal: `v / k = (1/k) * v`
///
/// # Warning
///
/// Division by zero will propagate as per the scalar type's behavior
/// (typically infinity or NaN for floating-point types).
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, divide_by_scalar};
///
/// let v = vec3(10.0, 20.0, 30.0);
/// let divided = divide_by_scalar(&v, 5.0);
///
/// assert_eq!(divided[0], 2.0);
/// assert_eq!(divided[1], 4.0);
/// assert_eq!(divided[2], 6.0);
/// ```
#[inline]
pub fn divide_by_scalar<T, const N: usize>(v: &Vector<T, N>, scalar: T) -> Vector<T, N>
where
    T: Scalar,
{
    v.map(|x| x / scalar)
}

// Operator trait implementations

// Vector * Scalar
impl<T, const N: usize> Mul<T> for Vector<T, N>
where
    T: Scalar,
{
    type Output = Self;

    #[inline]
    fn mul(self, scalar: T) -> Self::Output {
        scale(&self, scalar)
    }
}

impl<T, const N: usize> Mul<T> for &Vector<T, N>
where
    T: Scalar,
{
    type Output = Vector<T, N>;

    #[inline]
    fn mul(self, scalar: T) -> Self::Output {
        scale(self, scalar)
    }
}

// Scalar * Vector
impl<const N: usize> Mul<Vector<f32, N>> for f32 {
    type Output = Vector<f32, N>;

    #[inline]
    fn mul(self, v: Vector<f32, N>) -> Self::Output {
        scale(&v, self)
    }
}

impl<const N: usize> Mul<&Vector<f32, N>> for f32 {
    type Output = Vector<f32, N>;

    #[inline]
    fn mul(self, v: &Vector<f32, N>) -> Self::Output {
        scale(v, self)
    }
}

impl<const N: usize> Mul<Vector<f64, N>> for f64 {
    type Output = Vector<f64, N>;

    #[inline]
    fn mul(self, v: Vector<f64, N>) -> Self::Output {
        scale(&v, self)
    }
}

impl<const N: usize> Mul<&Vector<f64, N>> for f64 {
    type Output = Vector<f64, N>;

    #[inline]
    fn mul(self, v: &Vector<f64, N>) -> Self::Output {
        scale(v, self)
    }
}

// Vector / Scalar
impl<T, const N: usize> Div<T> for Vector<T, N>
where
    T: Scalar,
{
    type Output = Self;

    #[inline]
    fn div(self, scalar: T) -> Self::Output {
        divide_by_scalar(&self, scalar)
    }
}

impl<T, const N: usize> Div<T> for &Vector<T, N>
where
    T: Scalar,
{
    type Output = Vector<T, N>;

    #[inline]
    fn div(self, scalar: T) -> Self::Output {
        divide_by_scalar(self, scalar)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vector::vec3;

    #[test]
    fn test_scale() {
        let v = vec3(1.0, 2.0, 3.0);
        let scaled = scale(&v, 2.5);

        assert_eq!(scaled[0], 2.5);
        assert_eq!(scaled[1], 5.0);
        assert_eq!(scaled[2], 7.5);
    }

    #[test]
    fn test_scale_operator_vector_scalar() {
        let v = vec3(1.0, 2.0, 3.0);
        let scaled = v * 2.5;

        assert_eq!(scaled[0], 2.5);
        assert_eq!(scaled[1], 5.0);
        assert_eq!(scaled[2], 7.5);
    }

    #[test]
    fn test_scale_operator_scalar_vector() {
        let v = vec3(1.0, 2.0, 3.0);
        let scaled: Vector<f64, 3> = 2.5 * v;

        assert_eq!(scaled[0], 2.5);
        assert_eq!(scaled[1], 5.0);
        assert_eq!(scaled[2], 7.5);
    }

    #[test]
    fn test_divide_by_scalar() {
        let v = vec3(10.0, 20.0, 30.0);
        let divided = divide_by_scalar(&v, 5.0);

        assert_eq!(divided[0], 2.0);
        assert_eq!(divided[1], 4.0);
        assert_eq!(divided[2], 6.0);
    }

    #[test]
    fn test_divide_operator() {
        let v = vec3(10.0, 20.0, 30.0);
        let divided = v / 5.0;

        assert_eq!(divided[0], 2.0);
        assert_eq!(divided[1], 4.0);
        assert_eq!(divided[2], 6.0);
    }

    #[test]
    fn test_scalar_multiplication_properties() {
        let v = vec3(1.0, 2.0, 3.0);
        let k1 = 2.0;
        let k2 = 3.0;

        // Associativity: k1 * (k2 * v) = (k1 * k2) * v
        let left: Vector<f64, 3> = k1 * (k2 * &v);
        let right: Vector<f64, 3> = (k1 * k2) * &v;
        assert_eq!(left[0], right[0]);
        assert_eq!(left[1], right[1]);
        assert_eq!(left[2], right[2]);

        // Identity: 1 * v = v
        let identity: Vector<f64, 3> = 1.0 * &v;
        assert_eq!(identity[0], v[0]);
        assert_eq!(identity[1], v[1]);
        assert_eq!(identity[2], v[2]);

        // Zero: 0 * v = 0
        let zero: Vector<f64, 3> = 0.0 * &v;
        assert_eq!(zero[0], 0.0);
        assert_eq!(zero[1], 0.0);
        assert_eq!(zero[2], 0.0);
    }

    #[test]
    fn test_distributivity() {
        let v = vec3(1.0, 2.0, 3.0);
        let w = vec3(4.0, 5.0, 6.0);
        let k = 2.0;

        // k * (v + w) = k*v + k*w
        let left: Vector<f64, 3> = k * (&v + &w);
        let right: Vector<f64, 3> = k * &v + k * &w;
        assert_eq!(left[0], right[0]);
        assert_eq!(left[1], right[1]);
        assert_eq!(left[2], right[2]);
    }

    #[test]
    fn test_scalar_division() {
        let v = vec3(2.0, 4.0, 6.0);
        let k = 2.0;

        // v / k = (1/k) * v
        let divided = &v / k;
        let scaled: Vector<f64, 3> = (1.0 / k) * &v;
        assert_eq!(divided[0], scaled[0]);
        assert_eq!(divided[1], scaled[1]);
        assert_eq!(divided[2], scaled[2]);
    }

    #[test]
    fn test_negative_scalar() {
        let v = vec3(1.0, 2.0, 3.0);
        let scaled: Vector<f64, 3> = -2.0 * v;
        assert_eq!(scaled[0], -2.0);
        assert_eq!(scaled[1], -4.0);
        assert_eq!(scaled[2], -6.0);
    }
}
