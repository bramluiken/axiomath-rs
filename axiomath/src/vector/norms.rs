//! Vector norms and magnitude operations.
//!
//! This module provides various norm calculations for vectors, including
//! the Euclidean norm (L2), Manhattan norm (L1), infinity norm (L∞), and
//! general p-norms.
//!
//! # Mathematical Background
//!
//! A norm is a function that assigns a positive length or size to a vector.
//! The most common norms are:
//!
//! - **L1 norm (Manhattan)**: ‖v‖₁ = Σ|vᵢ|
//! - **L2 norm (Euclidean)**: ‖v‖₂ = √(Σvᵢ²)
//! - **L∞ norm (Infinity)**: ‖v‖∞ = max|vᵢ|
//! - **Lp norm**: ‖v‖ₚ = (Σ|vᵢ|ᵖ)^(1/p)

use crate::core::traits::{Real, Scalar};
use crate::vector::Vector;

/// Computes the squared magnitude (L2 norm squared) of a vector.
///
/// This is more efficient than `magnitude()` as it avoids the square root operation.
///
/// # Mathematical Definition
///
/// ```text
/// ‖v‖² = v₀² + v₁² + ... + vₙ²
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, magnitude_squared};
///
/// let v = vec3(3.0, 4.0, 0.0);
/// assert_eq!(magnitude_squared(&v), 25.0);
/// ```
#[inline]
pub fn magnitude_squared<T, const N: usize>(v: &Vector<T, N>) -> T
where
    T: Scalar,
{
    v.iter().map(|&x| x * x).fold(T::zero(), |acc, x| acc + x)
}

/// Computes the magnitude (Euclidean norm, L2 norm) of a vector.
///
/// # Mathematical Definition
///
/// ```text
/// ‖v‖ = √(v₀² + v₁² + ... + vₙ²)
/// ```
///
/// This is the standard Euclidean distance from the origin.
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, magnitude};
///
/// let v = vec3(3.0, 4.0, 0.0);
/// assert_eq!(magnitude(&v), 5.0);
/// ```
#[inline]
pub fn magnitude<T, const N: usize>(v: &Vector<T, N>) -> T
where
    T: Real,
{
    magnitude_squared(v).sqrt()
}

/// Computes the Euclidean norm (same as magnitude, L2 norm).
///
/// This is an alias for `magnitude()` using standard mathematical terminology.
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, euclidean_norm};
///
/// let v = vec3(3.0, 4.0, 0.0);
/// assert_eq!(euclidean_norm(&v), 5.0);
/// ```
#[inline]
pub fn euclidean_norm<T, const N: usize>(v: &Vector<T, N>) -> T
where
    T: Real,
{
    magnitude(v)
}

/// Computes the Manhattan norm (L1 norm, taxicab norm) of a vector.
///
/// # Mathematical Definition
///
/// ```text
/// ‖v‖₁ = |v₀| + |v₁| + ... + |vₙ|
/// ```
///
/// The L1 norm represents the sum of absolute values of components.
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, manhattan_norm};
///
/// let v = vec3(3.0, -4.0, 5.0);
/// assert_eq!(manhattan_norm(&v), 12.0);
/// ```
#[inline]
pub fn manhattan_norm<T, const N: usize>(v: &Vector<T, N>) -> T
where
    T: Real,
{
    v.iter().map(|&x| x.abs()).fold(T::zero(), |acc, x| acc + x)
}

/// Computes the infinity norm (L∞ norm, maximum norm) of a vector.
///
/// # Mathematical Definition
///
/// ```text
/// ‖v‖∞ = max(|v₀|, |v₁|, ..., |vₙ|)
/// ```
///
/// The L∞ norm is the maximum absolute value among all components.
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, infinity_norm};
///
/// let v = vec3(3.0, -7.0, 5.0);
/// assert_eq!(infinity_norm(&v), 7.0);
/// ```
#[inline]
pub fn infinity_norm<T, const N: usize>(v: &Vector<T, N>) -> T
where
    T: Real,
{
    v.iter()
        .map(|&x| x.abs())
        .fold(T::zero(), |max, x| if x > max { x } else { max })
}

/// Computes the p-norm (Lp norm) of a vector.
///
/// # Mathematical Definition
///
/// ```text
/// ‖v‖ₚ = (|v₀|ᵖ + |v₁|ᵖ + ... + |vₙ|ᵖ)^(1/p)
/// ```
///
/// # Special Cases
///
/// - p = 1: Manhattan norm
/// - p = 2: Euclidean norm
/// - p → ∞: Infinity norm
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, p_norm};
///
/// let v = vec3(3.0, 4.0, 0.0);
///
/// // L2 norm (p=2)
/// let l2: f64 = p_norm(&v, 2.0);
/// assert!((l2 - 5.0).abs() < 1e-10);
///
/// // L1 norm (p=1)
/// let l1: f64 = p_norm(&v, 1.0);
/// assert!((l1 - 7.0).abs() < 1e-10);
/// ```
#[inline]
pub fn p_norm<T, const N: usize>(v: &Vector<T, N>, p: T) -> T
where
    T: crate::core::traits::Float,
{
    let sum = v
        .iter()
        .map(|&x| x.abs().powf(p))
        .fold(T::zero(), |acc, x| acc + x);
    sum.powf(T::one() / p)
}

/// Computes the squared distance between two vectors.
///
/// This is more efficient than `distance()` as it avoids the square root operation.
///
/// # Mathematical Definition
///
/// ```text
/// d²(u, v) = ‖u - v‖² = Σ(uᵢ - vᵢ)²
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, distance_squared};
///
/// let a = vec3(1.0, 2.0, 3.0);
/// let b = vec3(4.0, 6.0, 3.0);
/// assert_eq!(distance_squared(&a, &b), 25.0);
/// ```
#[inline]
pub fn distance_squared<T, const N: usize>(a: &Vector<T, N>, b: &Vector<T, N>) -> T
where
    T: Scalar,
{
    a.iter()
        .zip(b.iter())
        .map(|(&ai, &bi)| {
            let diff = ai - bi;
            diff * diff
        })
        .fold(T::zero(), |acc, x| acc + x)
}

/// Computes the Euclidean distance between two vectors.
///
/// # Mathematical Definition
///
/// ```text
/// d(u, v) = ‖u - v‖ = √(Σ(uᵢ - vᵢ)²)
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, distance};
///
/// let a = vec3(1.0, 2.0, 3.0);
/// let b = vec3(4.0, 6.0, 3.0);
/// assert_eq!(distance(&a, &b), 5.0);
/// ```
#[inline]
pub fn distance<T, const N: usize>(a: &Vector<T, N>, b: &Vector<T, N>) -> T
where
    T: Real,
{
    distance_squared(a, b).sqrt()
}

/// Normalizes a vector to unit length.
///
/// Returns `None` if the vector has zero or near-zero magnitude.
///
/// # Mathematical Definition
///
/// ```text
/// normalize(v) = v / ‖v‖
/// ```
///
/// The resulting vector has magnitude 1 and points in the same direction as the original.
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, normalize, magnitude};
///
/// let v = vec3(3.0, 4.0, 0.0);
/// let normalized = normalize(&v).unwrap();
///
/// assert_eq!(normalized[0], 0.6);
/// assert_eq!(normalized[1], 0.8);
/// let mag: f64 = magnitude(&normalized);
/// assert!((mag - 1.0).abs() < 1e-10);
///
/// // Zero vector cannot be normalized
/// let zero = vec3(0.0, 0.0, 0.0);
/// assert!(normalize(&zero).is_none());
/// ```
#[inline]
pub fn normalize<T, const N: usize>(v: &Vector<T, N>) -> Option<Vector<T, N>>
where
    T: crate::core::traits::Float,
{
    let mag = magnitude(v);
    let epsilon = T::epsilon() * T::from(10).unwrap();

    // Check if magnitude is effectively zero
    if mag <= epsilon {
        return None;
    }

    Some(v.map(|x| x / mag))
}

/// Attempts to normalize a vector, returning a default direction if the magnitude is zero.
///
/// If the input vector has zero or near-zero magnitude, returns the provided default
/// instead of None.
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, normalize_or};
///
/// let v = vec3(3.0, 4.0, 0.0);
/// let normalized = normalize_or(&v, vec3(1.0, 0.0, 0.0));
/// assert_eq!(normalized[0], 0.6);
///
/// let zero = vec3(0.0, 0.0, 0.0);
/// let default = vec3(1.0, 0.0, 0.0);
/// let result = normalize_or(&zero, default);
/// assert_eq!(result[0], 1.0);
/// ```
#[inline]
pub fn normalize_or<T, const N: usize>(
    v: &Vector<T, N>,
    default: Vector<T, N>,
) -> Vector<T, N>
where
    T: crate::core::traits::Float,
{
    normalize(v).unwrap_or(default)
}

/// Checks if a vector is normalized (has unit length).
///
/// Uses a tolerance for floating-point comparison.
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, is_normalized};
///
/// let normalized = vec3(0.6, 0.8, 0.0);
/// assert!(is_normalized(&normalized, 1e-10));
///
/// let not_normalized = vec3(3.0, 4.0, 0.0);
/// assert!(!is_normalized(&not_normalized, 1e-10));
/// ```
#[inline]
pub fn is_normalized<T, const N: usize>(v: &Vector<T, N>, tolerance: T) -> bool
where
    T: crate::core::traits::Float,
{
    let mag_sq = magnitude_squared(v);
    let diff = (mag_sq - T::one()).abs();
    diff <= tolerance
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vector::vec3;

    #[test]
    fn test_magnitude_squared() {
        let v = vec3(3.0, 4.0, 0.0);
        assert_eq!(magnitude_squared(&v), 25.0);

        let v2 = vec3(1.0, 2.0, 2.0);
        assert_eq!(magnitude_squared(&v2), 9.0);
    }

    #[test]
    fn test_magnitude() {
        let v = vec3(3.0, 4.0, 0.0);
        assert_eq!(magnitude(&v), 5.0);

        let v2 = vec3(1.0, 2.0, 2.0);
        assert_eq!(magnitude(&v2), 3.0);
    }

    #[test]
    fn test_euclidean_norm() {
        let v = vec3(3.0, 4.0, 0.0);
        assert_eq!(euclidean_norm(&v), 5.0);
        assert_eq!(euclidean_norm(&v), magnitude(&v));
    }

    #[test]
    fn test_manhattan_norm() {
        let v = vec3(3.0, -4.0, 5.0);
        assert_eq!(manhattan_norm(&v), 12.0);

        let v2 = vec3(1.0, 1.0, 1.0);
        assert_eq!(manhattan_norm(&v2), 3.0);
    }

    #[test]
    fn test_infinity_norm() {
        let v = vec3(3.0, -7.0, 5.0);
        assert_eq!(infinity_norm(&v), 7.0);

        let v2 = vec3(-2.0, 1.0, -3.0);
        assert_eq!(infinity_norm(&v2), 3.0);
    }

    #[test]
    fn test_p_norm() {
        let v = vec3(3.0, 4.0, 0.0);

        // L2 norm
        let l2 = p_norm(&v, 2.0);
        assert!((l2 - 5.0).abs() < 1e-10);

        // L1 norm
        let l1 = p_norm(&v, 1.0);
        assert!((l1 - 7.0).abs() < 1e-10);

        // Should match specialized functions
        assert!((l1 - manhattan_norm(&v)).abs() < 1e-10);
        assert!((l2 - euclidean_norm(&v)).abs() < 1e-10);
    }

    #[test]
    fn test_distance_squared() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(4.0, 6.0, 3.0);
        assert_eq!(distance_squared(&a, &b), 25.0);

        // Distance from vector to itself is zero
        assert_eq!(distance_squared(&a, &a), 0.0);
    }

    #[test]
    fn test_distance() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(4.0, 6.0, 3.0);
        assert_eq!(distance(&a, &b), 5.0);

        // Distance is symmetric
        assert_eq!(distance(&a, &b), distance(&b, &a));

        // Distance from vector to itself is zero
        assert_eq!(distance(&a, &a), 0.0);
    }

    #[test]
    fn test_normalize() {
        let v = vec3(3.0, 4.0, 0.0);
        let normalized = normalize(&v).unwrap();

        assert!((normalized[0] - 0.6).abs() < 1e-10);
        assert!((normalized[1] - 0.8).abs() < 1e-10);
        assert_eq!(normalized[2], 0.0);

        // Check that magnitude is 1
        let mag = magnitude(&normalized);
        assert!((mag - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_normalize_zero_vector() {
        let zero = vec3(0.0, 0.0, 0.0);
        assert!(normalize(&zero).is_none());

        let near_zero = vec3(1e-20, 1e-20, 1e-20);
        assert!(normalize(&near_zero).is_none());
    }

    #[test]
    fn test_normalize_or() {
        let v = vec3(3.0, 4.0, 0.0);
        let default = vec3(1.0, 0.0, 0.0);
        let normalized = normalize_or(&v, default);
        assert!((normalized[0] - 0.6).abs() < 1e-10);

        let zero = vec3(0.0, 0.0, 0.0);
        let result = normalize_or(&zero, default);
        assert_eq!(result[0], 1.0);
        assert_eq!(result[1], 0.0);
        assert_eq!(result[2], 0.0);
    }

    #[test]
    fn test_is_normalized() {
        let normalized = vec3(0.6, 0.8, 0.0);
        assert!(is_normalized(&normalized, 1e-10));

        let not_normalized = vec3(3.0, 4.0, 0.0);
        assert!(!is_normalized(&not_normalized, 1e-10));

        let unit_x = vec3(1.0, 0.0, 0.0);
        assert!(is_normalized(&unit_x, 1e-10));
    }

    #[test]
    fn test_norm_properties() {
        let v = vec3(3.0, 4.0, 5.0);

        // Non-negativity: ‖v‖ ≥ 0
        assert!(magnitude(&v) >= 0.0);
        assert!(manhattan_norm(&v) >= 0.0);
        assert!(infinity_norm(&v) >= 0.0);

        // Zero vector has zero norm
        let zero = vec3(0.0, 0.0, 0.0);
        assert_eq!(magnitude(&zero), 0.0);
        assert_eq!(manhattan_norm(&zero), 0.0);
        assert_eq!(infinity_norm(&zero), 0.0);

        // Norm inequality: ‖v‖∞ ≤ ‖v‖₂ ≤ ‖v‖₁
        let inf_norm = infinity_norm(&v);
        let l2_norm = euclidean_norm(&v);
        let l1_norm = manhattan_norm(&v);

        assert!(inf_norm <= l2_norm + 1e-10);
        assert!(l2_norm <= l1_norm + 1e-10);
    }

    #[test]
    fn test_triangle_inequality() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(4.0, 5.0, 6.0);

        // Triangle inequality: ‖a + b‖ ≤ ‖a‖ + ‖b‖
        let sum = &a + &b;
        let left = magnitude(&sum);
        let right = magnitude(&a) + magnitude(&b);

        assert!(left <= right + 1e-10);
    }

    #[test]
    fn test_scalar_multiplication_norm() {
        let v = vec3(3.0, 4.0, 0.0);
        let k = 2.5;

        // ‖k·v‖ = |k|·‖v‖
        let scaled: Vector<f64, 3> = k * v;
        let left = magnitude(&scaled);
        let right = k.abs() * magnitude(&v);

        assert!((left - right).abs() < 1e-10);
    }
}
