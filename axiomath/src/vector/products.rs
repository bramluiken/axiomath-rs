//! Vector product operations.
//!
//! This module provides various vector products including dot product,
//! cross product (3D only), and scalar triple product.
//!
//! # Mathematical Background
//!
//! Vector products are fundamental operations in linear algebra and physics:
//!
//! - **Dot Product (Inner Product)**: Projects one vector onto another, returns a scalar
//! - **Cross Product**: Produces a vector perpendicular to both inputs (3D only)
//! - **Triple Product**: Combines cross and dot products to compute volume

use crate::core::traits::Scalar;
use crate::vector::{Vec3, Vector};

/// Computes the dot product (inner product) of two vectors.
///
/// # Mathematical Definition
///
/// ```text
/// a · b = a₀b₀ + a₁b₁ + ... + aₙbₙ = Σaᵢbᵢ
/// ```
///
/// # Properties
///
/// - **Commutativity**: `a · b = b · a`
/// - **Distributivity**: `a · (b + c) = a · b + a · c`
/// - **Scalar multiplication**: `(ka) · b = k(a · b)`
/// - **Relationship to magnitude**: `a · a = ‖a‖²`
/// - **Relationship to angle**: `a · b = ‖a‖‖b‖cos(θ)`
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, dot_product};
///
/// let a = vec3(1.0, 2.0, 3.0);
/// let b = vec3(4.0, 5.0, 6.0);
/// let dot = dot_product(&a, &b);
///
/// assert_eq!(dot, 32.0); // 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
/// ```
#[inline]
pub fn dot_product<T, const N: usize>(a: &Vector<T, N>, b: &Vector<T, N>) -> T
where
    T: Scalar,
{
    a.iter()
        .zip(b.iter())
        .map(|(&ai, &bi)| ai * bi)
        .fold(T::zero(), |acc, x| acc + x)
}

/// Computes the cross product of two 3D vectors.
///
/// The cross product is only defined for 3-dimensional vectors.
///
/// # Mathematical Definition
///
/// ```text
/// a × b = |  i    j    k  |
///         | a₀   a₁   a₂ |
///         | b₀   b₁   b₂ |
///
/// = (a₁b₂ - a₂b₁)i - (a₀b₂ - a₂b₀)j + (a₀b₁ - a₁b₀)k
/// = [a₁b₂ - a₂b₁, a₂b₀ - a₀b₂, a₀b₁ - a₁b₀]
/// ```
///
/// # Properties
///
/// - **Anticommutativity**: `a × b = -(b × a)`
/// - **Distributivity**: `a × (b + c) = a × b + a × c`
/// - **Scalar multiplication**: `(ka) × b = k(a × b)`
/// - **Perpendicularity**: `(a × b) · a = 0` and `(a × b) · b = 0`
/// - **Magnitude**: `‖a × b‖ = ‖a‖‖b‖sin(θ)`
/// - **Parallel vectors**: If `a ∥ b`, then `a × b = 0`
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, cross_product};
///
/// let i = vec3(1.0, 0.0, 0.0);
/// let j = vec3(0.0, 1.0, 0.0);
/// let k = cross_product(&i, &j);
///
/// assert_eq!(k[0], 0.0);
/// assert_eq!(k[1], 0.0);
/// assert_eq!(k[2], 1.0); // i × j = k
/// ```
#[inline]
pub fn cross_product<T>(a: &Vec3<T>, b: &Vec3<T>) -> Vec3<T>
where
    T: Scalar,
{
    Vector::new([
        a[1] * b[2] - a[2] * b[1], // i component
        a[2] * b[0] - a[0] * b[2], // j component
        a[0] * b[1] - a[1] * b[0], // k component
    ])
}

/// Computes the scalar triple product of three 3D vectors.
///
/// The scalar triple product represents the signed volume of the parallelepiped
/// formed by the three vectors.
///
/// # Mathematical Definition
///
/// ```text
/// a · (b × c) = det([a b c])
/// ```
///
/// This equals the determinant of the 3×3 matrix formed by the three vectors.
///
/// # Properties
///
/// - **Cyclic permutation**: `a · (b × c) = b · (c × a) = c · (a × b)`
/// - **Anticyclic permutation**: `a · (b × c) = -a · (c × b)`
/// - **Geometric meaning**: Absolute value gives volume of parallelepiped
/// - **Coplanarity test**: If result is zero, vectors are coplanar
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, triple_product};
///
/// // Unit cube vectors
/// let i = vec3(1.0, 0.0, 0.0);
/// let j = vec3(0.0, 1.0, 0.0);
/// let k = vec3(0.0, 0.0, 1.0);
///
/// let volume = triple_product(&i, &j, &k);
/// assert_eq!(volume, 1.0); // Volume of unit cube
/// ```
#[inline]
pub fn triple_product<T>(a: &Vec3<T>, b: &Vec3<T>, c: &Vec3<T>) -> T
where
    T: Scalar,
{
    let cross = cross_product(b, c);
    dot_product(a, &cross)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vector::{magnitude, vec3};

    #[test]
    fn test_dot_product() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(4.0, 5.0, 6.0);
        let dot = dot_product(&a, &b);

        assert_eq!(dot, 32.0); // 1*4 + 2*5 + 3*6 = 32
    }

    #[test]
    fn test_dot_product_commutative() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(4.0, 5.0, 6.0);

        // a · b = b · a
        assert_eq!(dot_product(&a, &b), dot_product(&b, &a));
    }

    #[test]
    fn test_dot_product_distributive() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(4.0, 5.0, 6.0);
        let c = vec3(7.0, 8.0, 9.0);

        // a · (b + c) = a · b + a · c
        let left = dot_product(&a, &(&b + &c));
        let right = dot_product(&a, &b) + dot_product(&a, &c);

        assert!((left - right).abs() < 1e-10);
    }

    #[test]
    fn test_dot_product_with_self() {
        let v = vec3(3.0, 4.0, 0.0);

        // v · v = ‖v‖²
        let dot = dot_product(&v, &v);
        let mag_sq = magnitude(&v) * magnitude(&v);

        assert!((dot - mag_sq).abs() < 1e-10);
    }

    #[test]
    fn test_dot_product_orthogonal() {
        let i = vec3(1.0, 0.0, 0.0);
        let j = vec3(0.0, 1.0, 0.0);

        // Orthogonal vectors have zero dot product
        assert_eq!(dot_product(&i, &j), 0.0);
    }

    #[test]
    fn test_cross_product() {
        let i = vec3(1.0, 0.0, 0.0);
        let j = vec3(0.0, 1.0, 0.0);
        let k = cross_product(&i, &j);

        assert_eq!(k[0], 0.0);
        assert_eq!(k[1], 0.0);
        assert_eq!(k[2], 1.0); // i × j = k
    }

    #[test]
    fn test_cross_product_anticommutative() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(4.0, 5.0, 6.0);

        // a × b = -(b × a)
        let ab = cross_product(&a, &b);
        let ba = cross_product(&b, &a);

        assert_eq!(ab[0], -ba[0]);
        assert_eq!(ab[1], -ba[1]);
        assert_eq!(ab[2], -ba[2]);
    }

    #[test]
    fn test_cross_product_perpendicular() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(4.0, 5.0, 6.0);
        let cross = cross_product(&a, &b);

        // (a × b) · a = 0
        assert!((dot_product(&cross, &a)).abs() < 1e-10);

        // (a × b) · b = 0
        assert!((dot_product(&cross, &b)).abs() < 1e-10);
    }

    #[test]
    fn test_cross_product_parallel() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(2.0, 4.0, 6.0); // b = 2a

        // Parallel vectors have zero cross product
        let cross = cross_product(&a, &b);
        assert!(cross[0].abs() < 1e-10);
        assert!(cross[1].abs() < 1e-10);
        assert!(cross[2].abs() < 1e-10);
    }

    #[test]
    fn test_cross_product_magnitude() {
        let a = vec3(3.0, 0.0, 0.0);
        let b = vec3(0.0, 4.0, 0.0);
        let cross = cross_product(&a, &b);

        // ‖a × b‖ = ‖a‖‖b‖sin(90°) = ‖a‖‖b‖ = 3*4 = 12
        let cross_mag = magnitude(&cross);
        assert!((cross_mag - 12.0).abs() < 1e-10);
    }

    #[test]
    fn test_triple_product() {
        // Unit cube vectors
        let i = vec3(1.0, 0.0, 0.0);
        let j = vec3(0.0, 1.0, 0.0);
        let k = vec3(0.0, 0.0, 1.0);

        let volume = triple_product(&i, &j, &k);
        assert_eq!(volume, 1.0); // Volume of unit cube
    }

    #[test]
    fn test_triple_product_cyclic() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(4.0, 5.0, 6.0);
        let c = vec3(7.0, 8.0, 9.0);

        // a · (b × c) = b · (c × a) = c · (a × b)
        let abc = triple_product(&a, &b, &c);
        let bca = triple_product(&b, &c, &a);
        let cab = triple_product(&c, &a, &b);

        assert!((abc - bca).abs() < 1e-10);
        assert!((abc - cab).abs() < 1e-10);
    }

    #[test]
    fn test_triple_product_anticyclic() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(4.0, 5.0, 6.0);
        let c = vec3(7.0, 8.0, 9.0);

        // a · (b × c) = -a · (c × b)
        let abc = triple_product(&a, &b, &c);
        let acb = triple_product(&a, &c, &b);

        assert!((abc + acb).abs() < 1e-10);
    }

    #[test]
    fn test_triple_product_coplanar() {
        // Three coplanar vectors (in xy-plane)
        let a = vec3(1.0, 0.0, 0.0);
        let b = vec3(0.0, 1.0, 0.0);
        let c = vec3(1.0, 1.0, 0.0);

        let volume = triple_product(&a, &b, &c);
        assert_eq!(volume, 0.0); // Coplanar vectors have zero volume
    }

    #[test]
    fn test_cross_product_basis_vectors() {
        let i = vec3(1.0, 0.0, 0.0);
        let j = vec3(0.0, 1.0, 0.0);
        let k = vec3(0.0, 0.0, 1.0);

        // i × j = k
        let ij = cross_product(&i, &j);
        assert_eq!(ij[0], k[0]);
        assert_eq!(ij[1], k[1]);
        assert_eq!(ij[2], k[2]);

        // j × k = i
        let jk = cross_product(&j, &k);
        assert_eq!(jk[0], i[0]);
        assert_eq!(jk[1], i[1]);
        assert_eq!(jk[2], i[2]);

        // k × i = j
        let ki = cross_product(&k, &i);
        assert_eq!(ki[0], j[0]);
        assert_eq!(ki[1], j[1]);
        assert_eq!(ki[2], j[2]);
    }

    #[test]
    fn test_lagrange_identity() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(4.0, 5.0, 6.0);

        // ‖a × b‖² = ‖a‖²‖b‖² - (a·b)²
        let cross = cross_product(&a, &b);
        let left = magnitude(&cross) * magnitude(&cross);
        let right = magnitude(&a) * magnitude(&a) * magnitude(&b) * magnitude(&b)
            - dot_product(&a, &b) * dot_product(&a, &b);

        assert!((left - right).abs() < 1e-10);
    }
}
