//! Component-wise arithmetic operations for vectors.
//!
//! This module provides component-wise (element-wise or Hadamard) arithmetic operations:
//! addition, subtraction, multiplication, and division. These operations work on
//! corresponding components of two vectors.
//!
//! # Mathematical Background
//!
//! Component-wise operations apply the operation to each pair of corresponding elements:
//! - Addition: `(a + b)[i] = a[i] + b[i]`
//! - Subtraction: `(a - b)[i] = a[i] - b[i]`
//! - Multiplication (Hadamard): `(a ⊙ b)[i] = a[i] * b[i]`
//! - Division: `(a / b)[i] = a[i] / b[i]`

use crate::core::traits::Scalar;
use crate::vector::Vector;
use core::ops::{Add, Div, Mul, Neg, Sub};

/// Adds two vectors component-wise.
///
/// # Mathematical Definition
///
/// ```text
/// c = a + b  where  c[i] = a[i] + b[i]  for all i
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, add};
///
/// let a = vec3(1.0, 2.0, 3.0);
/// let b = vec3(4.0, 5.0, 6.0);
/// let c = add(&a, &b);
///
/// assert_eq!(c[0], 5.0);
/// assert_eq!(c[1], 7.0);
/// assert_eq!(c[2], 9.0);
/// ```
#[inline]
pub fn add<T, const N: usize>(a: &Vector<T, N>, b: &Vector<T, N>) -> Vector<T, N>
where
    T: Scalar,
{
    a.zip_with(b, |x, y| x + y)
}

/// Subtracts two vectors component-wise.
///
/// # Mathematical Definition
///
/// ```text
/// c = a - b  where  c[i] = a[i] - b[i]  for all i
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, subtract};
///
/// let a = vec3(5.0, 7.0, 9.0);
/// let b = vec3(1.0, 2.0, 3.0);
/// let c = subtract(&a, &b);
///
/// assert_eq!(c[0], 4.0);
/// assert_eq!(c[1], 5.0);
/// assert_eq!(c[2], 6.0);
/// ```
#[inline]
pub fn subtract<T, const N: usize>(a: &Vector<T, N>, b: &Vector<T, N>) -> Vector<T, N>
where
    T: Scalar,
{
    a.zip_with(b, |x, y| x - y)
}

/// Multiplies two vectors component-wise (Hadamard product).
///
/// # Mathematical Definition
///
/// The Hadamard product (also known as element-wise or Schur product):
/// ```text
/// c = a ⊙ b  where  c[i] = a[i] * b[i]  for all i
/// ```
///
/// Note: This is different from the dot product which returns a scalar.
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, multiply};
///
/// let a = vec3(2.0, 3.0, 4.0);
/// let b = vec3(5.0, 6.0, 7.0);
/// let c = multiply(&a, &b);
///
/// assert_eq!(c[0], 10.0);
/// assert_eq!(c[1], 18.0);
/// assert_eq!(c[2], 28.0);
/// ```
#[inline]
pub fn multiply<T, const N: usize>(a: &Vector<T, N>, b: &Vector<T, N>) -> Vector<T, N>
where
    T: Scalar,
{
    a.zip_with(b, |x, y| x * y)
}

/// Divides two vectors component-wise.
///
/// # Mathematical Definition
///
/// ```text
/// c = a / b  where  c[i] = a[i] / b[i]  for all i
/// ```
///
/// # Warning
///
/// Division by zero will propagate as per the scalar type's behavior
/// (typically infinity or NaN for floating-point types).
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, divide};
///
/// let a = vec3(10.0, 18.0, 28.0);
/// let b = vec3(2.0, 3.0, 4.0);
/// let c = divide(&a, &b);
///
/// assert_eq!(c[0], 5.0);
/// assert_eq!(c[1], 6.0);
/// assert_eq!(c[2], 7.0);
/// ```
#[inline]
pub fn divide<T, const N: usize>(a: &Vector<T, N>, b: &Vector<T, N>) -> Vector<T, N>
where
    T: Scalar,
{
    a.zip_with(b, |x, y| x / y)
}

/// Negates a vector (returns its additive inverse).
///
/// # Mathematical Definition
///
/// ```text
/// b = -a  where  b[i] = -a[i]  for all i
/// ```
///
/// # Properties
///
/// - `a + (-a) = 0` (additive inverse)
/// - `-(-a) = a` (double negation)
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, negate};
///
/// let a = vec3(1.0, -2.0, 3.0);
/// let b = negate(&a);
///
/// assert_eq!(b[0], -1.0);
/// assert_eq!(b[1], 2.0);
/// assert_eq!(b[2], -3.0);
/// ```
#[inline]
pub fn negate<T, const N: usize>(a: &Vector<T, N>) -> Vector<T, N>
where
    T: Scalar + core::ops::Neg<Output = T>,
{
    a.map(|x| -x)
}

/// Computes the component-wise reciprocal of a vector.
///
/// # Mathematical Definition
///
/// ```text
/// b = 1/a  where  b[i] = 1/a[i]  for all i
/// ```
///
/// # Warning
///
/// Reciprocal of zero will produce infinity or NaN for floating-point types.
///
/// # Examples
///
/// ```
/// use axiomath::vector::{vec3, reciprocal};
///
/// let a = vec3(2.0, 4.0, 5.0);
/// let b = reciprocal(&a);
///
/// assert_eq!(b[0], 0.5);
/// assert_eq!(b[1], 0.25);
/// assert_eq!(b[2], 0.2);
/// ```
#[inline]
pub fn reciprocal<T, const N: usize>(a: &Vector<T, N>) -> Vector<T, N>
where
    T: Scalar,
{
    a.map(|x| T::one() / x)
}

// Operator trait implementations

impl<T, const N: usize> Add for Vector<T, N>
where
    T: Scalar,
{
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        add(&self, &rhs)
    }
}

impl<T, const N: usize> Add<&Vector<T, N>> for Vector<T, N>
where
    T: Scalar,
{
    type Output = Self;

    #[inline]
    fn add(self, rhs: &Vector<T, N>) -> Self::Output {
        add(&self, rhs)
    }
}

impl<T, const N: usize> Add<Vector<T, N>> for &Vector<T, N>
where
    T: Scalar,
{
    type Output = Vector<T, N>;

    #[inline]
    fn add(self, rhs: Vector<T, N>) -> Self::Output {
        add(self, &rhs)
    }
}

impl<T, const N: usize> Add<&Vector<T, N>> for &Vector<T, N>
where
    T: Scalar,
{
    type Output = Vector<T, N>;

    #[inline]
    fn add(self, rhs: &Vector<T, N>) -> Self::Output {
        add(self, rhs)
    }
}

impl<T, const N: usize> Sub for Vector<T, N>
where
    T: Scalar,
{
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        subtract(&self, &rhs)
    }
}

impl<T, const N: usize> Sub<&Vector<T, N>> for Vector<T, N>
where
    T: Scalar,
{
    type Output = Self;

    #[inline]
    fn sub(self, rhs: &Vector<T, N>) -> Self::Output {
        subtract(&self, rhs)
    }
}

impl<T, const N: usize> Sub<Vector<T, N>> for &Vector<T, N>
where
    T: Scalar,
{
    type Output = Vector<T, N>;

    #[inline]
    fn sub(self, rhs: Vector<T, N>) -> Self::Output {
        subtract(self, &rhs)
    }
}

impl<T, const N: usize> Sub<&Vector<T, N>> for &Vector<T, N>
where
    T: Scalar,
{
    type Output = Vector<T, N>;

    #[inline]
    fn sub(self, rhs: &Vector<T, N>) -> Self::Output {
        subtract(self, rhs)
    }
}

impl<T, const N: usize> Mul for Vector<T, N>
where
    T: Scalar,
{
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        multiply(&self, &rhs)
    }
}

impl<T, const N: usize> Mul<&Vector<T, N>> for Vector<T, N>
where
    T: Scalar,
{
    type Output = Self;

    #[inline]
    fn mul(self, rhs: &Vector<T, N>) -> Self::Output {
        multiply(&self, rhs)
    }
}

impl<T, const N: usize> Mul<Vector<T, N>> for &Vector<T, N>
where
    T: Scalar,
{
    type Output = Vector<T, N>;

    #[inline]
    fn mul(self, rhs: Vector<T, N>) -> Self::Output {
        multiply(self, &rhs)
    }
}

impl<T, const N: usize> Mul<&Vector<T, N>> for &Vector<T, N>
where
    T: Scalar,
{
    type Output = Vector<T, N>;

    #[inline]
    fn mul(self, rhs: &Vector<T, N>) -> Self::Output {
        multiply(self, rhs)
    }
}

impl<T, const N: usize> Div for Vector<T, N>
where
    T: Scalar,
{
    type Output = Self;

    #[inline]
    fn div(self, rhs: Self) -> Self::Output {
        divide(&self, &rhs)
    }
}

impl<T, const N: usize> Div<&Vector<T, N>> for Vector<T, N>
where
    T: Scalar,
{
    type Output = Self;

    #[inline]
    fn div(self, rhs: &Vector<T, N>) -> Self::Output {
        divide(&self, rhs)
    }
}

impl<T, const N: usize> Div<Vector<T, N>> for &Vector<T, N>
where
    T: Scalar,
{
    type Output = Vector<T, N>;

    #[inline]
    fn div(self, rhs: Vector<T, N>) -> Self::Output {
        divide(self, &rhs)
    }
}

impl<T, const N: usize> Div<&Vector<T, N>> for &Vector<T, N>
where
    T: Scalar,
{
    type Output = Vector<T, N>;

    #[inline]
    fn div(self, rhs: &Vector<T, N>) -> Self::Output {
        divide(self, rhs)
    }
}

impl<T, const N: usize> Neg for Vector<T, N>
where
    T: Scalar + core::ops::Neg<Output = T>,
{
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        negate(&self)
    }
}

impl<T, const N: usize> Neg for &Vector<T, N>
where
    T: Scalar + core::ops::Neg<Output = T>,
{
    type Output = Vector<T, N>;

    #[inline]
    fn neg(self) -> Self::Output {
        negate(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vector::vec3;

    #[test]
    fn test_add() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(4.0, 5.0, 6.0);
        let c = add(&a, &b);

        assert_eq!(c[0], 5.0);
        assert_eq!(c[1], 7.0);
        assert_eq!(c[2], 9.0);
    }

    #[test]
    fn test_add_operator() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(4.0, 5.0, 6.0);

        // Test all combinations
        let c1 = a + b;
        assert_eq!(c1[0], 5.0);

        let c2 = a + &b;
        assert_eq!(c2[0], 5.0);

        let c3 = &a + b;
        assert_eq!(c3[0], 5.0);

        let c4 = &a + &b;
        assert_eq!(c4[0], 5.0);
    }

    #[test]
    fn test_subtract() {
        let a = vec3(5.0, 7.0, 9.0);
        let b = vec3(1.0, 2.0, 3.0);
        let c = subtract(&a, &b);

        assert_eq!(c[0], 4.0);
        assert_eq!(c[1], 5.0);
        assert_eq!(c[2], 6.0);
    }

    #[test]
    fn test_subtract_operator() {
        let a = vec3(5.0, 7.0, 9.0);
        let b = vec3(1.0, 2.0, 3.0);
        let c = a - b;

        assert_eq!(c[0], 4.0);
        assert_eq!(c[1], 5.0);
        assert_eq!(c[2], 6.0);
    }

    #[test]
    fn test_multiply() {
        let a = vec3(2.0, 3.0, 4.0);
        let b = vec3(5.0, 6.0, 7.0);
        let c = multiply(&a, &b);

        assert_eq!(c[0], 10.0);
        assert_eq!(c[1], 18.0);
        assert_eq!(c[2], 28.0);
    }

    #[test]
    fn test_multiply_operator() {
        let a = vec3(2.0, 3.0, 4.0);
        let b = vec3(5.0, 6.0, 7.0);
        let c = a * b;

        assert_eq!(c[0], 10.0);
        assert_eq!(c[1], 18.0);
        assert_eq!(c[2], 28.0);
    }

    #[test]
    fn test_divide() {
        let a = vec3(10.0, 18.0, 28.0);
        let b = vec3(2.0, 3.0, 4.0);
        let c = divide(&a, &b);

        assert_eq!(c[0], 5.0);
        assert_eq!(c[1], 6.0);
        assert_eq!(c[2], 7.0);
    }

    #[test]
    fn test_divide_operator() {
        let a = vec3(10.0, 18.0, 28.0);
        let b = vec3(2.0, 3.0, 4.0);
        let c = a / b;

        assert_eq!(c[0], 5.0);
        assert_eq!(c[1], 6.0);
        assert_eq!(c[2], 7.0);
    }

    #[test]
    fn test_negate() {
        let a = vec3(1.0, -2.0, 3.0);
        let b = negate(&a);

        assert_eq!(b[0], -1.0);
        assert_eq!(b[1], 2.0);
        assert_eq!(b[2], -3.0);
    }

    #[test]
    fn test_negate_operator() {
        let a = vec3(1.0, -2.0, 3.0);
        let b = -a;

        assert_eq!(b[0], -1.0);
        assert_eq!(b[1], 2.0);
        assert_eq!(b[2], -3.0);
    }

    #[test]
    fn test_reciprocal() {
        let a = vec3(2.0, 4.0, 5.0);
        let b = reciprocal(&a);

        assert_eq!(b[0], 0.5);
        assert_eq!(b[1], 0.25);
        assert_eq!(b[2], 0.2);
    }

    #[test]
    fn test_addition_properties() {
        let a = vec3(1.0, 2.0, 3.0);
        let b = vec3(4.0, 5.0, 6.0);
        let c = vec3(7.0, 8.0, 9.0);

        // Commutativity: a + b = b + a
        let ab = &a + &b;
        let ba = &b + &a;
        assert_eq!(ab[0], ba[0]);
        assert_eq!(ab[1], ba[1]);
        assert_eq!(ab[2], ba[2]);

        // Associativity: (a + b) + c = a + (b + c)
        let abc1 = (&a + &b) + &c;
        let abc2 = &a + (&b + &c);
        assert_eq!(abc1[0], abc2[0]);
        assert_eq!(abc1[1], abc2[1]);
        assert_eq!(abc1[2], abc2[2]);

        // Identity: a + 0 = a
        let zero = Vector::zero();
        let a_plus_zero = &a + &zero;
        assert_eq!(a_plus_zero[0], a[0]);
        assert_eq!(a_plus_zero[1], a[1]);
        assert_eq!(a_plus_zero[2], a[2]);

        // Inverse: a + (-a) = 0
        let neg_a = -&a;
        let result = &a + &neg_a;
        assert_eq!(result[0], 0.0);
        assert_eq!(result[1], 0.0);
        assert_eq!(result[2], 0.0);
    }

    #[test]
    fn test_multiplication_properties() {
        let a = vec3(2.0, 3.0, 4.0);
        let b = vec3(5.0, 6.0, 7.0);

        // Commutativity: a * b = b * a
        let ab = &a * &b;
        let ba = &b * &a;
        assert_eq!(ab[0], ba[0]);
        assert_eq!(ab[1], ba[1]);
        assert_eq!(ab[2], ba[2]);
    }

    #[test]
    fn test_zero_vector() {
        let zero: Vector<f64, 3> = Vector::zero();
        let a = vec3(1.0, 2.0, 3.0);

        // Addition with zero
        let result = &a + &zero;
        assert_eq!(result[0], a[0]);
        assert_eq!(result[1], a[1]);
        assert_eq!(result[2], a[2]);

        // Multiplication by zero
        let result = &a * &zero;
        assert_eq!(result[0], 0.0);
        assert_eq!(result[1], 0.0);
        assert_eq!(result[2], 0.0);
    }
}
