//! Core vector type definition and basic operations.
//!
//! This module defines the fundamental `Vector<T, N>` type using const generics
//! for compile-time dimension safety.

use crate::core::traits::Scalar;
use core::fmt;
use core::ops::{Index, IndexMut};

/// A generic N-dimensional vector.
///
/// `Vector` is a fixed-size array-backed vector type that uses const generics
/// to ensure compile-time dimension checking. It is generic over the scalar type `T`
/// and the dimension `N`.
///
/// # Type Parameters
///
/// * `T` - The scalar type (must implement `Scalar` trait)
/// * `N` - The dimension of the vector (compile-time constant)
///
/// # Examples
///
/// ```
/// use axiomath::vector::Vector;
///
/// // Create a 3D vector
/// let v = Vector::new([1.0, 2.0, 3.0]);
/// assert_eq!(v[0], 1.0);
/// assert_eq!(v[1], 2.0);
/// assert_eq!(v[2], 3.0);
///
/// // Create a zero vector
/// let zero: Vector<f64, 3> = Vector::zero();
/// assert_eq!(zero[0], 0.0);
/// ```
#[derive(Copy, Clone, PartialEq, Eq, Hash)]
pub struct Vector<T, const N: usize> {
    pub(crate) data: [T; N],
}

impl<T, const N: usize> Vector<T, N>
where
    T: Scalar,
{
    /// Creates a new vector from an array.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::vector::Vector;
    ///
    /// let v = Vector::new([1.0, 2.0, 3.0]);
    /// assert_eq!(v[0], 1.0);
    /// ```
    #[inline]
    pub const fn new(data: [T; N]) -> Self {
        Self { data }
    }

    /// Creates a zero vector (all components are zero).
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::vector::Vector;
    ///
    /// let v: Vector<f64, 3> = Vector::zero();
    /// assert_eq!(v[0], 0.0);
    /// assert_eq!(v[1], 0.0);
    /// assert_eq!(v[2], 0.0);
    /// ```
    #[inline]
    pub fn zero() -> Self {
        Self {
            data: [T::zero(); N],
        }
    }

    /// Creates a vector with all components set to one.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::vector::Vector;
    ///
    /// let v: Vector<f64, 3> = Vector::one();
    /// assert_eq!(v[0], 1.0);
    /// assert_eq!(v[1], 1.0);
    /// assert_eq!(v[2], 1.0);
    /// ```
    #[inline]
    pub fn one() -> Self {
        Self {
            data: [T::one(); N],
        }
    }

    /// Creates a vector from a function that maps each index to a value.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::vector::Vector;
    ///
    /// let v: Vector<f64, 5> = Vector::from_fn(|i| (i * 2) as f64);
    /// assert_eq!(v[0], 0.0);
    /// assert_eq!(v[1], 2.0);
    /// assert_eq!(v[2], 4.0);
    /// ```
    #[inline]
    pub fn from_fn<F>(mut f: F) -> Self
    where
        F: FnMut(usize) -> T,
    {
        let mut data = [T::zero(); N];
        for i in 0..N {
            data[i] = f(i);
        }
        Self { data }
    }

    /// Returns the dimension of the vector.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::vector::Vector;
    ///
    /// let v = Vector::new([1.0, 2.0, 3.0]);
    /// assert_eq!(v.dimension(), 3);
    /// ```
    #[inline]
    pub const fn dimension(&self) -> usize {
        N
    }

    /// Returns a reference to the underlying array.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::vector::Vector;
    ///
    /// let v = Vector::new([1.0, 2.0, 3.0]);
    /// let array = v.as_array();
    /// assert_eq!(array, &[1.0, 2.0, 3.0]);
    /// ```
    #[inline]
    pub const fn as_array(&self) -> &[T; N] {
        &self.data
    }

    /// Returns a mutable reference to the underlying array.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::vector::Vector;
    ///
    /// let mut v = Vector::new([1.0, 2.0, 3.0]);
    /// v.as_array_mut()[0] = 5.0;
    /// assert_eq!(v[0], 5.0);
    /// ```
    #[inline]
    pub fn as_array_mut(&mut self) -> &mut [T; N] {
        &mut self.data
    }

    /// Returns an iterator over the vector components.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::vector::Vector;
    ///
    /// let v = Vector::new([1.0, 2.0, 3.0]);
    /// let sum: f64 = v.iter().sum();
    /// assert_eq!(sum, 6.0);
    /// ```
    #[inline]
    pub fn iter(&self) -> core::slice::Iter<'_, T> {
        self.data.iter()
    }

    /// Returns a mutable iterator over the vector components.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::vector::Vector;
    ///
    /// let mut v = Vector::new([1.0, 2.0, 3.0]);
    /// for component in v.iter_mut() {
    ///     *component *= 2.0;
    /// }
    /// assert_eq!(v[0], 2.0);
    /// assert_eq!(v[1], 4.0);
    /// ```
    #[inline]
    pub fn iter_mut(&mut self) -> core::slice::IterMut<'_, T> {
        self.data.iter_mut()
    }

    /// Checks if all components are zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::vector::Vector;
    ///
    /// let zero: Vector<f64, 3> = Vector::zero();
    /// assert!(zero.is_zero());
    ///
    /// let v = Vector::new([1.0, 0.0, 0.0]);
    /// assert!(!v.is_zero());
    /// ```
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.data.iter().all(|&x| x.is_zero())
    }

    /// Applies a function to each component of the vector.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::vector::Vector;
    ///
    /// let v = Vector::new([1.0, 2.0, 3.0]);
    /// let doubled = v.map(|x| x * 2.0);
    /// assert_eq!(doubled[0], 2.0);
    /// assert_eq!(doubled[1], 4.0);
    /// assert_eq!(doubled[2], 6.0);
    /// ```
    #[inline]
    pub fn map<F>(&self, mut f: F) -> Self
    where
        F: FnMut(T) -> T,
    {
        let mut result = *self;
        for i in 0..N {
            result.data[i] = f(self.data[i]);
        }
        result
    }

    /// Applies a function to corresponding components of two vectors.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::vector::Vector;
    ///
    /// let a = Vector::new([1.0, 2.0, 3.0]);
    /// let b = Vector::new([4.0, 5.0, 6.0]);
    /// let sum = a.zip_with(&b, |x, y| x + y);
    /// assert_eq!(sum[0], 5.0);
    /// assert_eq!(sum[1], 7.0);
    /// assert_eq!(sum[2], 9.0);
    /// ```
    #[inline]
    pub fn zip_with<F>(&self, other: &Self, mut f: F) -> Self
    where
        F: FnMut(T, T) -> T,
    {
        let mut result = *self;
        for i in 0..N {
            result.data[i] = f(self.data[i], other.data[i]);
        }
        result
    }
}

// Index trait implementation for element access
impl<T, const N: usize> Index<usize> for Vector<T, N>
where
    T: Scalar,
{
    type Output = T;

    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

// IndexMut trait implementation for mutable element access
impl<T, const N: usize> IndexMut<usize> for Vector<T, N>
where
    T: Scalar,
{
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}

// From trait for conversion from array
impl<T, const N: usize> From<[T; N]> for Vector<T, N>
where
    T: Scalar,
{
    #[inline]
    fn from(data: [T; N]) -> Self {
        Self::new(data)
    }
}

// Into trait for conversion to array
impl<T, const N: usize> From<Vector<T, N>> for [T; N]
where
    T: Scalar,
{
    #[inline]
    fn from(vector: Vector<T, N>) -> Self {
        vector.data
    }
}

// Debug trait implementation
impl<T, const N: usize> fmt::Debug for Vector<T, N>
where
    T: Scalar + fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Vector")
            .field("data", &self.data)
            .finish()
    }
}

// Display trait implementation
impl<T, const N: usize> fmt::Display for Vector<T, N>
where
    T: Scalar + fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[")?;
        for (i, component) in self.data.iter().enumerate() {
            if i > 0 {
                write!(f, ", ")?;
            }
            write!(f, "{}", component)?;
        }
        write!(f, "]")
    }
}

// Default trait implementation
impl<T, const N: usize> Default for Vector<T, N>
where
    T: Scalar,
{
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let v = Vector::new([1.0, 2.0, 3.0]);
        assert_eq!(v[0], 1.0);
        assert_eq!(v[1], 2.0);
        assert_eq!(v[2], 3.0);
    }

    #[test]
    fn test_zero() {
        let v: Vector<f64, 3> = Vector::zero();
        assert_eq!(v[0], 0.0);
        assert_eq!(v[1], 0.0);
        assert_eq!(v[2], 0.0);
        assert!(v.is_zero());
    }

    #[test]
    fn test_one() {
        let v: Vector<f64, 3> = Vector::one();
        assert_eq!(v[0], 1.0);
        assert_eq!(v[1], 1.0);
        assert_eq!(v[2], 1.0);
    }

    #[test]
    fn test_from_fn() {
        let v: Vector<f64, 5> = Vector::from_fn(|i| (i * 2) as f64);
        assert_eq!(v[0], 0.0);
        assert_eq!(v[1], 2.0);
        assert_eq!(v[2], 4.0);
        assert_eq!(v[3], 6.0);
        assert_eq!(v[4], 8.0);
    }

    #[test]
    fn test_dimension() {
        let v = Vector::new([1.0, 2.0, 3.0]);
        assert_eq!(v.dimension(), 3);

        let v2: Vector<f64, 5> = Vector::zero();
        assert_eq!(v2.dimension(), 5);
    }

    #[test]
    fn test_indexing() {
        let mut v = Vector::new([1.0, 2.0, 3.0]);
        assert_eq!(v[0], 1.0);
        assert_eq!(v[1], 2.0);
        assert_eq!(v[2], 3.0);

        v[0] = 10.0;
        assert_eq!(v[0], 10.0);
    }

    #[test]
    fn test_is_zero() {
        let zero: Vector<f64, 3> = Vector::zero();
        assert!(zero.is_zero());

        let v = Vector::new([1.0, 0.0, 0.0]);
        assert!(!v.is_zero());

        let v2 = Vector::new([0.0, 0.0, 0.0]);
        assert!(v2.is_zero());
    }

    #[test]
    fn test_map() {
        let v = Vector::new([1.0, 2.0, 3.0]);
        let doubled = v.map(|x| x * 2.0);
        assert_eq!(doubled[0], 2.0);
        assert_eq!(doubled[1], 4.0);
        assert_eq!(doubled[2], 6.0);
    }

    #[test]
    fn test_zip_with() {
        let a = Vector::new([1.0, 2.0, 3.0]);
        let b = Vector::new([4.0, 5.0, 6.0]);
        let sum = a.zip_with(&b, |x, y| x + y);
        assert_eq!(sum[0], 5.0);
        assert_eq!(sum[1], 7.0);
        assert_eq!(sum[2], 9.0);
    }

    #[test]
    fn test_from_array() {
        let v = Vector::from([1.0, 2.0, 3.0]);
        assert_eq!(v[0], 1.0);
        assert_eq!(v[1], 2.0);
        assert_eq!(v[2], 3.0);
    }

    #[test]
    fn test_into_array() {
        let v = Vector::new([1.0, 2.0, 3.0]);
        let array: [f64; 3] = v.into();
        assert_eq!(array, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_iter() {
        let v = Vector::new([1.0, 2.0, 3.0]);
        let sum: f64 = v.iter().sum();
        assert_eq!(sum, 6.0);
    }

    #[test]
    fn test_iter_mut() {
        let mut v = Vector::new([1.0, 2.0, 3.0]);
        for component in v.iter_mut() {
            *component *= 2.0;
        }
        assert_eq!(v[0], 2.0);
        assert_eq!(v[1], 4.0);
        assert_eq!(v[2], 6.0);
    }

    #[test]
    fn test_display() {
        let v = Vector::new([1.0, 2.0, 3.0]);
        let s = format!("{}", v);
        assert_eq!(s, "[1, 2, 3]");
    }

    #[test]
    fn test_copy_clone() {
        let v1 = Vector::new([1.0, 2.0, 3.0]);
        let v2 = v1; // Copy
        let v3 = v1.clone(); // Clone

        assert_eq!(v1[0], v2[0]);
        assert_eq!(v1[0], v3[0]);
    }

    #[test]
    fn test_different_dimensions() {
        let v2: Vector<f64, 2> = Vector::new([1.0, 2.0]);
        let v3: Vector<f64, 3> = Vector::new([1.0, 2.0, 3.0]);
        let v4: Vector<f64, 4> = Vector::new([1.0, 2.0, 3.0, 4.0]);

        assert_eq!(v2.dimension(), 2);
        assert_eq!(v3.dimension(), 3);
        assert_eq!(v4.dimension(), 4);
    }

    #[test]
    fn test_default() {
        let v: Vector<f64, 3> = Vector::default();
        assert!(v.is_zero());
    }
}
