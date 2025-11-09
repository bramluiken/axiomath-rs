//! Element-wise arithmetic operations on matrices.
//!
//! This module implements component-wise (Hadamard) operations for matrices:
//! - Addition and subtraction
//! - Element-wise multiplication (Hadamard product)
//! - Element-wise division
//! - Scalar multiplication and division
//! - Negation
//!
//! # Note on Matrix Multiplication
//!
//! Standard matrix multiplication is NOT element-wise and is implemented
//! in the `multiplication` module.
//!
//! # Examples
//!
//! ```
//! use axiomath::matrix::Matrix;
//!
//! let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
//! let b = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);
//!
//! // Element-wise addition
//! let sum = a + b;
//! assert_eq!(sum[(0, 0)], 6.0);
//!
//! // Scalar multiplication
//! let scaled = a * 2.0;
//! assert_eq!(scaled[(1, 1)], 8.0);
//! ```

use crate::core::Scalar;
use crate::matrix::Matrix;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Adds two matrices element-wise.
///
/// # Mathematical Definition
///
/// ```text
/// (A + B)_ij = A_ij + B_ij
/// ```
///
/// # Properties
///
/// - Commutative: A + B = B + A
/// - Associative: (A + B) + C = A + (B + C)
/// - Identity: A + O = A (where O is the zero matrix)
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, add};
///
/// let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
/// let b = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);
/// let sum = add(&a, &b);
///
/// assert_eq!(sum[(0, 0)], 6.0);
/// assert_eq!(sum[(1, 1)], 12.0);
/// ```
pub fn add<T: Scalar, const ROWS: usize, const COLS: usize>(
    a: &Matrix<T, ROWS, COLS>,
    b: &Matrix<T, ROWS, COLS>,
) -> Matrix<T, ROWS, COLS> {
    a.zip_with(*b, |x, y| x + y)
}

/// Subtracts one matrix from another element-wise.
///
/// # Mathematical Definition
///
/// ```text
/// (A - B)_ij = A_ij - B_ij
/// ```
///
/// # Properties
///
/// - Anti-commutative: A - B = -(B - A)
/// - Identity: A - O = A
/// - Self-cancellation: A - A = O
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, subtract};
///
/// let a = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);
/// let b = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
/// let diff = subtract(&a, &b);
///
/// assert_eq!(diff[(0, 0)], 4.0);
/// assert_eq!(diff[(1, 1)], 4.0);
/// ```
pub fn subtract<T: Scalar, const ROWS: usize, const COLS: usize>(
    a: &Matrix<T, ROWS, COLS>,
    b: &Matrix<T, ROWS, COLS>,
) -> Matrix<T, ROWS, COLS> {
    a.zip_with(*b, |x, y| x - y)
}

/// Multiplies two matrices element-wise (Hadamard product).
///
/// # Mathematical Definition
///
/// ```text
/// (A ∘ B)_ij = A_ij · B_ij
/// ```
///
/// This is NOT standard matrix multiplication.
///
/// # Properties
///
/// - Commutative: A ∘ B = B ∘ A
/// - Associative: (A ∘ B) ∘ C = A ∘ (B ∘ C)
/// - Distributive: A ∘ (B + C) = A ∘ B + A ∘ C
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, multiply_elementwise};
///
/// let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
/// let b = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);
/// let hadamard = multiply_elementwise(&a, &b);
///
/// assert_eq!(hadamard[(0, 0)], 5.0);
/// assert_eq!(hadamard[(0, 1)], 12.0);
/// assert_eq!(hadamard[(1, 0)], 21.0);
/// assert_eq!(hadamard[(1, 1)], 32.0);
/// ```
pub fn multiply_elementwise<T: Scalar, const ROWS: usize, const COLS: usize>(
    a: &Matrix<T, ROWS, COLS>,
    b: &Matrix<T, ROWS, COLS>,
) -> Matrix<T, ROWS, COLS> {
    a.zip_with(*b, |x, y| x * y)
}

/// Divides one matrix by another element-wise.
///
/// # Mathematical Definition
///
/// ```text
/// (A / B)_ij = A_ij / B_ij
/// ```
///
/// # Warning
///
/// This will panic or produce infinity/NaN if any element of B is zero.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, divide_elementwise};
///
/// let a = Matrix::<f64, 2, 2>::from_rows([[10.0, 20.0], [30.0, 40.0]]);
/// let b = Matrix::<f64, 2, 2>::from_rows([[2.0, 4.0], [5.0, 8.0]]);
/// let quotient = divide_elementwise(&a, &b);
///
/// assert_eq!(quotient[(0, 0)], 5.0);
/// assert_eq!(quotient[(0, 1)], 5.0);
/// assert_eq!(quotient[(1, 0)], 6.0);
/// assert_eq!(quotient[(1, 1)], 5.0);
/// ```
pub fn divide_elementwise<T: Scalar, const ROWS: usize, const COLS: usize>(
    a: &Matrix<T, ROWS, COLS>,
    b: &Matrix<T, ROWS, COLS>,
) -> Matrix<T, ROWS, COLS> {
    a.zip_with(*b, |x, y| x / y)
}

/// Multiplies a matrix by a scalar.
///
/// # Mathematical Definition
///
/// ```text
/// (λA)_ij = λ · A_ij
/// ```
///
/// # Properties
///
/// - Distributive: λ(A + B) = λA + λB
/// - Associative: (λμ)A = λ(μA)
/// - Identity: 1·A = A
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, scale};
///
/// let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
/// let scaled = scale(&a, 3.0);
///
/// assert_eq!(scaled[(0, 0)], 3.0);
/// assert_eq!(scaled[(1, 1)], 12.0);
/// ```
pub fn scale<T: Scalar, const ROWS: usize, const COLS: usize>(
    matrix: &Matrix<T, ROWS, COLS>,
    scalar: T,
) -> Matrix<T, ROWS, COLS> {
    matrix.map(|x| x * scalar)
}

/// Divides a matrix by a scalar.
///
/// # Mathematical Definition
///
/// ```text
/// (A/λ)_ij = A_ij / λ
/// ```
///
/// # Warning
///
/// This will panic or produce infinity/NaN if the scalar is zero.
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, divide_by_scalar};
///
/// let a = Matrix::<f64, 2, 2>::from_rows([[10.0, 20.0], [30.0, 40.0]]);
/// let divided = divide_by_scalar(&a, 10.0);
///
/// assert_eq!(divided[(0, 0)], 1.0);
/// assert_eq!(divided[(1, 1)], 4.0);
/// ```
pub fn divide_by_scalar<T: Scalar, const ROWS: usize, const COLS: usize>(
    matrix: &Matrix<T, ROWS, COLS>,
    scalar: T,
) -> Matrix<T, ROWS, COLS> {
    matrix.map(|x| x / scalar)
}

/// Negates all elements of a matrix.
///
/// # Mathematical Definition
///
/// ```text
/// (-A)_ij = -A_ij
/// ```
///
/// # Properties
///
/// - Double negation: -(-A) = A
/// - Relation to subtraction: -A = O - A
///
/// # Examples
///
/// ```
/// use axiomath::matrix::{Matrix, negate};
///
/// let a = Matrix::<f64, 2, 2>::from_rows([[1.0, -2.0], [3.0, -4.0]]);
/// let neg = negate(&a);
///
/// assert_eq!(neg[(0, 0)], -1.0);
/// assert_eq!(neg[(0, 1)], 2.0);
/// assert_eq!(neg[(1, 0)], -3.0);
/// assert_eq!(neg[(1, 1)], 4.0);
/// ```
pub fn negate<T: Scalar, const ROWS: usize, const COLS: usize>(
    matrix: &Matrix<T, ROWS, COLS>,
) -> Matrix<T, ROWS, COLS> {
    matrix.map(|x| T::zero() - x)
}

// Operator overloading: Matrix + Matrix
impl<T: Scalar, const ROWS: usize, const COLS: usize> Add for Matrix<T, ROWS, COLS> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        add(&self, &other)
    }
}

impl<T: Scalar, const ROWS: usize, const COLS: usize> Add<&Matrix<T, ROWS, COLS>> for Matrix<T, ROWS, COLS> {
    type Output = Matrix<T, ROWS, COLS>;

    fn add(self, other: &Matrix<T, ROWS, COLS>) -> Self::Output {
        add(&self, other)
    }
}

impl<T: Scalar, const ROWS: usize, const COLS: usize> Add<Matrix<T, ROWS, COLS>> for &Matrix<T, ROWS, COLS> {
    type Output = Matrix<T, ROWS, COLS>;

    fn add(self, other: Matrix<T, ROWS, COLS>) -> Self::Output {
        add(self, &other)
    }
}

impl<T: Scalar, const ROWS: usize, const COLS: usize> Add for &Matrix<T, ROWS, COLS> {
    type Output = Matrix<T, ROWS, COLS>;

    fn add(self, other: Self) -> Self::Output {
        add(self, other)
    }
}

// Operator overloading: Matrix - Matrix
impl<T: Scalar, const ROWS: usize, const COLS: usize> Sub for Matrix<T, ROWS, COLS> {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        subtract(&self, &other)
    }
}

impl<T: Scalar, const ROWS: usize, const COLS: usize> Sub<&Matrix<T, ROWS, COLS>> for Matrix<T, ROWS, COLS> {
    type Output = Matrix<T, ROWS, COLS>;

    fn sub(self, other: &Matrix<T, ROWS, COLS>) -> Self::Output {
        subtract(&self, other)
    }
}

impl<T: Scalar, const ROWS: usize, const COLS: usize> Sub<Matrix<T, ROWS, COLS>> for &Matrix<T, ROWS, COLS> {
    type Output = Matrix<T, ROWS, COLS>;

    fn sub(self, other: Matrix<T, ROWS, COLS>) -> Self::Output {
        subtract(self, &other)
    }
}

impl<T: Scalar, const ROWS: usize, const COLS: usize> Sub for &Matrix<T, ROWS, COLS> {
    type Output = Matrix<T, ROWS, COLS>;

    fn sub(self, other: Self) -> Self::Output {
        subtract(self, other)
    }
}

// Operator overloading: -Matrix
impl<T: Scalar, const ROWS: usize, const COLS: usize> Neg for Matrix<T, ROWS, COLS> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        negate(&self)
    }
}

impl<T: Scalar, const ROWS: usize, const COLS: usize> Neg for &Matrix<T, ROWS, COLS> {
    type Output = Matrix<T, ROWS, COLS>;

    fn neg(self) -> Self::Output {
        negate(self)
    }
}

// Operator overloading: Matrix * scalar
impl<T: Scalar, const ROWS: usize, const COLS: usize> Mul<T> for Matrix<T, ROWS, COLS> {
    type Output = Self;

    fn mul(self, scalar: T) -> Self::Output {
        scale(&self, scalar)
    }
}

impl<T: Scalar, const ROWS: usize, const COLS: usize> Mul<T> for &Matrix<T, ROWS, COLS> {
    type Output = Matrix<T, ROWS, COLS>;

    fn mul(self, scalar: T) -> Self::Output {
        scale(self, scalar)
    }
}

// Operator overloading: Matrix / scalar
impl<T: Scalar, const ROWS: usize, const COLS: usize> Div<T> for Matrix<T, ROWS, COLS> {
    type Output = Self;

    fn div(self, scalar: T) -> Self::Output {
        divide_by_scalar(&self, scalar)
    }
}

impl<T: Scalar, const ROWS: usize, const COLS: usize> Div<T> for &Matrix<T, ROWS, COLS> {
    type Output = Matrix<T, ROWS, COLS>;

    fn div(self, scalar: T) -> Self::Output {
        divide_by_scalar(self, scalar)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add() {
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let b = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);
        let sum = add(&a, &b);

        assert_eq!(sum[(0, 0)], 6.0);
        assert_eq!(sum[(0, 1)], 8.0);
        assert_eq!(sum[(1, 0)], 10.0);
        assert_eq!(sum[(1, 1)], 12.0);
    }

    #[test]
    fn test_add_operator() {
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let b = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);

        let sum1 = a + b;
        let sum2 = &a + &b;
        let sum3 = a + &b;
        let sum4 = &a + b;

        assert_eq!(sum1, sum2);
        assert_eq!(sum1, sum3);
        assert_eq!(sum1, sum4);
    }

    #[test]
    fn test_add_commutative() {
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let b = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);

        assert_eq!(add(&a, &b), add(&b, &a));
    }

    #[test]
    fn test_add_associative() {
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let b = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);
        let c = Matrix::<f64, 2, 2>::from_rows([[9.0, 10.0], [11.0, 12.0]]);

        let left = add(&add(&a, &b), &c);
        let right = add(&a, &add(&b, &c));
        assert_eq!(left, right);
    }

    #[test]
    fn test_add_identity() {
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let zero = Matrix::<f64, 2, 2>::zero();

        assert_eq!(add(&a, &zero), a);
    }

    #[test]
    fn test_subtract() {
        let a = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);
        let b = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let diff = subtract(&a, &b);

        assert_eq!(diff[(0, 0)], 4.0);
        assert_eq!(diff[(0, 1)], 4.0);
        assert_eq!(diff[(1, 0)], 4.0);
        assert_eq!(diff[(1, 1)], 4.0);
    }

    #[test]
    fn test_subtract_operator() {
        let a = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);
        let b = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);

        let diff1 = a - b;
        let diff2 = &a - &b;

        assert_eq!(diff1, diff2);
    }

    #[test]
    fn test_subtract_self() {
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let result = subtract(&a, &a);

        assert!(result.is_zero());
    }

    #[test]
    fn test_multiply_elementwise() {
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let b = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);
        let product = multiply_elementwise(&a, &b);

        assert_eq!(product[(0, 0)], 5.0);
        assert_eq!(product[(0, 1)], 12.0);
        assert_eq!(product[(1, 0)], 21.0);
        assert_eq!(product[(1, 1)], 32.0);
    }

    #[test]
    fn test_multiply_elementwise_commutative() {
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let b = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);

        assert_eq!(multiply_elementwise(&a, &b), multiply_elementwise(&b, &a));
    }

    #[test]
    fn test_divide_elementwise() {
        let a = Matrix::<f64, 2, 2>::from_rows([[10.0, 20.0], [30.0, 40.0]]);
        let b = Matrix::<f64, 2, 2>::from_rows([[2.0, 4.0], [5.0, 8.0]]);
        let quotient = divide_elementwise(&a, &b);

        assert_eq!(quotient[(0, 0)], 5.0);
        assert_eq!(quotient[(0, 1)], 5.0);
        assert_eq!(quotient[(1, 0)], 6.0);
        assert_eq!(quotient[(1, 1)], 5.0);
    }

    #[test]
    fn test_scale() {
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let scaled = scale(&a, 3.0);

        assert_eq!(scaled[(0, 0)], 3.0);
        assert_eq!(scaled[(0, 1)], 6.0);
        assert_eq!(scaled[(1, 0)], 9.0);
        assert_eq!(scaled[(1, 1)], 12.0);
    }

    #[test]
    fn test_scale_operator() {
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let scaled1 = a * 3.0;
        let scaled2 = &a * 3.0;

        assert_eq!(scaled1, scaled2);
    }

    #[test]
    fn test_scale_identity() {
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let scaled = scale(&a, 1.0);

        assert_eq!(scaled, a);
    }

    #[test]
    fn test_scale_distributive() {
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let b = Matrix::<f64, 2, 2>::from_rows([[5.0, 6.0], [7.0, 8.0]]);
        let scalar = 2.0;

        let left = scale(&add(&a, &b), scalar);
        let right = add(&scale(&a, scalar), &scale(&b, scalar));

        assert_eq!(left, right);
    }

    #[test]
    fn test_divide_by_scalar() {
        let a = Matrix::<f64, 2, 2>::from_rows([[10.0, 20.0], [30.0, 40.0]]);
        let divided = divide_by_scalar(&a, 10.0);

        assert_eq!(divided[(0, 0)], 1.0);
        assert_eq!(divided[(0, 1)], 2.0);
        assert_eq!(divided[(1, 0)], 3.0);
        assert_eq!(divided[(1, 1)], 4.0);
    }

    #[test]
    fn test_divide_by_scalar_operator() {
        let a = Matrix::<f64, 2, 2>::from_rows([[10.0, 20.0], [30.0, 40.0]]);
        let divided1 = a / 10.0;
        let divided2 = &a / 10.0;

        assert_eq!(divided1, divided2);
    }

    #[test]
    fn test_negate() {
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, -2.0], [3.0, -4.0]]);
        let neg = negate(&a);

        assert_eq!(neg[(0, 0)], -1.0);
        assert_eq!(neg[(0, 1)], 2.0);
        assert_eq!(neg[(1, 0)], -3.0);
        assert_eq!(neg[(1, 1)], 4.0);
    }

    #[test]
    fn test_negate_operator() {
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, -2.0], [3.0, -4.0]]);
        let neg1 = -a;
        let neg2 = -&a;

        assert_eq!(neg1, neg2);
    }

    #[test]
    fn test_double_negation() {
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let double_neg = negate(&negate(&a));

        assert_eq!(double_neg, a);
    }

    #[test]
    fn test_negate_and_subtract_relation() {
        let a = Matrix::<f64, 2, 2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let zero = Matrix::<f64, 2, 2>::zero();

        assert_eq!(negate(&a), subtract(&zero, &a));
    }
}
