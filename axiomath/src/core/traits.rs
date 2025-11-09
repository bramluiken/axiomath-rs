//! Core trait definitions for numerical types.
//!
//! This module defines the fundamental traits that all numerical types must implement
//! to work with axiomath. These traits build upon `num-traits` to provide a
//! comprehensive type system for mathematical operations.

use num_traits::{Float as NumFloat, Num, NumAssign, NumCast, Signed};
use std::fmt::{Debug, Display};
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Base trait for scalar numerical types.
///
/// This trait represents any numeric type that can be used as a scalar
/// in mathematical operations. It combines basic arithmetic, comparison,
/// and conversion capabilities.
///
/// # Mathematical Properties
///
/// A type implementing `Scalar` should satisfy:
/// - Additive identity: `a + 0 = a`
/// - Multiplicative identity: `a * 1 = a`
/// - Commutativity: `a + b = b + a`, `a * b = b * a`
/// - Associativity: `(a + b) + c = a + (b + c)`
/// - Distributivity: `a * (b + c) = a * b + a * c`
///
/// # Examples
///
/// ```
/// use axiomath::core::Scalar;
///
/// fn square<T: Scalar>(x: T) -> T {
///     x * x
/// }
///
/// assert_eq!(square(3.0_f64), 9.0);
/// assert_eq!(square(4.0_f32), 16.0);
/// ```
pub trait Scalar:
    Num
    + NumAssign
    + NumCast
    + PartialOrd
    + Copy
    + Clone
    + Debug
    + Display
    + Send
    + Sync
    + 'static
{
    /// Returns the absolute value of self.
    ///
    /// For signed types, this is |x|. For unsigned types, this is a no-op.
    fn abs(self) -> Self;

    /// Returns the sign of self.
    ///
    /// - Returns `1` if positive
    /// - Returns `-1` if negative
    /// - Returns `0` if zero
    fn signum(self) -> Self;

    /// Raises self to the power of `exp`.
    fn powi(self, exp: i32) -> Self;
}

/// Trait for real number types.
///
/// Real numbers include all rational and irrational numbers on the continuous
/// number line. This trait extends `Scalar` with operations specific to real numbers.
///
/// # Mathematical Properties
///
/// A type implementing `Real` should satisfy all `Scalar` properties plus:
/// - Total ordering: For any a, b either a < b, a = b, or a > b
/// - Completeness: Every bounded sequence has a least upper bound
/// - Archimedean property: For any x > 0, there exists n such that nx > 1
///
/// # Examples
///
/// ```
/// use axiomath::core::Real;
///
/// fn distance<T: Real>(a: T, b: T) -> T {
///     (a - b).abs()
/// }
///
/// assert_eq!(distance(5.0_f64, 2.0_f64), 3.0);
/// ```
pub trait Real: Scalar + Signed + PartialOrd {
    /// Returns the square root of self.
    ///
    /// # Mathematical Definition
    ///
    /// √x is defined as the unique non-negative number y such that y² = x.
    ///
    /// # Properties
    ///
    /// - √(x²) = |x|
    /// - √(xy) = √x * √y for x, y ≥ 0
    /// - Returns NaN for negative inputs (for float types)
    fn sqrt(self) -> Self;

    /// Returns the cube root of self.
    ///
    /// # Mathematical Definition
    ///
    /// ∛x is defined as the unique number y such that y³ = x.
    fn cbrt(self) -> Self;

    /// Computes the square of self.
    ///
    /// This is a convenience method that may be more efficient than `self * self`
    /// for some types.
    #[inline]
    fn square(self) -> Self {
        self * self
    }

    /// Computes the cube of self.
    #[inline]
    fn cube(self) -> Self {
        self * self * self
    }

    /// Returns the reciprocal (multiplicative inverse) of self.
    ///
    /// # Mathematical Definition
    ///
    /// 1/x is the unique number y such that xy = 1.
    ///
    /// # Panics
    ///
    /// May panic or return infinity if self is zero.
    #[inline]
    fn recip(self) -> Self {
        Self::one() / self
    }
}

/// Trait for floating-point types.
///
/// This trait represents IEEE 754 floating-point numbers and provides
/// access to transcendental functions, special values, and floating-point
/// specific operations.
///
/// # Mathematical Functions
///
/// This trait provides:
/// - Exponential and logarithmic functions
/// - Trigonometric functions
/// - Hyperbolic functions
/// - Special float values (NaN, infinity)
///
/// # Examples
///
/// ```
/// use axiomath::core::Float;
///
/// let x = 1.0_f64;
/// assert!(x.is_finite());
/// assert!(!x.is_nan());
/// ```
pub trait Float: Real {
    /// Returns the machine epsilon for this float type.
    ///
    /// This is the difference between 1.0 and the next larger representable number.
    fn epsilon() -> Self;

    /// Returns positive infinity.
    fn infinity() -> Self;

    /// Returns negative infinity.
    fn neg_infinity() -> Self;

    /// Returns NaN (Not a Number).
    fn nan() -> Self;

    /// Returns true if self is NaN.
    fn is_nan(self) -> bool;

    /// Returns true if self is infinite (positive or negative).
    fn is_infinite(self) -> bool;

    /// Returns true if self is finite (not NaN, not infinite).
    fn is_finite(self) -> bool;

    /// Returns true if self is neither zero, infinite, subnormal, or NaN.
    fn is_normal(self) -> bool;

    /// Computes e^x (exponential function).
    ///
    /// # Mathematical Definition
    ///
    /// exp(x) = e^x = Σ(n=0 to ∞) x^n / n!
    ///
    /// # Properties
    ///
    /// - exp(0) = 1
    /// - exp(ln(x)) = x
    /// - exp(x + y) = exp(x) * exp(y)
    fn exp(self) -> Self;

    /// Computes 2^x.
    fn exp2(self) -> Self;

    /// Computes the natural logarithm (base e).
    ///
    /// # Mathematical Definition
    ///
    /// ln(x) is the inverse of exp(x), defined for x > 0.
    ///
    /// # Properties
    ///
    /// - ln(1) = 0
    /// - ln(e) = 1
    /// - ln(xy) = ln(x) + ln(y)
    /// - ln(x^n) = n * ln(x)
    fn ln(self) -> Self;

    /// Computes the logarithm base 2.
    fn log2(self) -> Self;

    /// Computes the logarithm base 10.
    fn log10(self) -> Self;

    /// Computes the logarithm with an arbitrary base.
    fn log(self, base: Self) -> Self;

    /// Computes the sine function.
    ///
    /// # Mathematical Definition
    ///
    /// sin(x) = Σ(n=0 to ∞) (-1)^n * x^(2n+1) / (2n+1)!
    ///
    /// # Properties
    ///
    /// - sin(0) = 0
    /// - sin(π/2) = 1
    /// - sin(-x) = -sin(x) (odd function)
    /// - Period: 2π
    fn sin(self) -> Self;

    /// Computes the cosine function.
    ///
    /// # Mathematical Definition
    ///
    /// cos(x) = Σ(n=0 to ∞) (-1)^n * x^(2n) / (2n)!
    ///
    /// # Properties
    ///
    /// - cos(0) = 1
    /// - cos(π) = -1
    /// - cos(-x) = cos(x) (even function)
    /// - Period: 2π
    fn cos(self) -> Self;

    /// Computes the tangent function.
    ///
    /// # Mathematical Definition
    ///
    /// tan(x) = sin(x) / cos(x)
    ///
    /// # Properties
    ///
    /// - tan(0) = 0
    /// - tan(π/4) = 1
    /// - Period: π
    /// - Undefined at x = π/2 + nπ
    fn tan(self) -> Self;

    /// Computes the arcsine (inverse sine).
    ///
    /// Returns a value in [-π/2, π/2].
    fn asin(self) -> Self;

    /// Computes the arccosine (inverse cosine).
    ///
    /// Returns a value in [0, π].
    fn acos(self) -> Self;

    /// Computes the arctangent (inverse tangent).
    ///
    /// Returns a value in (-π/2, π/2).
    fn atan(self) -> Self;

    /// Computes the four-quadrant arctangent of self (y) and other (x).
    ///
    /// Returns a value in [-π, π].
    ///
    /// # Mathematical Definition
    ///
    /// atan2(y, x) gives the angle θ in polar coordinates corresponding
    /// to the Cartesian coordinates (x, y).
    fn atan2(self, other: Self) -> Self;

    /// Computes the hyperbolic sine.
    ///
    /// # Mathematical Definition
    ///
    /// sinh(x) = (e^x - e^(-x)) / 2
    fn sinh(self) -> Self;

    /// Computes the hyperbolic cosine.
    ///
    /// # Mathematical Definition
    ///
    /// cosh(x) = (e^x + e^(-x)) / 2
    fn cosh(self) -> Self;

    /// Computes the hyperbolic tangent.
    ///
    /// # Mathematical Definition
    ///
    /// tanh(x) = sinh(x) / cosh(x) = (e^x - e^(-x)) / (e^x + e^(-x))
    fn tanh(self) -> Self;

    /// Raises self to the power of `exp`.
    ///
    /// # Mathematical Definition
    ///
    /// x^y = e^(y * ln(x)) for x > 0
    fn powf(self, exp: Self) -> Self;

    /// Returns the largest integer less than or equal to self.
    fn floor(self) -> Self;

    /// Returns the smallest integer greater than or equal to self.
    fn ceil(self) -> Self;

    /// Returns the nearest integer to self. Rounds half-way cases away from 0.0.
    fn round(self) -> Self;

    /// Returns the integer part of self (truncates towards zero).
    fn trunc(self) -> Self;

    /// Returns the fractional part of self.
    fn fract(self) -> Self;

    /// Computes self * a + b (fused multiply-add).
    ///
    /// This may be more accurate and faster than separate operations.
    fn mul_add(self, a: Self, b: Self) -> Self;

    /// Returns the minimum of two values.
    fn min(self, other: Self) -> Self;

    /// Returns the maximum of two values.
    fn max(self, other: Self) -> Self;

    /// Restricts self to be within the range [min, max].
    fn clamp(self, min: Self, max: Self) -> Self;
}

// Implement Scalar for primitive float types
impl Scalar for f32 {
    #[inline]
    fn abs(self) -> Self {
        f32::abs(self)
    }

    #[inline]
    fn signum(self) -> Self {
        f32::signum(self)
    }

    #[inline]
    fn powi(self, exp: i32) -> Self {
        f32::powi(self, exp)
    }
}

impl Scalar for f64 {
    #[inline]
    fn abs(self) -> Self {
        f64::abs(self)
    }

    #[inline]
    fn signum(self) -> Self {
        f64::signum(self)
    }

    #[inline]
    fn powi(self, exp: i32) -> Self {
        f64::powi(self, exp)
    }
}

// Implement Real for primitive float types
impl Real for f32 {
    #[inline]
    fn sqrt(self) -> Self {
        f32::sqrt(self)
    }

    #[inline]
    fn cbrt(self) -> Self {
        f32::cbrt(self)
    }
}

impl Real for f64 {
    #[inline]
    fn sqrt(self) -> Self {
        f64::sqrt(self)
    }

    #[inline]
    fn cbrt(self) -> Self {
        f64::cbrt(self)
    }
}

// Implement Float for primitive float types
impl Float for f32 {
    #[inline]
    fn epsilon() -> Self {
        f32::EPSILON
    }

    #[inline]
    fn infinity() -> Self {
        f32::INFINITY
    }

    #[inline]
    fn neg_infinity() -> Self {
        f32::NEG_INFINITY
    }

    #[inline]
    fn nan() -> Self {
        f32::NAN
    }

    #[inline]
    fn is_nan(self) -> bool {
        f32::is_nan(self)
    }

    #[inline]
    fn is_infinite(self) -> bool {
        f32::is_infinite(self)
    }

    #[inline]
    fn is_finite(self) -> bool {
        f32::is_finite(self)
    }

    #[inline]
    fn is_normal(self) -> bool {
        f32::is_normal(self)
    }

    #[inline]
    fn exp(self) -> Self {
        f32::exp(self)
    }

    #[inline]
    fn exp2(self) -> Self {
        f32::exp2(self)
    }

    #[inline]
    fn ln(self) -> Self {
        f32::ln(self)
    }

    #[inline]
    fn log2(self) -> Self {
        f32::log2(self)
    }

    #[inline]
    fn log10(self) -> Self {
        f32::log10(self)
    }

    #[inline]
    fn sin(self) -> Self {
        f32::sin(self)
    }

    #[inline]
    fn cos(self) -> Self {
        f32::cos(self)
    }

    #[inline]
    fn asin(self) -> Self {
        f32::asin(self)
    }

    #[inline]
    fn acos(self) -> Self {
        f32::acos(self)
    }

    #[inline]
    fn atan(self) -> Self {
        f32::atan(self)
    }

    #[inline]
    fn atan2(self, other: Self) -> Self {
        f32::atan2(self, other)
    }

    #[inline]
    fn sinh(self) -> Self {
        f32::sinh(self)
    }

    #[inline]
    fn cosh(self) -> Self {
        f32::cosh(self)
    }

    #[inline]
    fn tanh(self) -> Self {
        f32::tanh(self)
    }

    #[inline]
    fn log(self, base: Self) -> Self {
        f32::ln(self) / f32::ln(base)
    }

    #[inline]
    fn tan(self) -> Self {
        f32::tan(self)
    }

    #[inline]
    fn powf(self, exp: Self) -> Self {
        f32::powf(self, exp)
    }

    #[inline]
    fn fract(self) -> Self {
        self - f32::trunc(self)
    }

    #[inline]
    fn clamp(self, min: Self, max: Self) -> Self {
        f32::max(f32::min(self, max), min)
    }

    #[inline]
    fn floor(self) -> Self {
        f32::floor(self)
    }

    #[inline]
    fn ceil(self) -> Self {
        f32::ceil(self)
    }

    #[inline]
    fn round(self) -> Self {
        f32::round(self)
    }

    #[inline]
    fn trunc(self) -> Self {
        f32::trunc(self)
    }

    #[inline]
    fn mul_add(self, a: Self, b: Self) -> Self {
        f32::mul_add(self, a, b)
    }

    #[inline]
    fn min(self, other: Self) -> Self {
        f32::min(self, other)
    }

    #[inline]
    fn max(self, other: Self) -> Self {
        f32::max(self, other)
    }
}

impl Float for f64 {
    #[inline]
    fn epsilon() -> Self {
        f64::EPSILON
    }

    #[inline]
    fn infinity() -> Self {
        f64::INFINITY
    }

    #[inline]
    fn neg_infinity() -> Self {
        f64::NEG_INFINITY
    }

    #[inline]
    fn nan() -> Self {
        f64::NAN
    }

    #[inline]
    fn is_nan(self) -> bool {
        f64::is_nan(self)
    }

    #[inline]
    fn is_infinite(self) -> bool {
        f64::is_infinite(self)
    }

    #[inline]
    fn is_finite(self) -> bool {
        f64::is_finite(self)
    }

    #[inline]
    fn is_normal(self) -> bool {
        f64::is_normal(self)
    }

    #[inline]
    fn exp(self) -> Self {
        f64::exp(self)
    }

    #[inline]
    fn exp2(self) -> Self {
        f64::exp2(self)
    }

    #[inline]
    fn ln(self) -> Self {
        f64::ln(self)
    }

    #[inline]
    fn log2(self) -> Self {
        f64::log2(self)
    }

    #[inline]
    fn log10(self) -> Self {
        f64::log10(self)
    }

    #[inline]
    fn sin(self) -> Self {
        f64::sin(self)
    }

    #[inline]
    fn cos(self) -> Self {
        f64::cos(self)
    }

    #[inline]
    fn asin(self) -> Self {
        f64::asin(self)
    }

    #[inline]
    fn acos(self) -> Self {
        f64::acos(self)
    }

    #[inline]
    fn atan(self) -> Self {
        f64::atan(self)
    }

    #[inline]
    fn atan2(self, other: Self) -> Self {
        f64::atan2(self, other)
    }

    #[inline]
    fn sinh(self) -> Self {
        f64::sinh(self)
    }

    #[inline]
    fn cosh(self) -> Self {
        f64::cosh(self)
    }

    #[inline]
    fn tanh(self) -> Self {
        f64::tanh(self)
    }

    #[inline]
    fn log(self, base: Self) -> Self {
        f64::ln(self) / f64::ln(base)
    }

    #[inline]
    fn tan(self) -> Self {
        f64::tan(self)
    }

    #[inline]
    fn powf(self, exp: Self) -> Self {
        f64::powf(self, exp)
    }

    #[inline]
    fn fract(self) -> Self {
        self - f64::trunc(self)
    }

    #[inline]
    fn clamp(self, min: Self, max: Self) -> Self {
        f64::max(f64::min(self, max), min)
    }

    #[inline]
    fn floor(self) -> Self {
        f64::floor(self)
    }

    #[inline]
    fn ceil(self) -> Self {
        f64::ceil(self)
    }

    #[inline]
    fn round(self) -> Self {
        f64::round(self)
    }

    #[inline]
    fn trunc(self) -> Self {
        f64::trunc(self)
    }

    #[inline]
    fn mul_add(self, a: Self, b: Self) -> Self {
        f64::mul_add(self, a, b)
    }

    #[inline]
    fn min(self, other: Self) -> Self {
        f64::min(self, other)
    }

    #[inline]
    fn max(self, other: Self) -> Self {
        f64::max(self, other)
    }
}

/// Trait alias for numeric operations without assignment.
pub trait NumberOps<Rhs = Self, Output = Self>:
    Add<Rhs, Output = Output> + Sub<Rhs, Output = Output> + Mul<Rhs, Output = Output> + Div<Rhs, Output = Output>
{
}

impl<T, Rhs, Output> NumberOps<Rhs, Output> for T where
    T: Add<Rhs, Output = Output>
        + Sub<Rhs, Output = Output>
        + Mul<Rhs, Output = Output>
        + Div<Rhs, Output = Output>
{
}

/// Trait alias for float operations including transcendental functions.
pub trait FloatOps: Float + Neg<Output = Self> {}

impl<T> FloatOps for T where T: Float + Neg<Output = Self> {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scalar_operations_f32() {
        let x = 4.0_f32;
        assert_eq!(x.abs(), 4.0);
        assert_eq!((-x).abs(), 4.0);
        assert_eq!(x.signum(), 1.0);
        assert_eq!(x.powi(2), 16.0);
    }

    #[test]
    fn test_scalar_operations_f64() {
        let x = 4.0_f64;
        assert_eq!(x.abs(), 4.0);
        assert_eq!((-x).abs(), 4.0);
        assert_eq!(x.signum(), 1.0);
        assert_eq!(x.powi(2), 16.0);
    }

    #[test]
    fn test_real_operations() {
        let x = 9.0_f64;
        assert_eq!(x.sqrt(), 3.0);
        assert_eq!(x.square(), 81.0);
        assert_eq!(x.recip(), 1.0 / 9.0);
    }

    #[test]
    fn test_float_special_values() {
        use num_traits::Float as NF;
        let nan: f64 = NF::nan();
        let inf: f64 = NF::infinity();
        let neg_inf: f64 = NF::neg_infinity();
        assert!(nan.is_nan());
        assert!(inf.is_infinite());
        assert!(neg_inf.is_infinite());
        assert!(!(0.0_f64).is_infinite());
        assert!((1.0_f64).is_finite());
    }

    #[test]
    fn test_trigonometric_identities() {
        use std::f64::consts::PI;

        let x = PI / 4.0;
        let sin_x = x.sin();
        let cos_x = x.cos();

        // sin²(x) + cos²(x) = 1
        let result = sin_x * sin_x + cos_x * cos_x;
        assert!((result - 1.0).abs() < 1e-10);

        // tan(x) = sin(x) / cos(x)
        assert!((x.tan() - sin_x / cos_x).abs() < 1e-10);
    }

    #[test]
    fn test_exponential_logarithm_inverse() {
        let x = 2.0_f64;
        assert!((x.exp().ln() - x).abs() < 1e-10);
        assert!((x.ln().exp() - x).abs() < 1e-10);
    }
}
