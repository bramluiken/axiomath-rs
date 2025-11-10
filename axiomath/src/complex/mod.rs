//! Complex number arithmetic and operations.
//!
//! This module provides a comprehensive implementation of complex numbers, including
//! arithmetic operations, polar form conversions, and complex-valued functions.
//!
//! # Theory
//!
//! A complex number is a number of the form z = a + bi, where:
//! - a is the real part
//! - b is the imaginary part
//! - i is the imaginary unit satisfying i² = -1
//!
//! Complex numbers extend the real number system and are essential in many areas
//! of mathematics, physics, and engineering, including:
//! - Signal processing and Fourier analysis
//! - Quantum mechanics
//! - Fluid dynamics
//! - Control theory
//! - Electrical engineering
//!
//! # Mathematical Properties
//!
//! Complex numbers form a field under addition and multiplication:
//! - **Addition**: (a + bi) + (c + di) = (a + c) + (b + d)i
//! - **Multiplication**: (a + bi)(c + di) = (ac - bd) + (ad + bc)i
//! - **Conjugate**: z̄ = a - bi
//! - **Magnitude**: |z| = √(a² + b²)
//! - **Argument**: arg(z) = atan2(b, a)
//!
//! # Polar Form
//!
//! Any complex number can be represented in polar form:
//! - z = r(cos θ + i sin θ) = re^(iθ)
//! where r = |z| is the magnitude and θ = arg(z) is the argument.
//!
//! # Examples
//!
//! ```
//! use axiomath::complex::{Complex, add, multiply, conjugate, magnitude};
//!
//! // Create complex numbers
//! let z1 = Complex::new(3.0, 4.0);  // 3 + 4i
//! let z2 = Complex::new(1.0, 2.0);  // 1 + 2i
//!
//! // Basic arithmetic
//! let sum = add(&z1, &z2);          // 4 + 6i
//! let product = multiply(&z1, &z2); // -5 + 10i
//!
//! // Properties
//! let conj = conjugate(&z1);        // 3 - 4i
//! let mag = magnitude(&z1);         // 5.0
//!
//! // Using operator overloads
//! let z3 = z1 + z2;
//! let z4 = z1 * z2;
//! ```
//!
//! # References
//!
//! - Brown, J. W., & Churchill, R. V. (2013). *Complex Variables and Applications*. McGraw-Hill.
//! - Needham, T. (1997). *Visual Complex Analysis*. Oxford University Press.
//! - Ahlfors, L. V. (1979). *Complex Analysis*. McGraw-Hill.

use crate::core::{Float, Scalar};
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// A complex number with real and imaginary parts.
///
/// The generic type `T` should implement [`Scalar`] for basic arithmetic operations,
/// or [`Float`] for transcendental functions.
///
/// # Examples
///
/// ```
/// use axiomath::complex::Complex;
///
/// let z = Complex::new(3.0, 4.0);
/// assert_eq!(z.real, 3.0);
/// assert_eq!(z.imag, 4.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Complex<T> {
    /// The real part of the complex number
    pub real: T,
    /// The imaginary part of the complex number
    pub imag: T,
}

impl<T: Scalar> Complex<T> {
    /// Creates a new complex number from real and imaginary parts.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::complex::Complex;
    ///
    /// let z = Complex::new(3.0, 4.0);
    /// assert_eq!(z.real, 3.0);
    /// assert_eq!(z.imag, 4.0);
    /// ```
    #[inline]
    pub fn new(real: T, imag: T) -> Self {
        Complex { real, imag }
    }

    /// Creates a complex number representing zero (0 + 0i).
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::complex::Complex;
    ///
    /// let z: Complex<f64> = Complex::zero();
    /// assert_eq!(z.real, 0.0);
    /// assert_eq!(z.imag, 0.0);
    /// ```
    #[inline]
    pub fn zero() -> Self {
        Complex {
            real: T::zero(),
            imag: T::zero(),
        }
    }

    /// Creates a complex number representing one (1 + 0i).
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::complex::Complex;
    ///
    /// let z: Complex<f64> = Complex::one();
    /// assert_eq!(z.real, 1.0);
    /// assert_eq!(z.imag, 0.0);
    /// ```
    #[inline]
    pub fn one() -> Self {
        Complex {
            real: T::one(),
            imag: T::zero(),
        }
    }

    /// Creates a complex number representing the imaginary unit (0 + 1i).
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::complex::Complex;
    ///
    /// let i: Complex<f64> = Complex::i();
    /// assert_eq!(i.real, 0.0);
    /// assert_eq!(i.imag, 1.0);
    /// ```
    #[inline]
    pub fn i() -> Self {
        Complex {
            real: T::zero(),
            imag: T::one(),
        }
    }

    /// Creates a complex number from a real value (real + 0i).
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::complex::Complex;
    ///
    /// let z = Complex::from_real(5.0);
    /// assert_eq!(z.real, 5.0);
    /// assert_eq!(z.imag, 0.0);
    /// ```
    #[inline]
    pub fn from_real(real: T) -> Self {
        Complex {
            real,
            imag: T::zero(),
        }
    }

    /// Creates a pure imaginary number (0 + imag*i).
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::complex::Complex;
    ///
    /// let z = Complex::from_imag(3.0);
    /// assert_eq!(z.real, 0.0);
    /// assert_eq!(z.imag, 3.0);
    /// ```
    #[inline]
    pub fn from_imag(imag: T) -> Self {
        Complex {
            real: T::zero(),
            imag,
        }
    }

    /// Checks if the complex number is zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::complex::Complex;
    ///
    /// let z = Complex::zero();
    /// assert!(z.is_zero());
    ///
    /// let w = Complex::new(1.0, 2.0);
    /// assert!(!w.is_zero());
    /// ```
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.real.is_zero() && self.imag.is_zero()
    }
}

impl<T: Float> Complex<T> {
    /// Checks if either component is NaN.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::complex::Complex;
    ///
    /// let z = Complex::new(1.0, f64::NAN);
    /// assert!(z.is_nan());
    /// ```
    #[inline]
    pub fn is_nan(&self) -> bool {
        self.real.is_nan() || self.imag.is_nan()
    }

    /// Checks if either component is infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::complex::Complex;
    ///
    /// let z = Complex::new(f64::INFINITY, 2.0);
    /// assert!(z.is_infinite());
    /// ```
    #[inline]
    pub fn is_infinite(&self) -> bool {
        self.real.is_infinite() || self.imag.is_infinite()
    }

    /// Checks if both components are finite (not infinite and not NaN).
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::complex::Complex;
    ///
    /// let z = Complex::new(1.0, 2.0);
    /// assert!(z.is_finite());
    ///
    /// let w = Complex::new(f64::INFINITY, 2.0);
    /// assert!(!w.is_finite());
    /// ```
    #[inline]
    pub fn is_finite(&self) -> bool {
        self.real.is_finite() && self.imag.is_finite()
    }
}

impl<T: fmt::Display> fmt::Display for Complex<T>
where
    T: Scalar,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.imag < T::zero() {
            write!(f, "{} - {}i", self.real, self.imag.abs())
        } else {
            write!(f, "{} + {}i", self.real, self.imag)
        }
    }
}

// ================================================================================================
// Basic Arithmetic Operations
// ================================================================================================

/// Adds two complex numbers.
///
/// # Mathematical Definition
///
/// ```text
/// (a + bi) + (c + di) = (a + c) + (b + d)i
/// ```
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::complex::{Complex, add};
///
/// let z1 = Complex::new(3.0, 4.0);
/// let z2 = Complex::new(1.0, 2.0);
/// let result = add(&z1, &z2);
///
/// assert_eq!(result.real, 4.0);
/// assert_eq!(result.imag, 6.0);
/// ```
#[inline]
pub fn add<T: Scalar>(a: &Complex<T>, b: &Complex<T>) -> Complex<T> {
    Complex {
        real: a.real + b.real,
        imag: a.imag + b.imag,
    }
}

/// Subtracts one complex number from another.
///
/// # Mathematical Definition
///
/// ```text
/// (a + bi) - (c + di) = (a - c) + (b - d)i
/// ```
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::complex::{Complex, subtract};
///
/// let z1 = Complex::new(5.0, 7.0);
/// let z2 = Complex::new(2.0, 3.0);
/// let result = subtract(&z1, &z2);
///
/// assert_eq!(result.real, 3.0);
/// assert_eq!(result.imag, 4.0);
/// ```
#[inline]
pub fn subtract<T: Scalar>(a: &Complex<T>, b: &Complex<T>) -> Complex<T> {
    Complex {
        real: a.real - b.real,
        imag: a.imag - b.imag,
    }
}

/// Multiplies two complex numbers.
///
/// # Mathematical Definition
///
/// ```text
/// (a + bi)(c + di) = (ac - bd) + (ad + bc)i
/// ```
///
/// This follows from the distributive property and the fact that i² = -1.
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::complex::{Complex, multiply};
///
/// let z1 = Complex::new(3.0, 4.0);
/// let z2 = Complex::new(1.0, 2.0);
/// let result = multiply(&z1, &z2);
///
/// // (3 + 4i)(1 + 2i) = 3 + 6i + 4i + 8i² = 3 + 10i - 8 = -5 + 10i
/// assert_eq!(result.real, -5.0);
/// assert_eq!(result.imag, 10.0);
/// ```
///
/// # References
///
/// - Brown & Churchill, "Complex Variables and Applications", Section 3
#[inline]
pub fn multiply<T: Scalar>(a: &Complex<T>, b: &Complex<T>) -> Complex<T> {
    Complex {
        real: a.real * b.real - a.imag * b.imag,
        imag: a.real * b.imag + a.imag * b.real,
    }
}

/// Divides one complex number by another.
///
/// # Mathematical Definition
///
/// ```text
/// (a + bi) / (c + di) = [(a + bi)(c - di)] / (c² + d²)
///                     = [(ac + bd) + (bc - ad)i] / (c² + d²)
/// ```
///
/// This uses Smith's formula for numerical stability, which avoids overflow
/// when |c| and |d| are large.
///
/// # Algorithm
///
/// 1. If |d| ≤ |c|:
///    - r = d/c
///    - den = c + d*r
///    - real = (a + b*r) / den
///    - imag = (b - a*r) / den
/// 2. Otherwise:
///    - r = c/d
///    - den = d + c*r
///    - real = (a*r + b) / den
///    - imag = (b*r - a) / den
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::complex::{Complex, divide};
/// use approx::assert_abs_diff_eq;
///
/// let z1 = Complex::new(1.0, 1.0);
/// let z2 = Complex::new(1.0, 0.0);
/// let result = divide(&z1, &z2);
///
/// assert_abs_diff_eq!(result.real, 1.0, epsilon = 1e-10);
/// assert_abs_diff_eq!(result.imag, 1.0, epsilon = 1e-10);
/// ```
///
/// # References
///
/// - Smith, R. L. (1962). "Algorithm 116: Complex division". Communications of the ACM.
pub fn divide<T: Float>(a: &Complex<T>, b: &Complex<T>) -> Complex<T> {
    // Handle special cases
    if b.is_zero() {
        return Complex {
            real: T::infinity(),
            imag: T::infinity(),
        };
    }

    if a.is_zero() {
        return Complex::zero();
    }

    // Smith's algorithm for numerical stability
    if b.imag.abs() <= b.real.abs() {
        let r = b.imag / b.real;
        let den = b.real + b.imag * r;
        Complex {
            real: (a.real + a.imag * r) / den,
            imag: (a.imag - a.real * r) / den,
        }
    } else {
        let r = b.real / b.imag;
        let den = b.imag + b.real * r;
        Complex {
            real: (a.real * r + a.imag) / den,
            imag: (a.imag * r - a.real) / den,
        }
    }
}

/// Negates a complex number.
///
/// # Mathematical Definition
///
/// ```text
/// -(a + bi) = -a - bi
/// ```
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::complex::{Complex, negate};
///
/// let z = Complex::new(3.0, 4.0);
/// let result = negate(&z);
///
/// assert_eq!(result.real, -3.0);
/// assert_eq!(result.imag, -4.0);
/// ```
#[inline]
pub fn negate<T: Scalar>(z: &Complex<T>) -> Complex<T>
where
    T: Neg<Output = T>,
{
    Complex {
        real: -z.real,
        imag: -z.imag,
    }
}

/// Returns the complex conjugate.
///
/// # Mathematical Definition
///
/// ```text
/// z̄ = a - bi  (for z = a + bi)
/// ```
///
/// The conjugate reflects the complex number across the real axis.
///
/// # Properties
///
/// - z + z̄ = 2*Re(z)
/// - z * z̄ = |z|²
/// - (z̄)̄ = z
/// - (z₁ + z₂)̄ = z̄₁ + z̄₂
/// - (z₁ * z₂)̄ = z̄₁ * z̄₂
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::complex::{Complex, conjugate};
///
/// let z = Complex::new(3.0, 4.0);
/// let conj = conjugate(&z);
///
/// assert_eq!(conj.real, 3.0);
/// assert_eq!(conj.imag, -4.0);
/// ```
#[inline]
pub fn conjugate<T: Scalar>(z: &Complex<T>) -> Complex<T>
where
    T: Neg<Output = T>,
{
    Complex {
        real: z.real,
        imag: -z.imag,
    }
}

/// Returns the reciprocal of a complex number (1/z).
///
/// # Mathematical Definition
///
/// ```text
/// 1/z = z̄ / |z|² = (a - bi) / (a² + b²)
/// ```
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::complex::{Complex, reciprocal};
/// use approx::assert_abs_diff_eq;
///
/// let z = Complex::new(3.0, 4.0);
/// let recip = reciprocal(&z);
///
/// // 1/(3+4i) = (3-4i)/25 = 0.12 - 0.16i
/// assert_abs_diff_eq!(recip.real, 0.12, epsilon = 1e-10);
/// assert_abs_diff_eq!(recip.imag, -0.16, epsilon = 1e-10);
/// ```
pub fn reciprocal<T: Float>(z: &Complex<T>) -> Complex<T> {
    let one = Complex::one();
    divide(&one, z)
}

// ================================================================================================
// Polar Form and Magnitude
// ================================================================================================

/// Returns the magnitude (absolute value, modulus) of a complex number.
///
/// # Mathematical Definition
///
/// ```text
/// |z| = √(a² + b²)  for z = a + bi
/// ```
///
/// The magnitude represents the distance from the origin in the complex plane.
///
/// # Algorithm
///
/// Uses a numerically stable computation similar to `hypot(a, b)` to avoid
/// overflow/underflow for large or small values:
///
/// ```text
/// if |a| ≥ |b|:
///     |z| = |a| * √(1 + (b/a)²)
/// else:
///     |z| = |b| * √(1 + (a/b)²)
/// ```
///
/// # Properties
///
/// - |z| ≥ 0
/// - |z| = 0 if and only if z = 0
/// - |z̄| = |z|
/// - |z₁ * z₂| = |z₁| * |z₂|
/// - |z₁ / z₂| = |z₁| / |z₂|
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::complex::{Complex, magnitude};
///
/// let z = Complex::new(3.0, 4.0);
/// let mag = magnitude(&z);
/// assert_eq!(mag, 5.0); // √(3² + 4²) = √25 = 5
/// ```
///
/// # References
///
/// - Brown & Churchill, "Complex Variables and Applications", Section 4
pub fn magnitude<T: Float>(z: &Complex<T>) -> T {
    // Handle special cases
    if z.is_zero() {
        return T::zero();
    }

    let a = z.real.abs();
    let b = z.imag.abs();

    // Use stable computation to avoid overflow
    if a >= b {
        let r = b / a;
        a * (T::one() + r * r).sqrt()
    } else {
        let r = a / b;
        b * (T::one() + r * r).sqrt()
    }
}

/// Returns the argument (phase angle) of a complex number in radians.
///
/// # Mathematical Definition
///
/// ```text
/// arg(z) = atan2(b, a)  for z = a + bi
/// ```
///
/// The argument is the angle θ in the polar representation z = re^(iθ).
/// The value returned is in the range (-π, π].
///
/// # Properties
///
/// - arg(0) is undefined (returns 0 by convention)
/// - arg(z̄) = -arg(z)
/// - arg(z₁ * z₂) = arg(z₁) + arg(z₂) (modulo 2π)
/// - arg(z₁ / z₂) = arg(z₁) - arg(z₂) (modulo 2π)
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::complex::{Complex, argument};
/// use std::f64::consts::PI;
///
/// let z = Complex::new(1.0, 1.0);
/// let arg = argument(&z);
/// assert!((arg - PI/4.0).abs() < 1e-10); // 45 degrees
/// ```
///
/// # References
///
/// - Brown & Churchill, "Complex Variables and Applications", Section 7
pub fn argument<T: Float>(z: &Complex<T>) -> T {
    z.imag.atan2(z.real)
}

/// Converts polar coordinates (r, θ) to a complex number.
///
/// # Mathematical Definition
///
/// ```text
/// z = r * e^(iθ) = r * (cos θ + i sin θ)
/// ```
///
/// where r is the magnitude and θ is the argument (phase angle).
///
/// # Arguments
///
/// * `r` - The magnitude (must be non-negative)
/// * `theta` - The argument in radians
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::complex::{from_polar, Complex};
/// use std::f64::consts::PI;
/// use approx::assert_abs_diff_eq;
///
/// let z = from_polar(5.0, PI/2.0);
/// assert_abs_diff_eq!(z.real, 0.0, epsilon = 1e-10);
/// assert_abs_diff_eq!(z.imag, 5.0, epsilon = 1e-10);
/// ```
///
/// # References
///
/// - Brown & Churchill, "Complex Variables and Applications", Section 7
pub fn from_polar<T: Float>(r: T, theta: T) -> Complex<T> {
    Complex {
        real: r * theta.cos(),
        imag: r * theta.sin(),
    }
}

/// Converts a complex number to polar form (magnitude, argument).
///
/// # Mathematical Definition
///
/// For z = a + bi, returns (r, θ) where:
/// - r = |z| = √(a² + b²)
/// - θ = arg(z) = atan2(b, a)
///
/// # Returns
///
/// A tuple `(magnitude, argument)` where:
/// - magnitude is the distance from the origin
/// - argument is the angle in radians, in the range (-π, π]
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::complex::{Complex, to_polar};
/// use std::f64::consts::PI;
/// use approx::assert_abs_diff_eq;
///
/// let z = Complex::new(3.0, 4.0);
/// let (r, theta) = to_polar(&z);
///
/// assert_abs_diff_eq!(r, 5.0, epsilon = 1e-10);
/// assert_abs_diff_eq!(theta, 0.927295218, epsilon = 1e-9);
/// ```
pub fn to_polar<T: Float>(z: &Complex<T>) -> (T, T) {
    (magnitude(z), argument(z))
}

/// Multiplies two complex numbers in polar form (more efficient).
///
/// # Mathematical Definition
///
/// For z₁ = r₁e^(iθ₁) and z₂ = r₂e^(iθ₂):
///
/// ```text
/// z₁ * z₂ = (r₁ * r₂) * e^(i(θ₁ + θ₂))
/// ```
///
/// This function is more efficient than regular multiplication when the
/// complex numbers are already in polar form, as it avoids conversion.
///
/// # Arguments
///
/// * `r1` - Magnitude of first complex number
/// * `theta1` - Argument of first complex number (radians)
/// * `r2` - Magnitude of second complex number
/// * `theta2` - Argument of second complex number (radians)
///
/// # Returns
///
/// The product as a complex number in rectangular form.
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::complex::multiply_polar;
/// use std::f64::consts::PI;
/// use approx::assert_abs_diff_eq;
///
/// // Multiply two unit complex numbers at 45° each = 90°
/// let result = multiply_polar(1.0, PI/4.0, 1.0, PI/4.0);
///
/// assert_abs_diff_eq!(result.real, 0.0, epsilon = 1e-10);
/// assert_abs_diff_eq!(result.imag, 1.0, epsilon = 1e-10);
/// ```
pub fn multiply_polar<T: Float>(r1: T, theta1: T, r2: T, theta2: T) -> Complex<T> {
    from_polar(r1 * r2, theta1 + theta2)
}

// ================================================================================================
// Complex Functions
// ================================================================================================

/// Computes the complex exponential function e^z.
///
/// # Mathematical Definition
///
/// For z = a + bi:
///
/// ```text
/// e^z = e^(a+bi) = e^a * e^(bi) = e^a * (cos b + i sin b)
/// ```
///
/// This is Euler's formula extended to complex exponents.
///
/// # Properties
///
/// - e^0 = 1
/// - e^(z₁ + z₂) = e^z₁ * e^z₂
/// - |e^(iy)| = 1 for real y
/// - e^(iπ) = -1 (Euler's identity)
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::complex::{Complex, complex_exp};
/// use std::f64::consts::PI;
/// use approx::assert_abs_diff_eq;
///
/// // e^(iπ) = -1
/// let z = Complex::new(0.0, PI);
/// let result = complex_exp(&z);
/// assert_abs_diff_eq!(result.real, -1.0, epsilon = 1e-10);
/// assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-10);
/// ```
///
/// # References
///
/// - NIST DLMF, Section 4.2: "Exponential Function"
/// - Needham, "Visual Complex Analysis", Chapter 1
pub fn complex_exp<T: Float>(z: &Complex<T>) -> Complex<T> {
    let exp_real = z.real.exp();
    Complex {
        real: exp_real * z.imag.cos(),
        imag: exp_real * z.imag.sin(),
    }
}

/// Computes the principal branch of the complex natural logarithm.
///
/// # Mathematical Definition
///
/// For z = a + bi:
///
/// ```text
/// ln(z) = ln|z| + i*arg(z)
/// ```
///
/// where arg(z) is the principal argument in the range (-π, π].
///
/// # Branch Cut
///
/// The principal branch has a branch cut along the negative real axis.
/// The logarithm is discontinuous across this cut.
///
/// # Properties
///
/// - ln(e^z) = z + 2πik for some integer k
/// - ln(z₁ * z₂) = ln(z₁) + ln(z₂) + 2πik
/// - ln(1) = 0
/// - ln(e) = 1
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::complex::{Complex, complex_ln};
/// use approx::assert_abs_diff_eq;
///
/// let z = Complex::new(std::f64::consts::E, 0.0);
/// let result = complex_ln(&z);
/// assert_abs_diff_eq!(result.real, 1.0, epsilon = 1e-10);
/// assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-10);
/// ```
///
/// # References
///
/// - NIST DLMF, Section 4.2: "Logarithm"
/// - Brown & Churchill, "Complex Variables and Applications", Section 31
pub fn complex_ln<T: Float>(z: &Complex<T>) -> Complex<T> {
    Complex {
        real: magnitude(z).ln(),
        imag: argument(z),
    }
}

/// Computes the principal square root of a complex number.
///
/// # Mathematical Definition
///
/// For z = a + bi, the principal square root is the square root with
/// non-negative real part (or positive imaginary part if real part is zero).
///
/// # Algorithm
///
/// Using the formula:
///
/// ```text
/// √z = √[(|z| + a)/2] + i * sign(b) * √[(|z| - a)/2]
/// ```
///
/// where sign(b) = 1 if b ≥ 0, else -1.
///
/// # Branch Cut
///
/// The principal branch has a branch cut along the negative real axis.
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::complex::{Complex, complex_sqrt};
/// use approx::assert_abs_diff_eq;
///
/// let z = Complex::new(-1.0, 0.0);
/// let result = complex_sqrt(&z);
/// assert_abs_diff_eq!(result.real, 0.0, epsilon = 1e-10);
/// assert_abs_diff_eq!(result.imag, 1.0, epsilon = 1e-10);
/// ```
///
/// # References
///
/// - NIST DLMF, Section 4.7: "Square Root"
pub fn complex_sqrt<T: Float>(z: &Complex<T>) -> Complex<T> {
    if z.is_zero() {
        return Complex::zero();
    }

    let mag = magnitude(z);
    let two = T::one() + T::one();

    let real_part = ((mag + z.real) / two).sqrt();
    let imag_part = ((mag - z.real) / two).sqrt();

    if z.imag >= T::zero() {
        Complex {
            real: real_part,
            imag: imag_part,
        }
    } else {
        Complex {
            real: real_part,
            imag: -imag_part,
        }
    }
}

/// Computes z raised to the power w: z^w.
///
/// # Mathematical Definition
///
/// For complex z and w:
///
/// ```text
/// z^w = e^(w * ln(z))
/// ```
///
/// Uses the principal branch of the logarithm.
///
/// # Special Cases
///
/// - 0^0 returns 1 (by convention)
/// - 0^w returns 0 for Re(w) > 0
/// - 0^w returns infinity for Re(w) < 0
///
/// # Branch Cut
///
/// Inherits the branch cut from complex_ln along the negative real axis for z.
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::complex::{Complex, complex_pow};
/// use approx::assert_abs_diff_eq;
///
/// // i^i = e^(-π/2)
/// let i = Complex::new(0.0, 1.0);
/// let result = complex_pow(&i, &i);
/// assert_abs_diff_eq!(result.real, 0.207879576, epsilon = 1e-9);
/// assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-9);
/// ```
///
/// # References
///
/// - NIST DLMF, Section 4.2: "Powers"
pub fn complex_pow<T: Float>(z: &Complex<T>, w: &Complex<T>) -> Complex<T> {
    // Handle special cases
    if z.is_zero() {
        if w.is_zero() {
            return Complex::one(); // 0^0 = 1 by convention
        } else if w.real > T::zero() {
            return Complex::zero();
        } else {
            return Complex {
                real: T::infinity(),
                imag: T::zero(),
            };
        }
    }

    if w.is_zero() {
        return Complex::one();
    }

    // z^w = e^(w * ln(z))
    let ln_z = complex_ln(z);
    let w_ln_z = multiply(w, &ln_z);
    complex_exp(&w_ln_z)
}

/// Computes the complex sine function.
///
/// # Mathematical Definition
///
/// For z = a + bi:
///
/// ```text
/// sin(z) = (e^(iz) - e^(-iz)) / (2i)
/// ```
///
/// Equivalently:
///
/// ```text
/// sin(a + bi) = sin(a)cosh(b) + i*cos(a)sinh(b)
/// ```
///
/// # Properties
///
/// - sin(-z) = -sin(z) (odd function)
/// - sin(z + 2π) = sin(z) (periodic)
/// - |sin(z)| can exceed 1 for complex z
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::complex::{Complex, complex_sin};
/// use std::f64::consts::PI;
/// use approx::assert_abs_diff_eq;
///
/// let z = Complex::new(PI / 2.0, 0.0);
/// let result = complex_sin(&z);
/// assert_abs_diff_eq!(result.real, 1.0, epsilon = 1e-10);
/// assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-10);
/// ```
///
/// # References
///
/// - NIST DLMF, Section 4.14: "Trigonometric Functions"
pub fn complex_sin<T: Float>(z: &Complex<T>) -> Complex<T> {
    let sin_a = z.real.sin();
    let cos_a = z.real.cos();
    let sinh_b = z.imag.sinh();
    let cosh_b = z.imag.cosh();

    Complex {
        real: sin_a * cosh_b,
        imag: cos_a * sinh_b,
    }
}

/// Computes the complex cosine function.
///
/// # Mathematical Definition
///
/// For z = a + bi:
///
/// ```text
/// cos(z) = (e^(iz) + e^(-iz)) / 2
/// ```
///
/// Equivalently:
///
/// ```text
/// cos(a + bi) = cos(a)cosh(b) - i*sin(a)sinh(b)
/// ```
///
/// # Properties
///
/// - cos(-z) = cos(z) (even function)
/// - cos(z + 2π) = cos(z) (periodic)
/// - |cos(z)| can exceed 1 for complex z
/// - sin²(z) + cos²(z) = 1 for all complex z
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::complex::{Complex, complex_cos};
/// use std::f64::consts::PI;
/// use approx::assert_abs_diff_eq;
///
/// let z = Complex::new(0.0, 0.0);
/// let result = complex_cos(&z);
/// assert_abs_diff_eq!(result.real, 1.0, epsilon = 1e-10);
/// assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-10);
/// ```
///
/// # References
///
/// - NIST DLMF, Section 4.14: "Trigonometric Functions"
pub fn complex_cos<T: Float>(z: &Complex<T>) -> Complex<T> {
    let sin_a = z.real.sin();
    let cos_a = z.real.cos();
    let sinh_b = z.imag.sinh();
    let cosh_b = z.imag.cosh();

    Complex {
        real: cos_a * cosh_b,
        imag: -sin_a * sinh_b,
    }
}

// ================================================================================================
// Operator Overloading
// ================================================================================================

// Addition: Complex + Complex
impl<T: Scalar> Add for Complex<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        add(&self, &rhs)
    }
}

impl<T: Scalar> Add<&Complex<T>> for Complex<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: &Complex<T>) -> Self::Output {
        add(&self, rhs)
    }
}

impl<T: Scalar> Add<Complex<T>> for &Complex<T> {
    type Output = Complex<T>;

    #[inline]
    fn add(self, rhs: Complex<T>) -> Self::Output {
        add(self, &rhs)
    }
}

impl<T: Scalar> Add<&Complex<T>> for &Complex<T> {
    type Output = Complex<T>;

    #[inline]
    fn add(self, rhs: &Complex<T>) -> Self::Output {
        add(self, rhs)
    }
}

// Subtraction: Complex - Complex
impl<T: Scalar> Sub for Complex<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        subtract(&self, &rhs)
    }
}

impl<T: Scalar> Sub<&Complex<T>> for Complex<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: &Complex<T>) -> Self::Output {
        subtract(&self, rhs)
    }
}

impl<T: Scalar> Sub<Complex<T>> for &Complex<T> {
    type Output = Complex<T>;

    #[inline]
    fn sub(self, rhs: Complex<T>) -> Self::Output {
        subtract(self, &rhs)
    }
}

impl<T: Scalar> Sub<&Complex<T>> for &Complex<T> {
    type Output = Complex<T>;

    #[inline]
    fn sub(self, rhs: &Complex<T>) -> Self::Output {
        subtract(self, rhs)
    }
}

// Multiplication: Complex * Complex
impl<T: Scalar> Mul for Complex<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        multiply(&self, &rhs)
    }
}

impl<T: Scalar> Mul<&Complex<T>> for Complex<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: &Complex<T>) -> Self::Output {
        multiply(&self, rhs)
    }
}

impl<T: Scalar> Mul<Complex<T>> for &Complex<T> {
    type Output = Complex<T>;

    #[inline]
    fn mul(self, rhs: Complex<T>) -> Self::Output {
        multiply(self, &rhs)
    }
}

impl<T: Scalar> Mul<&Complex<T>> for &Complex<T> {
    type Output = Complex<T>;

    #[inline]
    fn mul(self, rhs: &Complex<T>) -> Self::Output {
        multiply(self, rhs)
    }
}

// Division: Complex / Complex
impl<T: Float> Div for Complex<T> {
    type Output = Self;

    #[inline]
    fn div(self, rhs: Self) -> Self::Output {
        divide(&self, &rhs)
    }
}

impl<T: Float> Div<&Complex<T>> for Complex<T> {
    type Output = Self;

    #[inline]
    fn div(self, rhs: &Complex<T>) -> Self::Output {
        divide(&self, rhs)
    }
}

impl<T: Float> Div<Complex<T>> for &Complex<T> {
    type Output = Complex<T>;

    #[inline]
    fn div(self, rhs: Complex<T>) -> Self::Output {
        divide(self, &rhs)
    }
}

impl<T: Float> Div<&Complex<T>> for &Complex<T> {
    type Output = Complex<T>;

    #[inline]
    fn div(self, rhs: &Complex<T>) -> Self::Output {
        divide(self, rhs)
    }
}

// Negation: -Complex
impl<T: Scalar> Neg for Complex<T>
where
    T: Neg<Output = T>,
{
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        negate(&self)
    }
}

impl<T: Scalar> Neg for &Complex<T>
where
    T: Neg<Output = T>,
{
    type Output = Complex<T>;

    #[inline]
    fn neg(self) -> Self::Output {
        negate(self)
    }
}

// ================================================================================================
// Tests
// ================================================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_new() {
        let z = Complex::new(3.0, 4.0);
        assert_eq!(z.real, 3.0);
        assert_eq!(z.imag, 4.0);
    }

    #[test]
    fn test_zero() {
        let z: Complex<f64> = Complex::zero();
        assert_eq!(z.real, 0.0);
        assert_eq!(z.imag, 0.0);
        assert!(z.is_zero());
    }

    #[test]
    fn test_one() {
        let z: Complex<f64> = Complex::one();
        assert_eq!(z.real, 1.0);
        assert_eq!(z.imag, 0.0);
    }

    #[test]
    fn test_i() {
        let z: Complex<f64> = Complex::i();
        assert_eq!(z.real, 0.0);
        assert_eq!(z.imag, 1.0);
    }

    #[test]
    fn test_from_real() {
        let z = Complex::from_real(5.0);
        assert_eq!(z.real, 5.0);
        assert_eq!(z.imag, 0.0);
    }

    #[test]
    fn test_from_imag() {
        let z = Complex::from_imag(3.0);
        assert_eq!(z.real, 0.0);
        assert_eq!(z.imag, 3.0);
    }

    #[test]
    fn test_add() {
        let z1 = Complex::new(3.0, 4.0);
        let z2 = Complex::new(1.0, 2.0);
        let result = add(&z1, &z2);
        assert_eq!(result.real, 4.0);
        assert_eq!(result.imag, 6.0);
    }

    #[test]
    fn test_add_operator() {
        let z1 = Complex::new(3.0, 4.0);
        let z2 = Complex::new(1.0, 2.0);
        let result = z1 + z2;
        assert_eq!(result.real, 4.0);
        assert_eq!(result.imag, 6.0);

        // Test all combinations
        let result = z1 + &z2;
        assert_eq!(result.real, 4.0);
        let result = &z1 + z2;
        assert_eq!(result.real, 4.0);
        let result = &z1 + &z2;
        assert_eq!(result.real, 4.0);
    }

    #[test]
    fn test_subtract() {
        let z1 = Complex::new(5.0, 7.0);
        let z2 = Complex::new(2.0, 3.0);
        let result = subtract(&z1, &z2);
        assert_eq!(result.real, 3.0);
        assert_eq!(result.imag, 4.0);
    }

    #[test]
    fn test_subtract_operator() {
        let z1 = Complex::new(5.0, 7.0);
        let z2 = Complex::new(2.0, 3.0);
        let result = z1 - z2;
        assert_eq!(result.real, 3.0);
        assert_eq!(result.imag, 4.0);
    }

    #[test]
    fn test_multiply() {
        let z1 = Complex::new(3.0, 4.0);
        let z2 = Complex::new(1.0, 2.0);
        let result = multiply(&z1, &z2);
        // (3 + 4i)(1 + 2i) = 3 + 6i + 4i + 8i² = 3 + 10i - 8 = -5 + 10i
        assert_eq!(result.real, -5.0);
        assert_eq!(result.imag, 10.0);
    }

    #[test]
    fn test_multiply_operator() {
        let z1 = Complex::new(3.0, 4.0);
        let z2 = Complex::new(1.0, 2.0);
        let result = z1 * z2;
        assert_eq!(result.real, -5.0);
        assert_eq!(result.imag, 10.0);
    }

    #[test]
    fn test_multiply_i_squared() {
        // i * i = -1
        let i: Complex<f64> = Complex::i();
        let result = multiply(&i, &i);
        assert_abs_diff_eq!(result.real, -1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_divide() {
        let z1 = Complex::new(1.0, 1.0);
        let z2 = Complex::new(1.0, 0.0);
        let result = divide(&z1, &z2);
        assert_abs_diff_eq!(result.real, 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_divide_operator() {
        let z1 = Complex::new(1.0, 1.0);
        let z2 = Complex::new(1.0, 0.0);
        let result = z1 / z2;
        assert_abs_diff_eq!(result.real, 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_divide_complex() {
        // (4 + 2i) / (3 - i) = (4 + 2i)(3 + i) / ((3-i)(3+i))
        //                     = (12 + 4i + 6i + 2i²) / (9 + 1)
        //                     = (12 + 10i - 2) / 10
        //                     = (10 + 10i) / 10 = 1 + i
        let z1 = Complex::new(4.0, 2.0);
        let z2 = Complex::new(3.0, -1.0);
        let result = divide(&z1, &z2);
        assert_abs_diff_eq!(result.real, 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_divide_by_zero() {
        let z = Complex::new(1.0, 0.0);
        let zero = Complex::zero();
        let result = divide(&z, &zero);
        assert!(result.real.is_infinite());
    }

    #[test]
    fn test_negate() {
        let z = Complex::new(3.0, 4.0);
        let result = negate(&z);
        assert_eq!(result.real, -3.0);
        assert_eq!(result.imag, -4.0);
    }

    #[test]
    fn test_negate_operator() {
        let z = Complex::new(3.0, 4.0);
        let result = -z;
        assert_eq!(result.real, -3.0);
        assert_eq!(result.imag, -4.0);

        let result = -&z;
        assert_eq!(result.real, -3.0);
    }

    #[test]
    fn test_conjugate() {
        let z = Complex::new(3.0, 4.0);
        let conj = conjugate(&z);
        assert_eq!(conj.real, 3.0);
        assert_eq!(conj.imag, -4.0);
    }

    #[test]
    fn test_conjugate_twice() {
        // (z̄)̄ = z
        let z = Complex::new(3.0, 4.0);
        let conj = conjugate(&conjugate(&z));
        assert_eq!(conj.real, z.real);
        assert_eq!(conj.imag, z.imag);
    }

    #[test]
    fn test_conjugate_multiplication() {
        // z * z̄ should be real and equal to |z|²
        let z = Complex::new(3.0, 4.0);
        let conj = conjugate(&z);
        let product = multiply(&z, &conj);
        // |z|² = 3² + 4² = 25
        assert_abs_diff_eq!(product.real, 25.0, epsilon = 1e-10);
        assert_abs_diff_eq!(product.imag, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_reciprocal() {
        let z = Complex::new(3.0, 4.0);
        let recip = reciprocal(&z);
        // 1/(3+4i) = (3-4i)/25 = 0.12 - 0.16i
        assert_abs_diff_eq!(recip.real, 0.12, epsilon = 1e-10);
        assert_abs_diff_eq!(recip.imag, -0.16, epsilon = 1e-10);
    }

    #[test]
    fn test_reciprocal_identity() {
        // z * (1/z) = 1
        let z = Complex::new(3.0, 4.0);
        let recip = reciprocal(&z);
        let product = multiply(&z, &recip);
        assert_abs_diff_eq!(product.real, 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(product.imag, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_is_nan() {
        let z = Complex::new(1.0, f64::NAN);
        assert!(z.is_nan());

        let w = Complex::new(f64::NAN, 1.0);
        assert!(w.is_nan());

        let u = Complex::new(1.0, 2.0);
        assert!(!u.is_nan());
    }

    #[test]
    fn test_is_infinite() {
        let z = Complex::new(f64::INFINITY, 2.0);
        assert!(z.is_infinite());

        let w = Complex::new(1.0, f64::INFINITY);
        assert!(w.is_infinite());

        let u = Complex::new(1.0, 2.0);
        assert!(!u.is_infinite());
    }

    #[test]
    fn test_is_finite() {
        let z = Complex::new(1.0, 2.0);
        assert!(z.is_finite());

        let w = Complex::new(f64::INFINITY, 2.0);
        assert!(!w.is_finite());

        let u = Complex::new(1.0, f64::NAN);
        assert!(!u.is_finite());
    }

    #[test]
    fn test_display() {
        let z1 = Complex::new(3.0, 4.0);
        assert_eq!(format!("{}", z1), "3 + 4i");

        let z2 = Complex::new(3.0, -4.0);
        assert_eq!(format!("{}", z2), "3 - 4i");

        let z3 = Complex::new(0.0, 1.0);
        assert_eq!(format!("{}", z3), "0 + 1i");
    }

    #[test]
    fn test_associativity() {
        // (a + b) + c = a + (b + c)
        let a = Complex::new(1.0, 2.0);
        let b = Complex::new(3.0, 4.0);
        let c = Complex::new(5.0, 6.0);

        let left = add(&add(&a, &b), &c);
        let right = add(&a, &add(&b, &c));

        assert_abs_diff_eq!(left.real, right.real, epsilon = 1e-10);
        assert_abs_diff_eq!(left.imag, right.imag, epsilon = 1e-10);
    }

    #[test]
    fn test_commutativity() {
        // a + b = b + a
        let a = Complex::new(1.0, 2.0);
        let b = Complex::new(3.0, 4.0);

        let left = add(&a, &b);
        let right = add(&b, &a);

        assert_eq!(left.real, right.real);
        assert_eq!(left.imag, right.imag);
    }

    #[test]
    fn test_distributivity() {
        // a * (b + c) = a * b + a * c
        let a = Complex::new(1.0, 2.0);
        let b = Complex::new(3.0, 4.0);
        let c = Complex::new(5.0, 6.0);

        let left = multiply(&a, &add(&b, &c));
        let right = add(&multiply(&a, &b), &multiply(&a, &c));

        assert_abs_diff_eq!(left.real, right.real, epsilon = 1e-10);
        assert_abs_diff_eq!(left.imag, right.imag, epsilon = 1e-10);
    }

    #[test]
    fn test_additive_identity() {
        // z + 0 = z
        let z = Complex::new(3.0, 4.0);
        let zero = Complex::zero();
        let result = add(&z, &zero);

        assert_eq!(result.real, z.real);
        assert_eq!(result.imag, z.imag);
    }

    #[test]
    fn test_multiplicative_identity() {
        // z * 1 = z
        let z = Complex::new(3.0, 4.0);
        let one = Complex::one();
        let result = multiply(&z, &one);

        assert_eq!(result.real, z.real);
        assert_eq!(result.imag, z.imag);
    }

    // Polar form tests
    #[test]
    fn test_magnitude() {
        let z = Complex::new(3.0, 4.0);
        let mag = magnitude(&z);
        assert_abs_diff_eq!(mag, 5.0, epsilon = 1e-10);
    }

    #[test]
    fn test_magnitude_zero() {
        let z: Complex<f64> = Complex::zero();
        let mag = magnitude(&z);
        assert_eq!(mag, 0.0);
    }

    #[test]
    fn test_magnitude_real() {
        let z = Complex::from_real(5.0);
        let mag = magnitude(&z);
        assert_abs_diff_eq!(mag, 5.0, epsilon = 1e-10);
    }

    #[test]
    fn test_magnitude_imaginary() {
        let z = Complex::from_imag(5.0);
        let mag = magnitude(&z);
        assert_abs_diff_eq!(mag, 5.0, epsilon = 1e-10);
    }

    #[test]
    fn test_magnitude_conjugate_property() {
        // |z̄| = |z|
        let z = Complex::new(3.0, 4.0);
        let conj = conjugate(&z);
        assert_abs_diff_eq!(magnitude(&z), magnitude(&conj), epsilon = 1e-10);
    }

    #[test]
    fn test_magnitude_multiplicative() {
        // |z₁ * z₂| = |z₁| * |z₂|
        let z1 = Complex::new(3.0, 4.0);
        let z2 = Complex::new(5.0, 12.0);
        let product = multiply(&z1, &z2);

        let mag_product = magnitude(&product);
        let product_of_mags = magnitude(&z1) * magnitude(&z2);

        assert_abs_diff_eq!(mag_product, product_of_mags, epsilon = 1e-10);
    }

    #[test]
    fn test_argument() {
        use std::f64::consts::PI;

        // 45 degrees
        let z = Complex::new(1.0, 1.0);
        let arg = argument(&z);
        assert_abs_diff_eq!(arg, PI / 4.0, epsilon = 1e-10);

        // 90 degrees
        let z = Complex::new(0.0, 1.0);
        let arg = argument(&z);
        assert_abs_diff_eq!(arg, PI / 2.0, epsilon = 1e-10);

        // 180 degrees
        let z = Complex::new(-1.0, 0.0);
        let arg = argument(&z);
        assert_abs_diff_eq!(arg.abs(), PI, epsilon = 1e-10);

        // -90 degrees
        let z = Complex::new(0.0, -1.0);
        let arg = argument(&z);
        assert_abs_diff_eq!(arg, -PI / 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_argument_conjugate_property() {
        // arg(z̄) = -arg(z)
        let z = Complex::new(3.0, 4.0);
        let conj = conjugate(&z);
        assert_abs_diff_eq!(argument(&conj), -argument(&z), epsilon = 1e-10);
    }

    #[test]
    fn test_from_polar() {
        use std::f64::consts::PI;

        let z = from_polar(5.0, PI / 2.0);
        assert_abs_diff_eq!(z.real, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(z.imag, 5.0, epsilon = 1e-10);
    }

    #[test]
    fn test_to_polar() {
        let z = Complex::new(3.0, 4.0);
        let (r, theta) = to_polar(&z);

        assert_abs_diff_eq!(r, 5.0, epsilon = 1e-10);
        assert_abs_diff_eq!(theta, 0.927295218, epsilon = 1e-9);
    }

    #[test]
    fn test_polar_round_trip() {
        // Convert to polar and back
        let z = Complex::new(3.0, 4.0);
        let (r, theta) = to_polar(&z);
        let z2 = from_polar(r, theta);

        assert_abs_diff_eq!(z2.real, z.real, epsilon = 1e-10);
        assert_abs_diff_eq!(z2.imag, z.imag, epsilon = 1e-10);
    }

    #[test]
    fn test_multiply_polar() {
        use std::f64::consts::PI;

        // Multiply two unit complex numbers at 45° each = 90°
        let result = multiply_polar(1.0, PI / 4.0, 1.0, PI / 4.0);

        assert_abs_diff_eq!(result.real, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_multiply_polar_magnitudes() {
        use std::f64::consts::PI;

        // 2 * e^(i π/3) × 3 * e^(i π/6) = 6 * e^(i π/2)
        let result = multiply_polar(2.0, PI / 3.0, 3.0, PI / 6.0);
        let (r, theta) = to_polar(&result);

        assert_abs_diff_eq!(r, 6.0, epsilon = 1e-10);
        assert_abs_diff_eq!(theta, PI / 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_magnitude_squared_via_conjugate() {
        // z * z̄ = |z|²
        let z = Complex::new(3.0, 4.0);
        let conj = conjugate(&z);
        let product = multiply(&z, &conj);

        let mag_squared = magnitude(&z).powi(2);

        assert_abs_diff_eq!(product.real, mag_squared, epsilon = 1e-10);
        assert_abs_diff_eq!(product.imag, 0.0, epsilon = 1e-10);
    }

    // Complex functions tests
    #[test]
    fn test_complex_exp() {
        use std::f64::consts::{E, PI};

        // e^0 = 1
        let z: Complex<f64> = Complex::zero();
        let result = complex_exp(&z);
        assert_abs_diff_eq!(result.real, 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-10);

        // e^1 = e
        let z = Complex::from_real(1.0);
        let result = complex_exp(&z);
        assert_abs_diff_eq!(result.real, E, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-10);

        // e^(iπ) = -1 (Euler's identity)
        let z = Complex::new(0.0, PI);
        let result = complex_exp(&z);
        assert_abs_diff_eq!(result.real, -1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_complex_exp_additivity() {
        // e^(z1+z2) = e^z1 * e^z2
        let z1 = Complex::new(1.0, 2.0);
        let z2 = Complex::new(3.0, 4.0);

        let sum = add(&z1, &z2);
        let exp_sum = complex_exp(&sum);

        let exp_z1 = complex_exp(&z1);
        let exp_z2 = complex_exp(&z2);
        let product = multiply(&exp_z1, &exp_z2);

        assert_abs_diff_eq!(exp_sum.real, product.real, epsilon = 1e-10);
        assert_abs_diff_eq!(exp_sum.imag, product.imag, epsilon = 1e-10);
    }

    #[test]
    fn test_complex_ln() {
        use std::f64::consts::E;

        // ln(1) = 0
        let z: Complex<f64> = Complex::one();
        let result = complex_ln(&z);
        assert_abs_diff_eq!(result.real, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-10);

        // ln(e) = 1
        let z = Complex::from_real(E);
        let result = complex_ln(&z);
        assert_abs_diff_eq!(result.real, 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_complex_ln_exp_inverse() {
        // ln(e^z) ≈ z (for principal branch)
        let z = Complex::new(1.0, 0.5); // Stay away from branch cut
        let exp_z = complex_exp(&z);
        let ln_exp_z = complex_ln(&exp_z);

        assert_abs_diff_eq!(ln_exp_z.real, z.real, epsilon = 1e-10);
        assert_abs_diff_eq!(ln_exp_z.imag, z.imag, epsilon = 1e-10);
    }

    #[test]
    fn test_complex_sqrt() {
        // √(-1) = i
        let z = Complex::new(-1.0, 0.0);
        let result = complex_sqrt(&z);
        assert_abs_diff_eq!(result.real, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 1.0, epsilon = 1e-10);

        // √0 = 0
        let z: Complex<f64> = Complex::zero();
        let result = complex_sqrt(&z);
        assert_abs_diff_eq!(result.real, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-10);

        // √4 = 2
        let z = Complex::from_real(4.0);
        let result = complex_sqrt(&z);
        assert_abs_diff_eq!(result.real, 2.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_complex_sqrt_property() {
        // (√z)² = z
        let z = Complex::new(3.0, 4.0);
        let sqrt_z = complex_sqrt(&z);
        let squared = multiply(&sqrt_z, &sqrt_z);

        assert_abs_diff_eq!(squared.real, z.real, epsilon = 1e-10);
        assert_abs_diff_eq!(squared.imag, z.imag, epsilon = 1e-10);
    }

    #[test]
    fn test_complex_pow() {
        // i^2 = -1
        let i = Complex::new(0.0, 1.0);
        let two = Complex::from_real(2.0);
        let result = complex_pow(&i, &two);
        assert_abs_diff_eq!(result.real, -1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-10);

        // 2^3 = 8
        let two = Complex::from_real(2.0);
        let three = Complex::from_real(3.0);
        let result = complex_pow(&two, &three);
        assert_abs_diff_eq!(result.real, 8.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-10);

        // 0^0 = 1 (by convention)
        let zero: Complex<f64> = Complex::zero();
        let result = complex_pow(&zero, &zero);
        assert_abs_diff_eq!(result.real, 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_complex_pow_i_to_i() {
        // i^i = e^(-π/2)
        let i = Complex::new(0.0, 1.0);
        let result = complex_pow(&i, &i);
        assert_abs_diff_eq!(result.real, 0.207879576, epsilon = 1e-9);
        assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-9);
    }

    #[test]
    fn test_complex_sin() {
        use std::f64::consts::PI;

        // sin(0) = 0
        let z: Complex<f64> = Complex::zero();
        let result = complex_sin(&z);
        assert_abs_diff_eq!(result.real, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-10);

        // sin(π/2) = 1
        let z = Complex::new(PI / 2.0, 0.0);
        let result = complex_sin(&z);
        assert_abs_diff_eq!(result.real, 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-10);

        // sin(π) ≈ 0
        let z = Complex::new(PI, 0.0);
        let result = complex_sin(&z);
        assert_abs_diff_eq!(result.real, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_complex_sin_odd() {
        // sin(-z) = -sin(z)
        let z = Complex::new(1.0, 2.0);
        let neg_z = negate(&z);

        let sin_z = complex_sin(&z);
        let sin_neg_z = complex_sin(&neg_z);

        assert_abs_diff_eq!(sin_neg_z.real, -sin_z.real, epsilon = 1e-10);
        assert_abs_diff_eq!(sin_neg_z.imag, -sin_z.imag, epsilon = 1e-10);
    }

    #[test]
    fn test_complex_cos() {
        use std::f64::consts::PI;

        // cos(0) = 1
        let z: Complex<f64> = Complex::zero();
        let result = complex_cos(&z);
        assert_abs_diff_eq!(result.real, 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-10);

        // cos(π/2) ≈ 0
        let z = Complex::new(PI / 2.0, 0.0);
        let result = complex_cos(&z);
        assert_abs_diff_eq!(result.real, 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-10);

        // cos(π) = -1
        let z = Complex::new(PI, 0.0);
        let result = complex_cos(&z);
        assert_abs_diff_eq!(result.real, -1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(result.imag, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_complex_cos_even() {
        // cos(-z) = cos(z)
        let z = Complex::new(1.0, 2.0);
        let neg_z = negate(&z);

        let cos_z = complex_cos(&z);
        let cos_neg_z = complex_cos(&neg_z);

        assert_abs_diff_eq!(cos_neg_z.real, cos_z.real, epsilon = 1e-10);
        assert_abs_diff_eq!(cos_neg_z.imag, cos_z.imag, epsilon = 1e-10);
    }

    #[test]
    fn test_sin_cos_identity() {
        // sin²(z) + cos²(z) = 1
        let z = Complex::new(1.0, 2.0);

        let sin_z = complex_sin(&z);
        let cos_z = complex_cos(&z);

        let sin_squared = multiply(&sin_z, &sin_z);
        let cos_squared = multiply(&cos_z, &cos_z);
        let sum = add(&sin_squared, &cos_squared);

        assert_abs_diff_eq!(sum.real, 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(sum.imag, 0.0, epsilon = 1e-10);
    }
}
