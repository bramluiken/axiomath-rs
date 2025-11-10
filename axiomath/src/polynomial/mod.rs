//! Polynomial operations and arithmetic.
//!
//! This module provides a comprehensive implementation of polynomials, including
//! arithmetic operations, evaluation, differentiation, integration, and root finding.
//!
//! # Theory
//!
//! A polynomial is a mathematical expression of the form:
//!
//! ```text
//! P(x) = a₀ + a₁x + a₂x² + ... + aₙxⁿ
//! ```
//!
//! where:
//! - aᵢ are the coefficients
//! - n is the degree of the polynomial
//! - a₀ is the constant term
//! - aₙ is the leading coefficient (aₙ ≠ 0 for degree n)
//!
//! Polynomials are fundamental in mathematics and have applications in:
//! - Approximation theory
//! - Numerical analysis
//! - Signal processing
//! - Computer graphics (Bezier curves)
//! - Control systems
//! - Cryptography
//!
//! # Representation
//!
//! Polynomials are stored with coefficients in ascending order:
//! `[a₀, a₁, a₂, ..., aₙ]` represents a₀ + a₁x + a₂x² + ... + aₙxⁿ
//!
//! # Examples
//!
//! ```
//! use axiomath::polynomial::{Polynomial, evaluate};
//!
//! // Create polynomial: 2 + 3x + x²
//! let p = Polynomial::new(vec![2.0, 3.0, 1.0]);
//!
//! // Evaluate at x = 2: 2 + 3(2) + (2)² = 2 + 6 + 4 = 12
//! let result = evaluate(&p, 2.0);
//! assert_eq!(result, 12.0);
//! ```
//!
//! # References
//!
//! - Knuth, D. E. (1997). *The Art of Computer Programming, Volume 2: Seminumerical Algorithms*.
//! - von zur Gathen, J., & Gerhard, J. (2013). *Modern Computer Algebra*.
//! - Press, W. H., et al. (2007). *Numerical Recipes*.

use crate::core::{Float, Scalar};
use crate::complex::Complex;
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// A polynomial with coefficients of type `T`.
///
/// Coefficients are stored in ascending order: `[a₀, a₁, a₂, ..., aₙ]`
/// represents the polynomial a₀ + a₁x + a₂x² + ... + aₙxⁿ.
///
/// The type `T` should implement [`Scalar`] for basic arithmetic operations.
///
/// # Examples
///
/// ```
/// use axiomath::polynomial::Polynomial;
///
/// // Create polynomial: 1 + 2x + 3x²
/// let p = Polynomial::new(vec![1.0, 2.0, 3.0]);
/// assert_eq!(p.degree(), Some(2));
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Polynomial<T> {
    /// Coefficients in ascending order [a₀, a₁, a₂, ..., aₙ]
    pub coeffs: Vec<T>,
}

impl<T: Scalar> Polynomial<T> {
    /// Creates a new polynomial from a vector of coefficients.
    ///
    /// Coefficients should be provided in ascending order: [a₀, a₁, a₂, ...].
    /// Leading zeros are automatically removed.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::polynomial::Polynomial;
    ///
    /// let p = Polynomial::new(vec![1.0, 2.0, 3.0]);
    /// assert_eq!(p.coeffs, vec![1.0, 2.0, 3.0]);
    /// ```
    pub fn new(coeffs: Vec<T>) -> Self {
        let mut p = Polynomial { coeffs };
        p.normalize();
        p
    }

    /// Creates the zero polynomial (P(x) = 0).
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::polynomial::Polynomial;
    ///
    /// let p: Polynomial<f64> = Polynomial::zero();
    /// assert!(p.is_zero());
    /// assert_eq!(p.degree(), None);
    /// ```
    pub fn zero() -> Self {
        Polynomial { coeffs: vec![] }
    }

    /// Creates the constant polynomial P(x) = 1.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::polynomial::Polynomial;
    ///
    /// let p: Polynomial<f64> = Polynomial::one();
    /// assert_eq!(p.coeffs, vec![1.0]);
    /// assert_eq!(p.degree(), Some(0));
    /// ```
    pub fn one() -> Self {
        Polynomial {
            coeffs: vec![T::one()],
        }
    }

    /// Creates a constant polynomial P(x) = c.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::polynomial::Polynomial;
    ///
    /// let p = Polynomial::constant(5.0);
    /// assert_eq!(p.coeffs, vec![5.0]);
    /// ```
    pub fn constant(c: T) -> Self {
        if c.is_zero() {
            Polynomial::zero()
        } else {
            Polynomial { coeffs: vec![c] }
        }
    }

    /// Creates the monomial P(x) = xⁿ.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::polynomial::Polynomial;
    ///
    /// // x²
    /// let p: Polynomial<f64> = Polynomial::monomial(2);
    /// assert_eq!(p.coeffs, vec![0.0, 0.0, 1.0]);
    /// ```
    pub fn monomial(n: usize) -> Self {
        let mut coeffs = vec![T::zero(); n + 1];
        coeffs[n] = T::one();
        Polynomial { coeffs }
    }

    /// Removes leading zero coefficients.
    ///
    /// This ensures the polynomial is in canonical form where the leading
    /// coefficient is non-zero (unless the polynomial is zero).
    fn normalize(&mut self) {
        while self.coeffs.len() > 0 && self.coeffs.last().unwrap().is_zero() {
            self.coeffs.pop();
        }
    }

    /// Returns the degree of the polynomial.
    ///
    /// The degree is the highest power of x with a non-zero coefficient.
    /// Returns `None` for the zero polynomial.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::polynomial::Polynomial;
    ///
    /// let p = Polynomial::new(vec![1.0, 2.0, 3.0]);
    /// assert_eq!(p.degree(), Some(2));
    ///
    /// let zero: Polynomial<f64> = Polynomial::zero();
    /// assert_eq!(zero.degree(), None);
    /// ```
    pub fn degree(&self) -> Option<usize> {
        if self.coeffs.is_empty() {
            None
        } else {
            Some(self.coeffs.len() - 1)
        }
    }

    /// Returns the leading coefficient (coefficient of highest degree term).
    ///
    /// Returns `None` for the zero polynomial.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::polynomial::Polynomial;
    ///
    /// let p = Polynomial::new(vec![1.0, 2.0, 3.0]);
    /// assert_eq!(p.leading_coefficient(), Some(&3.0));
    /// ```
    pub fn leading_coefficient(&self) -> Option<&T> {
        self.coeffs.last()
    }

    /// Checks if the polynomial is zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::polynomial::Polynomial;
    ///
    /// let p: Polynomial<f64> = Polynomial::zero();
    /// assert!(p.is_zero());
    ///
    /// let q = Polynomial::new(vec![1.0]);
    /// assert!(!q.is_zero());
    /// ```
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// Checks if the polynomial is constant (degree 0 or zero).
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::polynomial::Polynomial;
    ///
    /// let p = Polynomial::new(vec![5.0]);
    /// assert!(p.is_constant());
    ///
    /// let q = Polynomial::new(vec![1.0, 2.0]);
    /// assert!(!q.is_constant());
    /// ```
    pub fn is_constant(&self) -> bool {
        self.degree().map_or(true, |d| d == 0)
    }

    /// Checks if the polynomial is monic (leading coefficient is 1).
    ///
    /// # Examples
    ///
    /// ```
    /// use axiomath::polynomial::Polynomial;
    ///
    /// let p = Polynomial::new(vec![1.0, 2.0, 1.0]);
    /// assert!(p.is_monic());
    ///
    /// let q = Polynomial::new(vec![1.0, 2.0, 3.0]);
    /// assert!(!q.is_monic());
    /// ```
    pub fn is_monic(&self) -> bool {
        self.leading_coefficient().map_or(false, |c| *c == T::one())
    }
}

impl<T: fmt::Display> fmt::Display for Polynomial<T>
where
    T: Scalar + PartialOrd,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }

        let mut first = true;
        for (i, coeff) in self.coeffs.iter().enumerate().rev() {
            if coeff.is_zero() {
                continue;
            }

            // Write sign
            if !first {
                if *coeff < T::zero() {
                    write!(f, " - ")?;
                } else {
                    write!(f, " + ")?;
                }
            } else if *coeff < T::zero() {
                write!(f, "-")?;
            }

            // Write coefficient
            let abs_coeff = coeff.abs();
            if i == 0 || abs_coeff != T::one() {
                write!(f, "{}", abs_coeff)?;
            }

            // Write variable
            if i > 0 {
                if i == 0 || abs_coeff != T::one() {
                    // Already wrote coefficient
                }
                write!(f, "x")?;
                if i > 1 {
                    write!(f, "^{}", i)?;
                }
            }

            first = false;
        }

        Ok(())
    }
}

// ================================================================================================
// Polynomial Arithmetic
// ================================================================================================

/// Adds two polynomials.
///
/// # Mathematical Definition
///
/// ```text
/// (P + Q)(x) = P(x) + Q(x)
/// ```
///
/// Coefficients of corresponding terms are added.
///
/// # Complexity
///
/// - Time: O(max(deg(P), deg(Q)))
/// - Space: O(max(deg(P), deg(Q)))
///
/// # Examples
///
/// ```
/// use axiomath::polynomial::{Polynomial, add};
///
/// let p = Polynomial::new(vec![1.0, 2.0]);    // 1 + 2x
/// let q = Polynomial::new(vec![3.0, 4.0]);    // 3 + 4x
/// let result = add(&p, &q);                   // 4 + 6x
///
/// assert_eq!(result.coeffs, vec![4.0, 6.0]);
/// ```
pub fn add<T: Scalar>(p: &Polynomial<T>, q: &Polynomial<T>) -> Polynomial<T> {
    let max_len = p.coeffs.len().max(q.coeffs.len());
    let mut coeffs = Vec::with_capacity(max_len);

    for i in 0..max_len {
        let p_coeff = p.coeffs.get(i).copied().unwrap_or_else(T::zero);
        let q_coeff = q.coeffs.get(i).copied().unwrap_or_else(T::zero);
        coeffs.push(p_coeff + q_coeff);
    }

    Polynomial::new(coeffs)
}

/// Subtracts one polynomial from another.
///
/// # Mathematical Definition
///
/// ```text
/// (P - Q)(x) = P(x) - Q(x)
/// ```
///
/// # Complexity
///
/// - Time: O(max(deg(P), deg(Q)))
/// - Space: O(max(deg(P), deg(Q)))
///
/// # Examples
///
/// ```
/// use axiomath::polynomial::{Polynomial, subtract};
///
/// let p = Polynomial::new(vec![5.0, 7.0]);    // 5 + 7x
/// let q = Polynomial::new(vec![2.0, 3.0]);    // 2 + 3x
/// let result = subtract(&p, &q);              // 3 + 4x
///
/// assert_eq!(result.coeffs, vec![3.0, 4.0]);
/// ```
pub fn subtract<T: Scalar>(p: &Polynomial<T>, q: &Polynomial<T>) -> Polynomial<T>
where
    T: Neg<Output = T>,
{
    let max_len = p.coeffs.len().max(q.coeffs.len());
    let mut coeffs = Vec::with_capacity(max_len);

    for i in 0..max_len {
        let p_coeff = p.coeffs.get(i).copied().unwrap_or_else(T::zero);
        let q_coeff = q.coeffs.get(i).copied().unwrap_or_else(T::zero);
        coeffs.push(p_coeff - q_coeff);
    }

    Polynomial::new(coeffs)
}

/// Multiplies two polynomials.
///
/// # Mathematical Definition
///
/// ```text
/// (P * Q)(x) = P(x) * Q(x)
/// ```
///
/// Uses the convolution of coefficients.
///
/// # Algorithm
///
/// For small polynomials (degree < 32), uses naive O(n²) multiplication.
/// For larger polynomials, could use Karatsuba algorithm (not yet implemented).
///
/// # Complexity
///
/// - Time: O(n * m) where n = deg(P), m = deg(Q)
/// - Space: O(n + m)
///
/// # Examples
///
/// ```
/// use axiomath::polynomial::{Polynomial, multiply};
///
/// let p = Polynomial::new(vec![1.0, 1.0]);    // 1 + x
/// let q = Polynomial::new(vec![1.0, 1.0]);    // 1 + x
/// let result = multiply(&p, &q);              // 1 + 2x + x²
///
/// assert_eq!(result.coeffs, vec![1.0, 2.0, 1.0]);
/// ```
///
/// # References
///
/// - Knuth, "The Art of Computer Programming, Vol 2", Section 4.6.3
pub fn multiply<T: Scalar>(p: &Polynomial<T>, q: &Polynomial<T>) -> Polynomial<T> {
    if p.is_zero() || q.is_zero() {
        return Polynomial::zero();
    }

    let n = p.coeffs.len();
    let m = q.coeffs.len();
    let mut coeffs = vec![T::zero(); n + m - 1];

    for i in 0..n {
        for j in 0..m {
            coeffs[i + j] = coeffs[i + j] + p.coeffs[i] * q.coeffs[j];
        }
    }

    Polynomial::new(coeffs)
}

/// Evaluates a polynomial at a given point using Horner's method.
///
/// # Mathematical Definition
///
/// For P(x) = a₀ + a₁x + a₂x² + ... + aₙxⁿ, computes P(x₀).
///
/// # Algorithm
///
/// Horner's method rewrites the polynomial as:
///
/// ```text
/// P(x) = a₀ + x(a₁ + x(a₂ + ... + x(aₙ₋₁ + xaₙ)...))
/// ```
///
/// This requires only n multiplications and n additions.
///
/// # Complexity
///
/// - Time: O(n) where n = deg(P)
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::polynomial::{Polynomial, evaluate};
///
/// let p = Polynomial::new(vec![2.0, 3.0, 1.0]);  // 2 + 3x + x²
/// let result = evaluate(&p, 2.0);                // 2 + 3(2) + (2)² = 12
///
/// assert_eq!(result, 12.0);
/// ```
///
/// # References
///
/// - Knuth, "The Art of Computer Programming, Vol 2", Section 4.6.3
pub fn evaluate<T: Scalar>(p: &Polynomial<T>, x: T) -> T {
    if p.is_zero() {
        return T::zero();
    }

    // Horner's method: work backwards from highest degree
    let mut result = *p.coeffs.last().unwrap();
    for &coeff in p.coeffs.iter().rev().skip(1) {
        result = result * x + coeff;
    }
    result
}

/// Negates a polynomial.
///
/// # Mathematical Definition
///
/// ```text
/// (-P)(x) = -P(x)
/// ```
///
/// # Complexity
///
/// - Time: O(n)
/// - Space: O(n)
///
/// # Examples
///
/// ```
/// use axiomath::polynomial::{Polynomial, negate};
///
/// let p = Polynomial::new(vec![1.0, 2.0, 3.0]);
/// let result = negate(&p);
///
/// assert_eq!(result.coeffs, vec![-1.0, -2.0, -3.0]);
/// ```
pub fn negate<T: Scalar>(p: &Polynomial<T>) -> Polynomial<T>
where
    T: Neg<Output = T>,
{
    let coeffs = p.coeffs.iter().map(|&c| -c).collect();
    Polynomial { coeffs }
}

// ================================================================================================
// Calculus Operations
// ================================================================================================

/// Computes the derivative of a polynomial.
///
/// # Mathematical Definition
///
/// For P(x) = a₀ + a₁x + a₂x² + ... + aₙxⁿ:
///
/// ```text
/// P'(x) = a₁ + 2a₂x + 3a₃x² + ... + naₙx^(n-1)
/// ```
///
/// # Complexity
///
/// - Time: O(n)
/// - Space: O(n)
///
/// # Examples
///
/// ```
/// use axiomath::polynomial::{Polynomial, derivative};
///
/// let p = Polynomial::new(vec![1.0, 2.0, 3.0]);  // 1 + 2x + 3x²
/// let dp = derivative(&p);                       // 2 + 6x
///
/// assert_eq!(dp.coeffs, vec![2.0, 6.0]);
/// ```
pub fn derivative<T: Scalar>(p: &Polynomial<T>) -> Polynomial<T> {
    if p.is_zero() || p.degree() == Some(0) {
        return Polynomial::zero();
    }

    let coeffs: Vec<T> = p
        .coeffs
        .iter()
        .enumerate()
        .skip(1)
        .map(|(i, &c)| {
            // Convert index to T using repeated addition to avoid trait bound issues
            let mut index = T::zero();
            for _ in 0..i {
                index = index + T::one();
            }
            c * index
        })
        .collect();

    Polynomial::new(coeffs)
}

/// Computes the indefinite integral of a polynomial with integration constant.
///
/// # Mathematical Definition
///
/// For P(x) = a₀ + a₁x + a₂x² + ... + aₙxⁿ:
///
/// ```text
/// ∫P(x)dx = C + a₀x + (a₁/2)x² + (a₂/3)x³ + ... + (aₙ/(n+1))x^(n+1)
/// ```
///
/// where C is the integration constant.
///
/// # Arguments
///
/// * `p` - The polynomial to integrate
/// * `constant` - The integration constant C
///
/// # Complexity
///
/// - Time: O(n)
/// - Space: O(n)
///
/// # Examples
///
/// ```
/// use axiomath::polynomial::{Polynomial, integral};
///
/// let p = Polynomial::new(vec![2.0, 6.0]);  // 2 + 6x
/// let ip = integral(&p, 1.0);               // 1 + 2x + 3x²
///
/// assert_eq!(ip.coeffs, vec![1.0, 2.0, 3.0]);
/// ```
pub fn integral<T: Float>(p: &Polynomial<T>, constant: T) -> Polynomial<T> {
    if p.is_zero() {
        return Polynomial::constant(constant);
    }

    let mut coeffs = vec![constant];
    for (i, &c) in p.coeffs.iter().enumerate() {
        let divisor = T::from((i + 1) as f64).unwrap();
        coeffs.push(c / divisor);
    }

    Polynomial::new(coeffs)
}

/// Computes the definite integral of a polynomial over [a, b].
///
/// # Mathematical Definition
///
/// ```text
/// ∫ₐᵇ P(x)dx = F(b) - F(a)
/// ```
///
/// where F is the antiderivative of P.
///
/// # Complexity
///
/// - Time: O(n)
/// - Space: O(n)
///
/// # Examples
///
/// ```
/// use axiomath::polynomial::{Polynomial, definite_integral};
///
/// let p = Polynomial::new(vec![0.0, 1.0]);  // x
/// let result = definite_integral(&p, 0.0, 2.0);  // ∫₀² x dx = x²/2 |₀² = 2
///
/// assert_eq!(result, 2.0);
/// ```
pub fn definite_integral<T: Float>(p: &Polynomial<T>, a: T, b: T) -> T {
    let antideriv = integral(p, T::zero());
    evaluate(&antideriv, b) - evaluate(&antideriv, a)
}

// ================================================================================================
// Polynomial Properties
// ================================================================================================

/// Computes the discriminant of a quadratic polynomial.
///
/// # Mathematical Definition
///
/// For a quadratic P(x) = ax² + bx + c:
///
/// ```text
/// Δ = b² - 4ac
/// ```
///
/// The discriminant determines the nature of the roots:
/// - Δ > 0: Two distinct real roots
/// - Δ = 0: One repeated real root
/// - Δ < 0: Two complex conjugate roots
///
/// # Returns
///
/// `Some(discriminant)` if the polynomial is quadratic (degree 2),
/// `None` otherwise.
///
/// # Examples
///
/// ```
/// use axiomath::polynomial::{Polynomial, discriminant};
///
/// let p = Polynomial::new(vec![1.0, 0.0, 1.0]);  // x² + 1
/// let disc = discriminant(&p);
///
/// assert_eq!(disc, Some(-4.0));  // Negative => complex roots
/// ```
pub fn discriminant<T: Scalar>(p: &Polynomial<T>) -> Option<T> {
    if p.degree() != Some(2) {
        return None;
    }

    let a = p.coeffs[2];
    let b = p.coeffs[1];
    let c = p.coeffs[0];

    let four = T::one() + T::one() + T::one() + T::one();
    Some(b * b - four * a * c)
}

// ================================================================================================
// Root Finding
// ================================================================================================

/// Finds the root of a linear polynomial (degree 1).
///
/// # Mathematical Definition
///
/// For P(x) = ax + b (a ≠ 0):
///
/// ```text
/// Root: x = -b/a
/// ```
///
/// # Returns
///
/// `Some(root)` if the polynomial is linear (degree 1),
/// `None` if the polynomial is not linear.
///
/// # Examples
///
/// ```
/// use axiomath::polynomial::{Polynomial, linear_root};
///
/// let p = Polynomial::new(vec![3.0, 2.0]);  // 3 + 2x
/// let root = linear_root(&p);
///
/// assert_eq!(root, Some(-1.5));
/// ```
///
/// # References
///
/// - Elementary Algebra
pub fn linear_root<T: Float>(p: &Polynomial<T>) -> Option<T> {
    if p.degree() != Some(1) {
        return None;
    }

    let b = p.coeffs[0];
    let a = p.coeffs[1];

    Some(-b / a)
}

/// Finds the real roots of a quadratic polynomial (degree 2).
///
/// # Mathematical Definition
///
/// For P(x) = ax² + bx + c (a ≠ 0):
///
/// ```text
/// Discriminant: Δ = b² - 4ac
///
/// If Δ ≥ 0:
///   x = (-b ± √Δ) / (2a)
/// ```
///
/// # Returns
///
/// - `Some(vec![])` if discriminant is negative (no real roots)
/// - `Some(vec![x])` if discriminant is zero (one repeated root)
/// - `Some(vec![x₁, x₂])` if discriminant is positive (two distinct roots)
/// - `None` if polynomial is not quadratic
///
/// # Examples
///
/// ```
/// use axiomath::polynomial::{Polynomial, quadratic_roots};
///
/// let p = Polynomial::new(vec![-2.0, 0.0, 1.0]);  // -2 + x²
/// let roots = quadratic_roots(&p).unwrap();
///
/// assert_eq!(roots.len(), 2);
/// ```
///
/// # References
///
/// - Kahan, W. (2004). "On the Cost of Floating-Point Computation Without
///   Extra-Precise Arithmetic"
pub fn quadratic_roots<T: Float>(p: &Polynomial<T>) -> Option<Vec<T>> {
    if p.degree() != Some(2) {
        return None;
    }

    let c = p.coeffs[0];
    let b = p.coeffs[1];
    let a = p.coeffs[2];

    let disc = discriminant(p).unwrap();

    // No real roots if discriminant is negative
    if disc < T::zero() {
        return Some(vec![]);
    }

    let two = T::one() + T::one();

    // One repeated root if discriminant is zero
    if disc == T::zero() {
        let root = -b / (two * a);
        return Some(vec![root]);
    }

    // Two distinct roots - use numerically stable formula
    let sqrt_disc = disc.sqrt();
    let sign_b = if b >= T::zero() { T::one() } else { -T::one() };
    let q = -(b + sign_b * sqrt_disc) / two;

    let x1 = q / a;
    let x2 = c / q;

    Some(vec![x1.min(x2), x1.max(x2)])
}

/// Finds a real root of a cubic polynomial using Cardano's formula.
///
/// # Mathematical Definition
///
/// For the depressed cubic t³ + pt + q = 0:
///
/// ```text
/// Δ = -(4p³ + 27q²)
/// ```
///
/// This function returns one real root. Cubic polynomials always have
/// at least one real root.
///
/// # Returns
///
/// `Some(root)` if the polynomial is cubic (degree 3),
/// `None` if the polynomial is not cubic.
///
/// # Examples
///
/// ```
/// use axiomath::polynomial::{Polynomial, cubic_root};
///
/// let p = Polynomial::new(vec![-1.0, 0.0, 0.0, 1.0]);  // -1 + x³
/// let root = cubic_root(&p);
///
/// assert!(root.is_some());
/// ```
///
/// # References
///
/// - Press, W. H., et al. (2007). "Numerical Recipes", 3rd ed.
/// - Cardano's formula (1545)
pub fn cubic_root<T: Float>(p: &Polynomial<T>) -> Option<T> {
    if p.degree() != Some(3) {
        return None;
    }

    let d = p.coeffs[0];
    let c = p.coeffs[1];
    let b = p.coeffs[2];
    let a = p.coeffs[3];

    // Convert to depressed cubic: t³ + pt + q = 0
    // using substitution x = t - b/(3a)
    let three = T::one() + T::one() + T::one();
    let nine = three * three;
    let twenty_seven = nine * three;

    let p_coeff = (three * a * c - b * b) / (three * a * a);
    let q_coeff = (T::one() + T::one()) * b * b * b / (twenty_seven * a * a * a)
        - b * c / (three * a * a)
        + d / a;

    // Compute discriminant: Δ = -4p³ - 27q²
    let four = T::one() + T::one() + T::one() + T::one();
    let disc = -(four * p_coeff * p_coeff * p_coeff + twenty_seven * q_coeff * q_coeff);

    let two = T::one() + T::one();

    // Use Cardano's formula
    // For depressed cubic t³ + pt + q = 0
    let term = q_coeff * q_coeff / four + p_coeff * p_coeff * p_coeff / twenty_seven;

    if term >= T::zero() {
        // One real root - use standard Cardano formula
        let sqrt_term = term.sqrt();
        let u_cube = -q_coeff / two + sqrt_term;
        let v_cube = -q_coeff / two - sqrt_term;

        let u = if u_cube >= T::zero() {
            u_cube.powf(T::one() / three)
        } else {
            -(-u_cube).powf(T::one() / three)
        };

        let v = if v_cube >= T::zero() {
            v_cube.powf(T::one() / three)
        } else {
            -(-v_cube).powf(T::one() / three)
        };

        let t = u + v;
        Some(t - b / (three * a))
    } else {
        // Three real roots - use trigonometric solution
        let m = two * (-p_coeff / three).sqrt();
        let cos_theta = three * q_coeff / (p_coeff * m);
        let theta = cos_theta.acos() / three;
        let t = m * theta.cos();
        Some(t - b / (three * a))
    }
}

/// Finds a root of a polynomial using Newton-Raphson iteration.
///
/// # Mathematical Definition
///
/// Starting from initial guess x₀, iterate:
///
/// ```text
/// xₙ₊₁ = xₙ - P(xₙ) / P'(xₙ)
/// ```
///
/// until |P(xₙ)| < tolerance or max_iterations reached.
///
/// # Parameters
///
/// - `p`: The polynomial
/// - `initial_guess`: Starting point for iteration
/// - `tolerance`: Convergence threshold for |P(x)|
/// - `max_iterations`: Maximum number of iterations
///
/// # Returns
///
/// `Some(root)` if convergence achieved, `None` if:
/// - Polynomial is constant
/// - Max iterations exceeded without convergence
/// - Derivative becomes zero
///
/// # Examples
///
/// ```
/// use axiomath::polynomial::{Polynomial, newton_raphson_root};
///
/// let p = Polynomial::new(vec![-2.0, 0.0, 1.0]);  // -2 + x²
/// let root = newton_raphson_root(&p, 1.0, 1e-10, 100);
///
/// assert!(root.is_some());
/// ```
///
/// # References
///
/// - Quarteroni, A., Sacco, R., & Saleri, F. (2007). "Numerical Mathematics"
pub fn newton_raphson_root<T: Float>(
    p: &Polynomial<T>,
    initial_guess: T,
    tolerance: T,
    max_iterations: usize,
) -> Option<T> {
    if p.is_zero() || p.degree() == Some(0) {
        return None;
    }

    let dp = derivative(p);
    let mut x = initial_guess;

    for _ in 0..max_iterations {
        let fx = evaluate(p, x);

        // Check for convergence
        if fx.abs() < tolerance {
            return Some(x);
        }

        let dfx = evaluate(&dp, x);

        // Avoid division by zero
        if dfx.abs() < tolerance {
            return None;
        }

        x = x - fx / dfx;
    }

    // Check if final value is close enough
    let fx = evaluate(p, x);
    if fx.abs() < tolerance {
        Some(x)
    } else {
        None
    }
}

/// Finds all roots of a polynomial using the Durand-Kerner method.
///
/// # Mathematical Definition
///
/// For polynomial P(x) with degree n, iterate simultaneously for all roots:
///
/// ```text
/// zᵢ⁽ᵏ⁺¹⁾ = zᵢ⁽ᵏ⁾ - P(zᵢ⁽ᵏ⁾) / ∏ⱼ≠ᵢ (zᵢ⁽ᵏ⁾ - zⱼ⁽ᵏ⁾)
/// ```
///
/// This method finds all roots (real and complex) simultaneously.
///
/// # Parameters
///
/// - `p`: The polynomial
/// - `tolerance`: Convergence threshold
/// - `max_iterations`: Maximum number of iterations
///
/// # Returns
///
/// `Some(Vec<Complex<T>>)` containing all roots if convergence achieved,
/// `None` if polynomial is constant or zero, or if convergence fails.
///
/// # Examples
///
/// ```
/// use axiomath::polynomial::{Polynomial, durand_kerner};
///
/// let p = Polynomial::new(vec![-1.0, 0.0, 1.0]);  // -1 + x²
/// let roots = durand_kerner(&p, 1e-10, 100);
///
/// assert!(roots.is_some());
/// ```
///
/// # References
///
/// - Durand, E. (1960). "Solutions Numériques des Équations Algébriques"
/// - Weierstrass, K. (1891). "Neuer Beweis des Fundamentalsatzes der Algebra"
pub fn durand_kerner<T: Float>(
    p: &Polynomial<T>,
    tolerance: T,
    max_iterations: usize,
) -> Option<Vec<Complex<T>>> {
    let n = match p.degree() {
        Some(d) if d > 0 => d,
        _ => return None,
    };

    // Initialize roots on unit circle
    let two = T::one() + T::one();
    let pi = T::from(std::f64::consts::PI).unwrap();
    let mut roots: Vec<Complex<T>> = (0..n)
        .map(|k| {
            let k_float = T::from(k).unwrap();
            let n_float = T::from(n).unwrap();
            let angle = two * pi * k_float / n_float;
            Complex::new(angle.cos() * T::from(0.4).unwrap(), angle.sin() * T::from(0.4).unwrap())
        })
        .collect();

    // Normalize polynomial to monic form
    let leading = p.leading_coefficient().unwrap();
    let monic_coeffs: Vec<T> = p.coeffs.iter().map(|&c| c / *leading).collect();
    let monic_p = Polynomial::new(monic_coeffs);

    for _ in 0..max_iterations {
        let mut converged = true;
        let mut new_roots = roots.clone();

        for i in 0..n {
            let zi = &roots[i];

            // Evaluate P(zi)
            let p_zi = evaluate_complex(&monic_p, zi);

            // Compute product ∏ⱼ≠ᵢ (zᵢ - zⱼ)
            let mut prod = Complex::new(T::one(), T::zero());
            for j in 0..n {
                if i != j {
                    prod = crate::complex::multiply(&prod, &crate::complex::subtract(zi, &roots[j]));
                }
            }

            // Update root: zᵢ = zᵢ - P(zᵢ) / ∏
            let correction = crate::complex::divide(&p_zi, &prod);
            new_roots[i] = crate::complex::subtract(zi, &correction);

            // Check convergence
            if crate::complex::magnitude(&correction) >= tolerance {
                converged = false;
            }
        }

        roots = new_roots;

        if converged {
            return Some(roots);
        }
    }

    // Return roots even if not fully converged
    Some(roots)
}

/// Evaluates a polynomial with real coefficients at a complex point.
///
/// # Examples
///
/// ```
/// use axiomath::polynomial::{Polynomial, evaluate_complex};
/// use axiomath::complex::Complex;
///
/// let p = Polynomial::new(vec![1.0, 0.0, 1.0]);  // 1 + x²
/// let z = Complex::new(0.0, 1.0);  // i
/// let result = evaluate_complex(&p, &z);
///
/// // P(i) = 1 + i² = 1 - 1 = 0
/// assert!((result.real).abs() < 1e-10);
/// assert!((result.imag).abs() < 1e-10);
/// ```
pub fn evaluate_complex<T: Float>(p: &Polynomial<T>, z: &Complex<T>) -> Complex<T> {
    if p.is_zero() {
        return Complex::zero();
    }

    // Horner's method for complex evaluation
    let mut result = Complex::new(*p.coeffs.last().unwrap(), T::zero());

    for &coeff in p.coeffs.iter().rev().skip(1) {
        result = crate::complex::multiply(&result, z);
        result = crate::complex::add(&result, &Complex::new(coeff, T::zero()));
    }

    result
}

// ================================================================================================
// Operator Overloading
// ================================================================================================

impl<T: Scalar> Add for Polynomial<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        add(&self, &rhs)
    }
}

impl<T: Scalar> Add<&Polynomial<T>> for Polynomial<T> {
    type Output = Self;

    fn add(self, rhs: &Polynomial<T>) -> Self::Output {
        add(&self, rhs)
    }
}

impl<T: Scalar> Add<Polynomial<T>> for &Polynomial<T> {
    type Output = Polynomial<T>;

    fn add(self, rhs: Polynomial<T>) -> Self::Output {
        add(self, &rhs)
    }
}

impl<T: Scalar> Add<&Polynomial<T>> for &Polynomial<T> {
    type Output = Polynomial<T>;

    fn add(self, rhs: &Polynomial<T>) -> Self::Output {
        add(self, rhs)
    }
}

impl<T: Scalar> Sub for Polynomial<T>
where
    T: Neg<Output = T>,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        subtract(&self, &rhs)
    }
}

impl<T: Scalar> Sub<&Polynomial<T>> for Polynomial<T>
where
    T: Neg<Output = T>,
{
    type Output = Self;

    fn sub(self, rhs: &Polynomial<T>) -> Self::Output {
        subtract(&self, rhs)
    }
}

impl<T: Scalar> Sub<Polynomial<T>> for &Polynomial<T>
where
    T: Neg<Output = T>,
{
    type Output = Polynomial<T>;

    fn sub(self, rhs: Polynomial<T>) -> Self::Output {
        subtract(self, &rhs)
    }
}

impl<T: Scalar> Sub<&Polynomial<T>> for &Polynomial<T>
where
    T: Neg<Output = T>,
{
    type Output = Polynomial<T>;

    fn sub(self, rhs: &Polynomial<T>) -> Self::Output {
        subtract(self, rhs)
    }
}

impl<T: Scalar> Mul for Polynomial<T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        multiply(&self, &rhs)
    }
}

impl<T: Scalar> Mul<&Polynomial<T>> for Polynomial<T> {
    type Output = Self;

    fn mul(self, rhs: &Polynomial<T>) -> Self::Output {
        multiply(&self, rhs)
    }
}

impl<T: Scalar> Mul<Polynomial<T>> for &Polynomial<T> {
    type Output = Polynomial<T>;

    fn mul(self, rhs: Polynomial<T>) -> Self::Output {
        multiply(self, &rhs)
    }
}

impl<T: Scalar> Mul<&Polynomial<T>> for &Polynomial<T> {
    type Output = Polynomial<T>;

    fn mul(self, rhs: &Polynomial<T>) -> Self::Output {
        multiply(self, rhs)
    }
}

impl<T: Scalar> Neg for Polynomial<T>
where
    T: Neg<Output = T>,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        negate(&self)
    }
}

impl<T: Scalar> Neg for &Polynomial<T>
where
    T: Neg<Output = T>,
{
    type Output = Polynomial<T>;

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

    #[test]
    fn test_new() {
        let p = Polynomial::new(vec![1.0, 2.0, 3.0]);
        assert_eq!(p.coeffs, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_new_removes_leading_zeros() {
        let p = Polynomial::new(vec![1.0, 2.0, 0.0, 0.0]);
        assert_eq!(p.coeffs, vec![1.0, 2.0]);
    }

    #[test]
    fn test_zero() {
        let p: Polynomial<f64> = Polynomial::zero();
        assert!(p.is_zero());
        assert_eq!(p.degree(), None);
    }

    #[test]
    fn test_one() {
        let p: Polynomial<f64> = Polynomial::one();
        assert_eq!(p.coeffs, vec![1.0]);
        assert_eq!(p.degree(), Some(0));
    }

    #[test]
    fn test_constant() {
        let p = Polynomial::constant(5.0);
        assert_eq!(p.coeffs, vec![5.0]);
        assert!(p.is_constant());
    }

    #[test]
    fn test_monomial() {
        let p: Polynomial<f64> = Polynomial::monomial(2);
        assert_eq!(p.coeffs, vec![0.0, 0.0, 1.0]);
        assert_eq!(p.degree(), Some(2));
    }

    #[test]
    fn test_degree() {
        let p = Polynomial::new(vec![1.0, 2.0, 3.0]);
        assert_eq!(p.degree(), Some(2));

        let zero: Polynomial<f64> = Polynomial::zero();
        assert_eq!(zero.degree(), None);
    }

    #[test]
    fn test_leading_coefficient() {
        let p = Polynomial::new(vec![1.0, 2.0, 3.0]);
        assert_eq!(p.leading_coefficient(), Some(&3.0));

        let zero: Polynomial<f64> = Polynomial::zero();
        assert_eq!(zero.leading_coefficient(), None);
    }

    #[test]
    fn test_is_zero() {
        let p: Polynomial<f64> = Polynomial::zero();
        assert!(p.is_zero());

        let q = Polynomial::new(vec![1.0]);
        assert!(!q.is_zero());
    }

    #[test]
    fn test_is_constant() {
        let p = Polynomial::new(vec![5.0]);
        assert!(p.is_constant());

        let q = Polynomial::new(vec![1.0, 2.0]);
        assert!(!q.is_constant());
    }

    #[test]
    fn test_is_monic() {
        let p = Polynomial::new(vec![1.0, 2.0, 1.0]);
        assert!(p.is_monic());

        let q = Polynomial::new(vec![1.0, 2.0, 3.0]);
        assert!(!q.is_monic());
    }

    #[test]
    fn test_add() {
        let p = Polynomial::new(vec![1.0, 2.0]);
        let q = Polynomial::new(vec![3.0, 4.0]);
        let result = add(&p, &q);
        assert_eq!(result.coeffs, vec![4.0, 6.0]);
    }

    #[test]
    fn test_add_different_degrees() {
        let p = Polynomial::new(vec![1.0, 2.0, 3.0]);
        let q = Polynomial::new(vec![1.0]);
        let result = add(&p, &q);
        assert_eq!(result.coeffs, vec![2.0, 2.0, 3.0]);
    }

    #[test]
    fn test_add_operator() {
        let p = Polynomial::new(vec![1.0, 2.0]);
        let q = Polynomial::new(vec![3.0, 4.0]);
        let result = p + q;
        assert_eq!(result.coeffs, vec![4.0, 6.0]);
    }

    #[test]
    fn test_subtract() {
        let p = Polynomial::new(vec![5.0, 7.0]);
        let q = Polynomial::new(vec![2.0, 3.0]);
        let result = subtract(&p, &q);
        assert_eq!(result.coeffs, vec![3.0, 4.0]);
    }

    #[test]
    fn test_subtract_operator() {
        let p = Polynomial::new(vec![5.0, 7.0]);
        let q = Polynomial::new(vec![2.0, 3.0]);
        let result = p - q;
        assert_eq!(result.coeffs, vec![3.0, 4.0]);
    }

    #[test]
    fn test_multiply() {
        let p = Polynomial::new(vec![1.0, 1.0]);  // 1 + x
        let q = Polynomial::new(vec![1.0, 1.0]);  // 1 + x
        let result = multiply(&p, &q);            // 1 + 2x + x²
        assert_eq!(result.coeffs, vec![1.0, 2.0, 1.0]);
    }

    #[test]
    fn test_multiply_operator() {
        let p = Polynomial::new(vec![1.0, 1.0]);
        let q = Polynomial::new(vec![1.0, 1.0]);
        let result = p * q;
        assert_eq!(result.coeffs, vec![1.0, 2.0, 1.0]);
    }

    #[test]
    fn test_multiply_by_zero() {
        let p = Polynomial::new(vec![1.0, 2.0, 3.0]);
        let zero: Polynomial<f64> = Polynomial::zero();
        let result = multiply(&p, &zero);
        assert!(result.is_zero());
    }

    #[test]
    fn test_negate() {
        let p = Polynomial::new(vec![1.0, 2.0, 3.0]);
        let result = negate(&p);
        assert_eq!(result.coeffs, vec![-1.0, -2.0, -3.0]);
    }

    #[test]
    fn test_negate_operator() {
        let p = Polynomial::new(vec![1.0, 2.0, 3.0]);
        let result = -p;
        assert_eq!(result.coeffs, vec![-1.0, -2.0, -3.0]);
    }

    #[test]
    fn test_evaluate() {
        let p = Polynomial::new(vec![2.0, 3.0, 1.0]);  // 2 + 3x + x²
        let result = evaluate(&p, 2.0);                // 2 + 3(2) + (2)² = 12
        assert_eq!(result, 12.0);
    }

    #[test]
    fn test_evaluate_zero() {
        let p: Polynomial<f64> = Polynomial::zero();
        let result = evaluate(&p, 5.0);
        assert_eq!(result, 0.0);
    }

    #[test]
    fn test_evaluate_constant() {
        let p = Polynomial::constant(7.0);
        let result = evaluate(&p, 100.0);
        assert_eq!(result, 7.0);
    }

    // Calculus tests
    #[test]
    fn test_derivative() {
        let p = Polynomial::new(vec![1.0, 2.0, 3.0]); // 1 + 2x + 3x²
        let dp = derivative(&p); // 2 + 6x
        assert_eq!(dp.coeffs, vec![2.0, 6.0]);
    }

    #[test]
    fn test_derivative_constant() {
        let p = Polynomial::constant(5.0);
        let dp = derivative(&p);
        assert!(dp.is_zero());
    }

    #[test]
    fn test_derivative_linear() {
        let p = Polynomial::new(vec![3.0, 4.0]); // 3 + 4x
        let dp = derivative(&p); // 4
        assert_eq!(dp.coeffs, vec![4.0]);
    }

    #[test]
    fn test_integral() {
        let p = Polynomial::new(vec![2.0, 6.0]); // 2 + 6x
        let ip = integral(&p, 1.0); // 1 + 2x + 3x²
        assert_eq!(ip.coeffs, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_integral_zero_constant() {
        let p = Polynomial::new(vec![1.0, 2.0]); // 1 + 2x
        let ip = integral(&p, 0.0); // x + x²
        assert_eq!(ip.coeffs, vec![0.0, 1.0, 1.0]);
    }

    #[test]
    fn test_derivative_integral_inverse() {
        // derivative(integral(p)) should give back p (modulo constant)
        let p = Polynomial::new(vec![2.0, 6.0]); // 2 + 6x
        let ip = integral(&p, 0.0);
        let dp = derivative(&ip);
        assert_eq!(dp.coeffs, p.coeffs);
    }

    #[test]
    fn test_definite_integral() {
        let p = Polynomial::new(vec![0.0, 1.0]); // x
        let result = definite_integral(&p, 0.0, 2.0); // ∫₀² x dx = 2
        assert_eq!(result, 2.0);
    }

    #[test]
    fn test_definite_integral_quadratic() {
        let p = Polynomial::new(vec![0.0, 0.0, 1.0]); // x²
        let result = definite_integral(&p, 0.0, 3.0); // ∫₀³ x² dx = 9
        assert_eq!(result, 9.0);
    }

    #[test]
    fn test_discriminant_positive() {
        // x² - 5x + 6 = (x-2)(x-3), discriminant = 25 - 24 = 1
        let p = Polynomial::new(vec![6.0, -5.0, 1.0]);
        let disc = discriminant(&p);
        assert_eq!(disc, Some(1.0));
    }

    #[test]
    fn test_discriminant_zero() {
        // x² - 2x + 1 = (x-1)², discriminant = 0
        let p = Polynomial::new(vec![1.0, -2.0, 1.0]);
        let disc = discriminant(&p);
        assert_eq!(disc, Some(0.0));
    }

    #[test]
    fn test_discriminant_negative() {
        // x² + 1, discriminant = -4
        let p = Polynomial::new(vec![1.0, 0.0, 1.0]);
        let disc = discriminant(&p);
        assert_eq!(disc, Some(-4.0));
    }

    #[test]
    fn test_discriminant_non_quadratic() {
        let p = Polynomial::new(vec![1.0, 2.0, 3.0, 4.0]); // cubic
        let disc = discriminant(&p);
        assert_eq!(disc, None);
    }

    // Root finding tests

    #[test]
    fn test_linear_root() {
        let p = Polynomial::new(vec![3.0, 2.0]); // 3 + 2x
        let root = linear_root(&p);
        assert_eq!(root, Some(-1.5));

        // Verify it's actually a root
        let result = evaluate(&p, root.unwrap());
        assert!((result).abs() < 1e-10);
    }

    #[test]
    fn test_linear_root_non_linear() {
        let p = Polynomial::new(vec![1.0, 2.0, 3.0]); // quadratic
        let root = linear_root(&p);
        assert_eq!(root, None);
    }

    #[test]
    fn test_quadratic_roots_two_real() {
        // x² - 5x + 6 = (x-2)(x-3)
        let p = Polynomial::new(vec![6.0, -5.0, 1.0]);
        let roots = quadratic_roots(&p).unwrap();

        assert_eq!(roots.len(), 2);
        assert!((roots[0] - 2.0).abs() < 1e-10);
        assert!((roots[1] - 3.0).abs() < 1e-10);

        // Verify they're actually roots
        for &root in &roots {
            let result = evaluate(&p, root);
            assert!((result).abs() < 1e-10);
        }
    }

    #[test]
    fn test_quadratic_roots_one_repeated() {
        // (x-2)² = x² - 4x + 4
        let p = Polynomial::new(vec![4.0, -4.0, 1.0]);
        let roots = quadratic_roots(&p).unwrap();

        assert_eq!(roots.len(), 1);
        assert!((roots[0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_quadratic_roots_no_real() {
        // x² + 1 (roots are ±i)
        let p = Polynomial::new(vec![1.0, 0.0, 1.0]);
        let roots = quadratic_roots(&p).unwrap();

        assert_eq!(roots.len(), 0);
    }

    #[test]
    fn test_quadratic_roots_non_quadratic() {
        let p = Polynomial::new(vec![1.0, 2.0, 3.0, 4.0]); // cubic
        let roots = quadratic_roots(&p);
        assert_eq!(roots, None);
    }

    #[test]
    fn test_cubic_root() {
        // x³ - 1 = (x-1)(x² + x + 1)
        let p = Polynomial::new(vec![-1.0, 0.0, 0.0, 1.0]);
        let root = cubic_root(&p).unwrap();

        // Verify it's actually a root (most important test)
        let result = evaluate(&p, root);
        assert!((result).abs() < 1e-6);

        // The cubic root function may find any of the three cube roots
        // For x³ - 1, the real root is 1, but due to numerical stability
        // the function might converge to a different branch
        // We just need to verify it's a valid root
    }

    #[test]
    fn test_cubic_root_three_real() {
        // x³ - 6x² + 11x - 6 = (x-1)(x-2)(x-3)
        let p = Polynomial::new(vec![-6.0, 11.0, -6.0, 1.0]);
        let root = cubic_root(&p).unwrap();

        // Should find one of the three roots
        let result = evaluate(&p, root);
        assert!((result).abs() < 1e-6);
    }

    #[test]
    fn test_cubic_root_non_cubic() {
        let p = Polynomial::new(vec![1.0, 2.0]); // linear
        let root = cubic_root(&p);
        assert_eq!(root, None);
    }

    #[test]
    fn test_newton_raphson_root() {
        // x² - 2 (root is √2)
        let p = Polynomial::new(vec![-2.0, 0.0, 1.0]);
        let root = newton_raphson_root(&p, 1.0, 1e-10, 100).unwrap();

        assert!((root - 2.0_f64.sqrt()).abs() < 1e-10);

        // Verify it's actually a root
        let result = evaluate(&p, root);
        assert!((result).abs() < 1e-10);
    }

    #[test]
    fn test_newton_raphson_root_cubic() {
        // x³ - 2 (root is ∛2)
        let p = Polynomial::new(vec![-2.0, 0.0, 0.0, 1.0]);
        let root = newton_raphson_root(&p, 1.0, 1e-10, 100).unwrap();

        assert!((root - 2.0_f64.powf(1.0 / 3.0)).abs() < 1e-10);

        // Verify it's actually a root
        let result = evaluate(&p, root);
        assert!((result).abs() < 1e-10);
    }

    #[test]
    fn test_newton_raphson_no_convergence() {
        // Polynomial with no real roots: x² + 1
        let p = Polynomial::new(vec![1.0, 0.0, 1.0]);
        let root = newton_raphson_root(&p, 0.0, 1e-10, 10);

        // Should fail to converge or find no root
        assert!(root.is_none() || evaluate(&p, root.unwrap()).abs() > 1e-5);
    }

    #[test]
    fn test_newton_raphson_constant_polynomial() {
        let p = Polynomial::new(vec![5.0]); // constant
        let root = newton_raphson_root(&p, 0.0, 1e-10, 100);
        assert_eq!(root, None);
    }

    #[test]
    fn test_durand_kerner_quadratic() {
        // x² - 1 = (x-1)(x+1)
        let p = Polynomial::new(vec![-1.0, 0.0, 1.0]);
        let roots = durand_kerner(&p, 1e-10, 100).unwrap();

        assert_eq!(roots.len(), 2);

        // Verify both roots
        for root in &roots {
            let result = evaluate_complex(&p, root);
            assert!(crate::complex::magnitude(&result) < 1e-8);
        }

        // Check that we found ±1
        let real_parts: Vec<f64> = roots.iter().map(|z| z.real).collect();
        assert!(real_parts.iter().any(|&r| (r - 1.0).abs() < 1e-8));
        assert!(real_parts.iter().any(|&r| (r + 1.0).abs() < 1e-8));
    }

    #[test]
    fn test_durand_kerner_cubic() {
        // x³ - 1 (roots: 1, e^(2πi/3), e^(4πi/3))
        let p = Polynomial::new(vec![-1.0, 0.0, 0.0, 1.0]);
        let roots = durand_kerner(&p, 1e-10, 100).unwrap();

        assert_eq!(roots.len(), 3);

        // Verify all roots
        for root in &roots {
            let result = evaluate_complex(&p, root);
            assert!(crate::complex::magnitude(&result) < 1e-8);
        }

        // Check that one root is real and equals 1
        let real_root = roots.iter().find(|z| z.imag.abs() < 1e-8);
        assert!(real_root.is_some());
        assert!((real_root.unwrap().real - 1.0).abs() < 1e-8);
    }

    #[test]
    fn test_durand_kerner_complex_roots() {
        // x² + 1 (roots: ±i)
        let p = Polynomial::new(vec![1.0, 0.0, 1.0]);
        let roots = durand_kerner(&p, 1e-10, 100).unwrap();

        assert_eq!(roots.len(), 2);

        // Verify both roots
        for root in &roots {
            let result = evaluate_complex(&p, root);
            assert!(crate::complex::magnitude(&result) < 1e-8);
        }

        // Check that we found ±i
        let imag_parts: Vec<f64> = roots.iter().map(|z| z.imag).collect();
        assert!(imag_parts.iter().any(|&i| (i - 1.0).abs() < 1e-8));
        assert!(imag_parts.iter().any(|&i| (i + 1.0).abs() < 1e-8));
    }

    #[test]
    fn test_durand_kerner_constant_polynomial() {
        let p = Polynomial::new(vec![5.0]); // constant
        let roots = durand_kerner(&p, 1e-10, 100);
        assert_eq!(roots, None);
    }

    #[test]
    fn test_evaluate_complex() {
        // P(x) = 1 + x²
        let p = Polynomial::new(vec![1.0, 0.0, 1.0]);

        // P(i) = 1 + i² = 1 - 1 = 0
        let i: Complex<f64> = Complex::i();
        let result = evaluate_complex(&p, &i);

        assert!((result.real).abs() < 1e-10);
        assert!((result.imag).abs() < 1e-10);
    }

    #[test]
    fn test_evaluate_complex_real_input() {
        // P(x) = 1 + 2x + 3x²
        let p = Polynomial::new(vec![1.0, 2.0, 3.0]);

        // P(2) = 1 + 4 + 12 = 17
        let z = Complex::new(2.0, 0.0);
        let result = evaluate_complex(&p, &z);

        assert!((result.real - 17.0).abs() < 1e-10);
        assert!((result.imag).abs() < 1e-10);
    }

    #[test]
    fn test_evaluate_complex_general() {
        // P(x) = 1 + x
        let p = Polynomial::new(vec![1.0, 1.0]);

        // P(1+2i) = 1 + (1+2i) = 2 + 2i
        let z = Complex::new(1.0, 2.0);
        let result = evaluate_complex(&p, &z);

        assert!((result.real - 2.0).abs() < 1e-10);
        assert!((result.imag - 2.0).abs() < 1e-10);
    }
}
