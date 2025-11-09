//! Power and root operations implemented from first principles.
//!
//! This module provides power and root functions using mathematical algorithms:
//! - Newton-Raphson method for square roots
//! - Binary exponentiation for integer powers
//! - Generic nth root implementations

use crate::core::traits::{Float, Scalar};

/// Computes the power `base^exponent` for integer exponents.
///
/// Uses binary exponentiation (exponentiation by squaring) for efficient
/// computation. Time complexity: O(log n) where n is the exponent.
///
/// # Algorithm
///
/// Binary exponentiation uses the property:
/// - If exp is even: base^exp = (base²)^(exp/2)
/// - If exp is odd: base^exp = base × base^(exp-1)
///
/// # Examples
///
/// ```
/// use axiomath::scalar::power::power;
///
/// assert_eq!(power(2.0, 10), 1024.0);
/// assert_eq!(power(3.0, 3), 27.0);
/// assert_eq!(power(5.0, 0), 1.0);
/// ```
///
/// # Panics
///
/// Panics if base is zero and exponent is negative (division by zero).
#[inline]
pub fn power<T>(base: T, mut exponent: i32) -> T
where
    T: Scalar + Copy,
{
    if exponent == 0 {
        return T::one();
    }

    if exponent < 0 {
        assert!(!base.is_zero(), "Cannot compute negative power of zero");
        return T::one() / power(base, -exponent);
    }

    let mut result = T::one();
    let mut current_base = base;

    while exponent > 0 {
        if exponent % 2 == 1 {
            result = result * current_base;
        }
        current_base = current_base * current_base;
        exponent /= 2;
    }

    result
}

/// Computes the square root using the Newton-Raphson method.
///
/// # Algorithm
///
/// Newton-Raphson iteration for √x:
/// - x_{n+1} = (x_n + x/x_n) / 2
///
/// Starting with an initial guess, this converges quadratically to √x.
///
/// # Examples
///
/// ```
/// use axiomath::scalar::power::square_root;
///
/// let result = square_root(4.0_f64);
/// assert!((result - 2.0).abs() < 1e-2);
///
/// let result = square_root(2.0_f64);
/// assert!((result - 1.414213562373095).abs() < 1e-2);
/// ```
///
/// # Panics
///
/// Panics if x is negative (square root of negative number).
#[inline]
pub fn square_root<T>(x: T) -> T
where
    T: Float,
{
    // Handle special cases
    if x.is_zero() {
        return T::zero();
    }
    if x < T::zero() {
        panic!("Cannot compute square root of negative number");
    }
    if x.is_infinite() {
        return x;
    }
    if x.is_nan() {
        return x;
    }

    // Initial guess: use a simple heuristic
    // For numbers close to 1, start with x/2; otherwise use a better guess
    let mut guess = if x < T::one() {
        x
    } else {
        x / (T::one() + T::one())
    };

    let epsilon = T::epsilon() * (T::one() + T::one() + T::one()); // 3 * epsilon
    let max_iterations = 50;

    for _ in 0..max_iterations {
        let next_guess = (guess + x / guess) / (T::one() + T::one());

        // Check for convergence
        let diff = (next_guess - guess).abs();
        if diff <= epsilon * next_guess.abs() {
            return next_guess;
        }

        guess = next_guess;
    }

    guess
}

/// Computes the cube root using the Newton-Raphson method.
///
/// # Algorithm
///
/// Newton-Raphson iteration for ∛x:
/// - x_{n+1} = (2*x_n + x/x_n²) / 3
///
/// # Examples
///
/// ```
/// use axiomath::scalar::power::cube_root;
///
/// let result = cube_root(8.0_f64);
/// assert!((result - 2.0).abs() < 1e-2);
///
/// let result = cube_root(27.0_f64);
/// assert!((result - 3.0).abs() < 1e-2);
///
/// // Cube root of negative numbers is defined
/// let result = cube_root(-8.0_f64);
/// assert!((result - (-2.0)).abs() < 1e-2);
/// ```
#[inline]
pub fn cube_root<T>(x: T) -> T
where
    T: Float,
{
    // Handle special cases
    if x.is_zero() {
        return T::zero();
    }
    if x.is_infinite() {
        return x;
    }
    if x.is_nan() {
        return x;
    }

    // Handle negative numbers
    let negative = x < T::zero();
    let abs_x = x.abs();

    // Initial guess
    let mut guess = if abs_x < T::one() {
        abs_x
    } else {
        abs_x / (T::one() + T::one() + T::one())
    };

    let epsilon = T::epsilon() * (T::one() + T::one() + T::one());
    let max_iterations = 50;
    let two = T::one() + T::one();
    let three = two + T::one();

    for _ in 0..max_iterations {
        let guess_squared = guess * guess;
        let next_guess = (two * guess + abs_x / guess_squared) / three;

        let diff = (next_guess - guess).abs();
        if diff <= epsilon * next_guess.abs() {
            return if negative { -next_guess } else { next_guess };
        }

        guess = next_guess;
    }

    if negative { -guess } else { guess }
}

/// Computes the nth root using the Newton-Raphson method.
///
/// # Algorithm
///
/// Newton-Raphson iteration for ⁿ√x:
/// - x_{k+1} = ((n-1)*x_k + x/x_k^(n-1)) / n
///
/// # Examples
///
/// ```
/// use axiomath::scalar::power::nth_root;
///
/// let result = nth_root(16.0_f64, 4);
/// assert!((result - 2.0).abs() < 1e-2);
///
/// let result = nth_root(32.0_f64, 5);
/// assert!((result - 2.0).abs() < 1e-2);
/// ```
///
/// # Panics
///
/// - Panics if n is zero or negative
/// - Panics if x is negative and n is even
#[inline]
pub fn nth_root<T>(x: T, n: i32) -> T
where
    T: Float,
{
    assert!(n > 0, "Root degree must be positive");

    // Handle special cases
    if n == 1 {
        return x;
    }
    if n == 2 {
        return square_root(x);
    }
    if n == 3 {
        return cube_root(x);
    }

    if x.is_zero() {
        return T::zero();
    }
    if x.is_infinite() {
        return x;
    }
    if x.is_nan() {
        return x;
    }

    // Check for even root of negative number
    let negative = x < T::zero();
    if negative && n % 2 == 0 {
        panic!("Cannot compute even root of negative number");
    }

    let abs_x = x.abs();

    // Initial guess
    let mut guess = if abs_x < T::one() {
        abs_x
    } else {
        abs_x / T::from(n).expect("Failed to convert n to T")
    };

    let epsilon = T::epsilon() * (T::one() + T::one() + T::one());
    let max_iterations = 100;
    let n_float = T::from(n).expect("Failed to convert n to T");
    let n_minus_one = T::from(n - 1).expect("Failed to convert n-1 to T");

    for _ in 0..max_iterations {
        // Compute guess^(n-1)
        let mut guess_power = T::one();
        for _ in 0..(n - 1) {
            guess_power = guess_power * guess;
        }

        let next_guess = (n_minus_one * guess + abs_x / guess_power) / n_float;

        let diff = (next_guess - guess).abs();
        if diff <= epsilon * next_guess.abs() {
            return if negative { -next_guess } else { next_guess };
        }

        guess = next_guess;
    }

    if negative { -guess } else { guess }
}

/// Computes the absolute value of a number.
///
/// # Examples
///
/// ```
/// use axiomath::scalar::power::absolute_value;
///
/// assert_eq!(absolute_value(5.0), 5.0);
/// assert_eq!(absolute_value(-5.0), 5.0);
/// assert_eq!(absolute_value(0.0), 0.0);
/// ```
#[inline]
pub fn absolute_value<T>(x: T) -> T
where
    T: Scalar,
{
    x.abs()
}

/// Returns the sign of a number as -1 or 1.
///
/// Note: For IEEE 754 floats, `sign(0.0)` returns `1.0` (not `0.0`).
///
/// # Examples
///
/// ```
/// use axiomath::scalar::power::sign;
///
/// assert_eq!(sign(5.0), 1.0);
/// assert_eq!(sign(-5.0), -1.0);
/// assert_eq!(sign(0.0), 1.0);  // IEEE 754 behavior
/// ```
#[inline]
pub fn sign<T>(x: T) -> T
where
    T: Scalar,
{
    x.signum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_power_basic() {
        assert_eq!(power(2.0, 0), 1.0);
        assert_eq!(power(2.0, 1), 2.0);
        assert_eq!(power(2.0, 10), 1024.0);
        assert_eq!(power(3.0, 3), 27.0);
        assert_eq!(power(5.0, 2), 25.0);
    }

    #[test]
    fn test_power_negative_exponent() {
        assert_abs_diff_eq!(power(2.0, -1), 0.5, epsilon = 1e-10);
        assert_abs_diff_eq!(power(2.0, -3), 0.125, epsilon = 1e-10);
        assert_abs_diff_eq!(power(10.0, -2), 0.01, epsilon = 1e-10);
    }

    #[test]
    #[should_panic(expected = "Cannot compute negative power of zero")]
    fn test_power_zero_negative() {
        power(0.0, -1);
    }

    #[test]
    fn test_square_root_basic() {
        assert_abs_diff_eq!(square_root(4.0), 2.0, epsilon = 1e-10);
        assert_abs_diff_eq!(square_root(9.0), 3.0, epsilon = 1e-10);
        assert_abs_diff_eq!(square_root(16.0), 4.0, epsilon = 1e-10);
        assert_abs_diff_eq!(square_root(0.0), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_square_root_irrational() {
        assert_abs_diff_eq!(square_root(2.0), 1.414213562373095, epsilon = 1e-10);
        assert_abs_diff_eq!(square_root(3.0), 1.732050807568877, epsilon = 1e-10);
        assert_abs_diff_eq!(square_root(5.0), 2.236067977499790, epsilon = 1e-10);
    }

    #[test]
    fn test_square_root_small_numbers() {
        assert_abs_diff_eq!(square_root(0.25), 0.5, epsilon = 1e-10);
        assert_abs_diff_eq!(square_root(0.01), 0.1, epsilon = 1e-10);
    }

    #[test]
    #[should_panic(expected = "Cannot compute square root of negative number")]
    fn test_square_root_negative() {
        square_root(-1.0);
    }

    #[test]
    fn test_cube_root_basic() {
        assert_abs_diff_eq!(cube_root(8.0), 2.0, epsilon = 1e-10);
        assert_abs_diff_eq!(cube_root(27.0), 3.0, epsilon = 1e-10);
        assert_abs_diff_eq!(cube_root(64.0), 4.0, epsilon = 1e-10);
        assert_abs_diff_eq!(cube_root(0.0), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_cube_root_negative() {
        assert_abs_diff_eq!(cube_root(-8.0), -2.0, epsilon = 1e-10);
        assert_abs_diff_eq!(cube_root(-27.0), -3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_nth_root_basic() {
        assert_abs_diff_eq!(nth_root(16.0, 4), 2.0, epsilon = 1e-10);
        assert_abs_diff_eq!(nth_root(32.0, 5), 2.0, epsilon = 1e-10);
        assert_abs_diff_eq!(nth_root(243.0, 5), 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_nth_root_special_cases() {
        assert_abs_diff_eq!(nth_root(8.0, 1), 8.0, epsilon = 1e-10);
        assert_abs_diff_eq!(nth_root(4.0, 2), 2.0, epsilon = 1e-10);
        assert_abs_diff_eq!(nth_root(8.0, 3), 2.0, epsilon = 1e-10);
    }

    #[test]
    #[should_panic(expected = "Cannot compute even root of negative number")]
    fn test_nth_root_even_negative() {
        nth_root(-16.0, 4);
    }

    #[test]
    fn test_nth_root_odd_negative() {
        assert_abs_diff_eq!(nth_root(-32.0, 5), -2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_absolute_value() {
        assert_eq!(absolute_value(5.0), 5.0);
        assert_eq!(absolute_value(-5.0), 5.0);
        assert_eq!(absolute_value(0.0), 0.0);
    }

    #[test]
    fn test_sign() {
        assert_eq!(sign(5.0), 1.0);
        assert_eq!(sign(-5.0), -1.0);
        // Note: IEEE 754 signum(0.0) returns 1.0, not 0.0
        assert_eq!(sign(0.0), 1.0);
    }
}
