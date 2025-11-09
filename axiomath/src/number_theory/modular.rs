//! Modular arithmetic operations.
//!
//! This module provides modular arithmetic functions:
//! - Modular addition and multiplication
//! - Modular exponentiation (square-and-multiply algorithm)
//! - Modular inverse using extended Euclidean algorithm

use crate::number_theory::divisibility::extended_gcd;
use num_traits::sign::Signed;
use num_traits::{CheckedAdd, CheckedMul, One, PrimInt, ToPrimitive, Zero};

/// Computes (a + b) mod modulus with overflow protection.
///
/// # Examples
///
/// ```
/// use axiomath::number_theory::modular::modular_add;
///
/// assert_eq!(modular_add(10, 15, 12), 1);
/// assert_eq!(modular_add(5, 3, 7), 1);
/// ```
///
/// # Panics
///
/// Panics if modulus is zero.
#[inline]
pub fn modular_add<T>(a: T, b: T, modulus: T) -> T
where
    T: PrimInt + Copy + CheckedAdd,
{
    if modulus.is_zero() {
        panic!("Modulus cannot be zero");
    }

    let a_mod = a % modulus;
    let b_mod = b % modulus;

    // Handle overflow by using checked addition
    match a_mod.checked_add(&b_mod) {
        Some(sum) => sum % modulus,
        None => {
            // If overflow would occur, compute differently
            // (a + b) mod m = ((a mod m) + (b mod m)) mod m
            let diff = modulus - a_mod;
            if b_mod >= diff {
                b_mod - diff
            } else {
                a_mod + b_mod
            }
        }
    }
}

/// Computes (a × b) mod modulus with overflow protection.
///
/// # Examples
///
/// ```
/// use axiomath::number_theory::modular::modular_multiply;
///
/// assert_eq!(modular_multiply(7, 8, 12), 8);
/// assert_eq!(modular_multiply(123, 456, 100), 88);
/// ```
///
/// # Panics
///
/// Panics if modulus is zero.
#[inline]
pub fn modular_multiply<T>(a: T, b: T, modulus: T) -> T
where
    T: PrimInt + Copy + CheckedMul,
{
    if modulus.is_zero() {
        panic!("Modulus cannot be zero");
    }

    let a_mod = a % modulus;
    let b_mod = b % modulus;

    // Try checked multiplication first
    match a_mod.checked_mul(&b_mod) {
        Some(product) => product % modulus,
        None => {
            // If overflow would occur, use repeated addition
            // This is slower but safe
            modular_multiply_safe(a_mod, b_mod, modulus)
        }
    }
}

/// Helper function for modular multiplication using repeated addition (slow but safe).
#[inline]
fn modular_multiply_safe<T>(a: T, b: T, modulus: T) -> T
where
    T: PrimInt + Copy,
{
    let mut result = T::zero();
    let mut base = a % modulus;
    let mut exp = b;

    while exp > T::zero() {
        if exp % (T::one() + T::one()) == T::one() {
            result = (result + base) % modulus;
        }
        base = (base + base) % modulus;
        exp = exp / (T::one() + T::one());
    }

    result
}

/// Computes (base^exponent) mod modulus using the square-and-multiply algorithm.
///
/// # Algorithm
///
/// Binary exponentiation (exponentiation by squaring) in modular arithmetic:
/// - If exp is even: base^exp mod m = (base²)^(exp/2) mod m
/// - If exp is odd: base^exp mod m = base × base^(exp-1) mod m
///
/// Time complexity: O(log exponent)
///
/// # Examples
///
/// ```
/// use axiomath::number_theory::modular::modular_power;
///
/// assert_eq!(modular_power(3, 4, 5), 1); // 3^4 = 81 ≡ 1 (mod 5)
/// assert_eq!(modular_power(2, 10, 1000), 24); // 2^10 = 1024 ≡ 24 (mod 1000)
/// assert_eq!(modular_power(7, 100, 13), 9);
/// ```
///
/// # Panics
///
/// Panics if modulus is zero or if exponent is negative.
#[inline]
pub fn modular_power<T>(base: T, exponent: T, modulus: T) -> T
where
    T: PrimInt + Copy + CheckedMul,
{
    if modulus.is_zero() {
        panic!("Modulus cannot be zero");
    }
    if modulus == T::one() {
        return T::zero();
    }

    // Handle negative exponents
    if exponent < T::zero() {
        panic!("Negative exponents not supported in modular_power. Use modular_inverse first.");
    }

    if exponent.is_zero() {
        return T::one();
    }

    let mut result = T::one();
    let mut base_mod = base % modulus;
    let mut exp = exponent;
    let two = T::one() + T::one();

    while exp > T::zero() {
        // If exp is odd, multiply result by base
        if exp % two == T::one() {
            result = modular_multiply(result, base_mod, modulus);
        }

        // Square the base
        base_mod = modular_multiply(base_mod, base_mod, modulus);

        // Divide exponent by 2
        exp = exp / two;
    }

    result
}

/// Computes the modular multiplicative inverse of a modulo modulus.
///
/// Returns x such that (a × x) ≡ 1 (mod modulus).
///
/// # Algorithm
///
/// Uses the extended Euclidean algorithm. The inverse exists if and only if
/// gcd(a, modulus) = 1.
///
/// # Examples
///
/// ```
/// use axiomath::number_theory::modular::{modular_inverse, modular_multiply};
///
/// let inv = modular_inverse(3, 11).unwrap();
/// assert_eq!(inv, 4); // 3 × 4 = 12 ≡ 1 (mod 11)
///
/// // Verify: a × a^(-1) ≡ 1 (mod m)
/// assert_eq!(modular_multiply(3, inv, 11), 1);
///
/// let inv = modular_inverse(7, 26).unwrap();
/// assert_eq!(modular_multiply(7, inv, 26), 1);
/// ```
///
/// # Returns
///
/// - `Some(inverse)` if the inverse exists (when gcd(a, modulus) = 1)
/// - `None` if the inverse does not exist
///
/// # Panics
///
/// Panics if modulus is zero or one.
#[inline]
pub fn modular_inverse<T>(a: T, modulus: T) -> Option<T>
where
    T: PrimInt + Signed + Copy + ToPrimitive,
{
    if modulus.is_zero() || modulus == T::one() {
        panic!("Modulus must be greater than 1");
    }

    let (gcd, x, _) = extended_gcd(a, modulus);

    // Inverse exists only if gcd(a, modulus) = 1
    if gcd != T::one() {
        return None;
    }

    // Ensure the result is positive
    let inverse = if x < T::zero() {
        x + modulus
    } else {
        x
    };

    Some(inverse % modulus)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_modular_add_basic() {
        assert_eq!(modular_add(10, 15, 12), 1);
        assert_eq!(modular_add(5, 3, 7), 1);
        assert_eq!(modular_add(0, 5, 10), 5);
        assert_eq!(modular_add(7, 7, 10), 4);
    }

    #[test]
    #[should_panic(expected = "Modulus cannot be zero")]
    fn test_modular_add_zero_modulus() {
        modular_add(5, 3, 0);
    }

    #[test]
    fn test_modular_multiply_basic() {
        assert_eq!(modular_multiply(7, 8, 12), 8);
        assert_eq!(modular_multiply(123, 456, 100), 88);
        assert_eq!(modular_multiply(5, 3, 7), 1);
        assert_eq!(modular_multiply(0, 5, 10), 0);
    }

    #[test]
    #[should_panic(expected = "Modulus cannot be zero")]
    fn test_modular_multiply_zero_modulus() {
        modular_multiply(5, 3, 0);
    }

    #[test]
    fn test_modular_power_basic() {
        assert_eq!(modular_power(3, 4, 5), 1); // 3^4 = 81 ≡ 1 (mod 5)
        assert_eq!(modular_power(2, 10, 1000), 24); // 2^10 = 1024 ≡ 24 (mod 1000)
        assert_eq!(modular_power(7, 100, 13), 9);
    }

    #[test]
    fn test_modular_power_zero_exponent() {
        assert_eq!(modular_power(5, 0, 13), 1);
        assert_eq!(modular_power(123, 0, 100), 1);
    }

    #[test]
    fn test_modular_power_one_exponent() {
        assert_eq!(modular_power(5, 1, 13), 5);
        assert_eq!(modular_power(123, 1, 100), 23);
    }

    #[test]
    fn test_modular_power_fermat() {
        // Fermat's Little Theorem: a^(p-1) ≡ 1 (mod p) for prime p and gcd(a, p) = 1
        let p = 13; // Prime
        for a in 1..p {
            assert_eq!(modular_power(a, p - 1, p), 1);
        }
    }

    #[test]
    #[should_panic(expected = "Modulus cannot be zero")]
    fn test_modular_power_zero_modulus() {
        modular_power(5, 3, 0);
    }

    #[test]
    fn test_modular_inverse_basic() {
        let inv = modular_inverse(3, 11).unwrap();
        assert_eq!(inv, 4); // 3 × 4 = 12 ≡ 1 (mod 11)
        assert_eq!(modular_multiply(3, inv, 11), 1);

        let inv = modular_inverse(7, 26).unwrap();
        assert_eq!(modular_multiply(7, inv, 26), 1);
    }

    #[test]
    fn test_modular_inverse_coprime() {
        // Test various coprime pairs
        let test_cases = [(3, 11), (5, 17), (7, 13), (11, 23)];

        for (a, m) in test_cases {
            if let Some(inv) = modular_inverse(a, m) {
                assert_eq!(modular_multiply(a, inv, m), 1);
            }
        }
    }

    #[test]
    fn test_modular_inverse_no_inverse() {
        // 6 and 9 are not coprime (gcd = 3)
        assert_eq!(modular_inverse(6, 9), None);

        // Even numbers don't have inverses modulo even numbers
        assert_eq!(modular_inverse(4, 8), None);
    }

    #[test]
    #[should_panic(expected = "Modulus must be greater than 1")]
    fn test_modular_inverse_invalid_modulus() {
        modular_inverse(5, 1);
    }

    #[test]
    fn test_modular_arithmetic_properties() {
        // Test: (a + b) mod m = ((a mod m) + (b mod m)) mod m
        let a = 123;
        let b = 456;
        let m = 100;

        let sum1 = (a + b) % m;
        let sum2 = modular_add(a, b, m);
        assert_eq!(sum1, sum2);

        // Test: (a × b) mod m = ((a mod m) × (b mod m)) mod m
        let prod1 = (a * b) % m;
        let prod2 = modular_multiply(a, b, m);
        assert_eq!(prod1, prod2);
    }

    #[test]
    fn test_modular_power_vs_naive() {
        // Compare modular_power with naive exponentiation for small values
        for base in 1..10 {
            for exp in 0..10 {
                for modulus in 2..20 {
                    let result1 = modular_power(base, exp, modulus);
                    let result2 = base.pow(exp as u32) % modulus;
                    assert_eq!(result1, result2);
                }
            }
        }
    }
}
