//! Divisibility operations including GCD and LCM.
//!
//! This module provides functions for divisibility testing and computation:
//! - Greatest Common Divisor (GCD) using Euclidean algorithm
//! - Least Common Multiple (LCM)
//! - Extended Euclidean algorithm

use num_traits::sign::Signed;
use num_traits::{CheckedMul, PrimInt, Zero};

/// Computes the greatest common divisor (GCD) of two integers using the Euclidean algorithm.
///
/// # Algorithm
///
/// The Euclidean algorithm is based on the principle that:
/// ```text
/// gcd(a, b) = gcd(b, a mod b)
/// ```
///
/// The algorithm repeatedly applies this until b = 0, at which point gcd(a, 0) = |a|.
///
/// Time complexity: O(log(min(a, b)))
///
/// # Examples
///
/// ```
/// use axiomath::number_theory::divisibility::greatest_common_divisor;
///
/// assert_eq!(greatest_common_divisor(48, 18), 6);
/// assert_eq!(greatest_common_divisor(100, 75), 25);
/// assert_eq!(greatest_common_divisor(17, 19), 1); // Coprime
/// assert_eq!(greatest_common_divisor(0, 5), 5);
/// ```
#[inline]
pub fn greatest_common_divisor<T>(mut a: T, mut b: T) -> T
where
    T: PrimInt + Signed,
{
    // Handle zero cases
    if a.is_zero() {
        return b.abs();
    }
    if b.is_zero() {
        return a.abs();
    }

    // Make both positive
    a = a.abs();
    b = b.abs();

    // Euclidean algorithm
    while !b.is_zero() {
        let temp = b;
        b = a % b;
        a = temp;
    }

    a
}

/// Computes the least common multiple (LCM) of two integers.
///
/// # Algorithm
///
/// Uses the relationship between GCD and LCM:
/// ```text
/// lcm(a, b) = |a × b| / gcd(a, b)
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::number_theory::divisibility::least_common_multiple;
///
/// assert_eq!(least_common_multiple(12, 18), 36);
/// assert_eq!(least_common_multiple(21, 6), 42);
/// assert_eq!(least_common_multiple(7, 5), 35);
/// ```
///
/// # Panics
///
/// Panics if both arguments are zero or if the result would overflow.
#[inline]
pub fn least_common_multiple<T>(a: T, b: T) -> T
where
    T: PrimInt + Signed + CheckedMul,
{
    if a.is_zero() && b.is_zero() {
        panic!("LCM of (0, 0) is undefined");
    }
    if a.is_zero() || b.is_zero() {
        return T::zero();
    }

    let gcd = greatest_common_divisor(a, b);
    let a_abs = a.abs();
    let b_abs = b.abs();

    // Calculate lcm = |a * b| / gcd
    // To avoid overflow, compute (a / gcd) * b
    let a_div_gcd = a_abs / gcd;

    a_div_gcd
        .checked_mul(&b_abs)
        .expect("LCM computation would overflow")
}

/// Computes the extended Euclidean algorithm.
///
/// Returns a tuple (gcd, x, y) such that:
/// ```text
/// gcd = a × x + b × y
/// ```
///
/// where gcd is the greatest common divisor of a and b.
///
/// This is useful for finding modular inverses and solving linear Diophantine equations.
///
/// # Algorithm
///
/// The extended Euclidean algorithm maintains the invariants:
/// - r_i = s_i × a + t_i × b
///
/// at each step of the Euclidean algorithm.
///
/// # Examples
///
/// ```
/// use axiomath::number_theory::divisibility::extended_gcd;
///
/// let (gcd, x, y) = extended_gcd(35, 15);
/// assert_eq!(gcd, 5);
/// assert_eq!(35 * x + 15 * y, gcd); // Bézout's identity
///
/// let (gcd, x, y) = extended_gcd(240, 46);
/// assert_eq!(gcd, 2);
/// assert_eq!(240 * x + 46 * y, gcd);
/// ```
#[inline]
pub fn extended_gcd<T>(a: T, b: T) -> (T, T, T)
where
    T: PrimInt + Signed,
{
    if b.is_zero() {
        // gcd(a, 0) = a = a × 1 + 0 × 0
        return (a, T::one(), T::zero());
    }

    let mut r0 = a;
    let mut r1 = b;
    let mut s0 = T::one();
    let mut s1 = T::zero();
    let mut t0 = T::zero();
    let mut t1 = T::one();

    while !r1.is_zero() {
        let quotient = r0 / r1;

        let r_temp = r0 - quotient * r1;
        r0 = r1;
        r1 = r_temp;

        let s_temp = s0 - quotient * s1;
        s0 = s1;
        s1 = s_temp;

        let t_temp = t0 - quotient * t1;
        t0 = t1;
        t1 = t_temp;
    }

    (r0, s0, t0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gcd_basic() {
        assert_eq!(greatest_common_divisor(48, 18), 6);
        assert_eq!(greatest_common_divisor(100, 75), 25);
        assert_eq!(greatest_common_divisor(54, 24), 6);
        assert_eq!(greatest_common_divisor(1071, 462), 21);
    }

    #[test]
    fn test_gcd_coprime() {
        // Coprime numbers have gcd = 1
        assert_eq!(greatest_common_divisor(17, 19), 1);
        assert_eq!(greatest_common_divisor(13, 28), 1);
        assert_eq!(greatest_common_divisor(35, 64), 1);
    }

    #[test]
    fn test_gcd_zero() {
        assert_eq!(greatest_common_divisor(0, 5), 5);
        assert_eq!(greatest_common_divisor(5, 0), 5);
        assert_eq!(greatest_common_divisor(0, 0), 0);
    }

    #[test]
    fn test_gcd_negative() {
        assert_eq!(greatest_common_divisor(-48, 18), 6);
        assert_eq!(greatest_common_divisor(48, -18), 6);
        assert_eq!(greatest_common_divisor(-48, -18), 6);
    }

    #[test]
    fn test_gcd_same() {
        assert_eq!(greatest_common_divisor(42, 42), 42);
        assert_eq!(greatest_common_divisor(7, 7), 7);
    }

    #[test]
    fn test_lcm_basic() {
        assert_eq!(least_common_multiple(12, 18), 36);
        assert_eq!(least_common_multiple(21, 6), 42);
        assert_eq!(least_common_multiple(4, 6), 12);
    }

    #[test]
    fn test_lcm_coprime() {
        // For coprime numbers, lcm(a, b) = a × b
        assert_eq!(least_common_multiple(7, 5), 35);
        assert_eq!(least_common_multiple(13, 17), 221);
    }

    #[test]
    fn test_lcm_zero() {
        assert_eq!(least_common_multiple(0, 5), 0);
        assert_eq!(least_common_multiple(5, 0), 0);
    }

    #[test]
    #[should_panic(expected = "LCM of (0, 0) is undefined")]
    fn test_lcm_both_zero() {
        least_common_multiple(0, 0);
    }

    #[test]
    fn test_lcm_gcd_relationship() {
        // For all a, b: gcd(a, b) × lcm(a, b) = |a × b|
        let test_pairs = [(12, 18), (21, 6), (35, 49), (100, 75)];

        for (a, b) in test_pairs {
            let gcd = greatest_common_divisor(a, b);
            let lcm = least_common_multiple(a, b);
            assert_eq!(gcd * lcm, a * b);
        }
    }

    #[test]
    fn test_extended_gcd_basic() {
        let (gcd, x, y) = extended_gcd(35, 15);
        assert_eq!(gcd, 5);
        assert_eq!(35 * x + 15 * y, gcd); // Bézout's identity

        let (gcd, x, y) = extended_gcd(240, 46);
        assert_eq!(gcd, 2);
        assert_eq!(240 * x + 46 * y, gcd);
    }

    #[test]
    fn test_extended_gcd_coprime() {
        let (gcd, x, y) = extended_gcd(17, 13);
        assert_eq!(gcd, 1);
        assert_eq!(17 * x + 13 * y, 1);

        let (gcd, x, y) = extended_gcd(35, 64);
        assert_eq!(gcd, 1);
        assert_eq!(35 * x + 64 * y, 1);
    }

    #[test]
    fn test_extended_gcd_zero() {
        let (gcd, x, y) = extended_gcd(42, 0);
        assert_eq!(gcd, 42);
        assert_eq!(x, 1);
        assert_eq!(y, 0);
    }

    #[test]
    fn test_extended_gcd_consistency() {
        // extended_gcd should give same gcd as greatest_common_divisor
        let test_pairs = [(48, 18), (100, 75), (1071, 462)];

        for (a, b) in test_pairs {
            let (gcd_ext, x, y) = extended_gcd(a, b);
            let gcd = greatest_common_divisor(a, b);
            assert_eq!(gcd_ext, gcd);
            assert_eq!(a * x + b * y, gcd);
        }
    }
}
