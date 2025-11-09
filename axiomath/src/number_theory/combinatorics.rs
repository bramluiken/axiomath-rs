//! Combinatorial functions including factorials and binomial coefficients.
//!
//! This module provides combinatorial counting functions:
//! - Factorial with overflow protection
//! - Binomial coefficients (n choose k)
//! - Permutations
//! - Stirling numbers

use num_traits::{CheckedMul, One, PrimInt, ToPrimitive, Zero};

/// Computes the factorial of n: n! = n × (n-1) × ... × 2 × 1
///
/// # Algorithm
///
/// Iterative computation with overflow checking.
///
/// # Examples
///
/// ```
/// use axiomath::number_theory::combinatorics::factorial;
///
/// assert_eq!(factorial(0_u64), 1);
/// assert_eq!(factorial(5_u64), 120);
/// assert_eq!(factorial(10_u64), 3628800);
/// ```
///
/// # Panics
///
/// Panics if the result would overflow the integer type.
#[inline]
pub fn factorial<T>(n: T) -> T
where
    T: PrimInt + CheckedMul,
{
    if n.is_zero() || n.is_one() {
        return T::one();
    }

    let mut result = T::one();
    let mut i = T::one() + T::one(); // Start from 2

    while i <= n {
        result = result
            .checked_mul(&i)
            .expect("Factorial computation would overflow");
        i = i + T::one();
    }

    result
}

/// Computes the binomial coefficient C(n, k) = n! / (k! × (n-k)!)
///
/// Also known as "n choose k", this counts the number of ways to choose
/// k items from n items without regard to order.
///
/// # Algorithm
///
/// Uses the multiplicative formula to avoid computing large factorials:
/// ```text
/// C(n, k) = ∏(i=1 to k) (n - k + i) / i
/// ```
///
/// This is more efficient and less prone to overflow than computing factorials directly.
///
/// # Examples
///
/// ```
/// use axiomath::number_theory::combinatorics::binomial_coefficient;
///
/// assert_eq!(binomial_coefficient(5, 2), 10);
/// assert_eq!(binomial_coefficient(10, 3), 120);
/// assert_eq!(binomial_coefficient(6, 0), 1);
/// assert_eq!(binomial_coefficient(6, 6), 1);
/// ```
///
/// # Panics
///
/// Panics if k > n or if the result would overflow.
#[inline]
pub fn binomial_coefficient<T>(n: T, k: T) -> T
where
    T: PrimInt + CheckedMul + Copy,
{
    if k > n {
        panic!("k must be less than or equal to n");
    }

    // C(n, 0) = C(n, n) = 1
    if k.is_zero() || k == n {
        return T::one();
    }

    // Use symmetry: C(n, k) = C(n, n-k)
    // Choose the smaller k for efficiency
    let k = if k > n - k { n - k } else { k };

    let mut result = T::one();
    let mut i = T::one();

    while i <= k {
        result = result
            .checked_mul(&(n - k + i))
            .expect("Binomial coefficient computation would overflow");
        result = result / i;
        i = i + T::one();
    }

    result
}

/// Computes the number of permutations P(n, k) = n! / (n-k)!
///
/// This counts the number of ways to arrange k items from n items
/// where order matters.
///
/// # Algorithm
///
/// ```text
/// P(n, k) = n × (n-1) × ... × (n-k+1)
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::number_theory::combinatorics::permutations;
///
/// assert_eq!(permutations(5, 2), 20);
/// assert_eq!(permutations(10, 3), 720);
/// assert_eq!(permutations(5, 5), 120); // Same as factorial(5)
/// assert_eq!(permutations(5, 0), 1);
/// ```
///
/// # Panics
///
/// Panics if k > n or if the result would overflow.
#[inline]
pub fn permutations<T>(n: T, k: T) -> T
where
    T: PrimInt + CheckedMul + Copy,
{
    if k > n {
        panic!("k must be less than or equal to n");
    }

    if k.is_zero() {
        return T::one();
    }

    let mut result = T::one();
    let mut i = T::zero();

    while i < k {
        result = result
            .checked_mul(&(n - i))
            .expect("Permutations computation would overflow");
        i = i + T::one();
    }

    result
}

/// Computes the Stirling number of the first kind s(n, k).
///
/// Stirling numbers of the first kind count the number of permutations
/// of n elements with k disjoint cycles.
///
/// # Algorithm
///
/// Uses the recurrence relation:
/// ```text
/// s(n, k) = s(n-1, k-1) - (n-1) × s(n-1, k)
/// ```
///
/// with base cases:
/// - s(0, 0) = 1
/// - s(n, 0) = 0 for n > 0
/// - s(0, k) = 0 for k > 0
///
/// # Examples
///
/// ```
/// use axiomath::number_theory::combinatorics::stirling_number_first_kind;
///
/// assert_eq!(stirling_number_first_kind(3, 2), 3);
/// assert_eq!(stirling_number_first_kind(4, 2), 11);
/// ```
///
/// # Note
///
/// This implementation computes unsigned Stirling numbers of the first kind.
/// The signed version alternates signs: (-1)^(n-k) × s(n, k).
#[inline]
pub fn stirling_number_first_kind<T>(n: T, k: T) -> T
where
    T: PrimInt + Copy,
{
    // Base cases
    if n.is_zero() && k.is_zero() {
        return T::one();
    }
    if n.is_zero() || k.is_zero() {
        return T::zero();
    }
    if k > n {
        return T::zero();
    }
    if k == n {
        return T::one();
    }

    // Build table using dynamic programming
    let n_usize = to_usize(n);
    let k_usize = to_usize(k);

    let mut table = vec![vec![T::zero(); k_usize + 1]; n_usize + 1];
    table[0][0] = T::one();

    for i in 1..=n_usize {
        for j in 1..=k_usize.min(i) {
            if j == 1 {
                table[i][j] = factorial(from_usize::<T>(i - 1));
            } else if i == j {
                table[i][j] = T::one();
            } else {
                let i_minus_1 = from_usize::<T>(i - 1);
                table[i][j] = table[i - 1][j - 1] + i_minus_1 * table[i - 1][j];
            }
        }
    }

    table[n_usize][k_usize]
}

/// Computes the Stirling number of the second kind S(n, k).
///
/// Stirling numbers of the second kind count the number of ways to
/// partition n elements into k non-empty subsets.
///
/// # Algorithm
///
/// Uses the recurrence relation:
/// ```text
/// S(n, k) = k × S(n-1, k) + S(n-1, k-1)
/// ```
///
/// with base cases:
/// - S(0, 0) = 1
/// - S(n, 0) = 0 for n > 0
/// - S(0, k) = 0 for k > 0
///
/// # Examples
///
/// ```
/// use axiomath::number_theory::combinatorics::stirling_number_second_kind;
///
/// assert_eq!(stirling_number_second_kind(3, 2), 3);
/// assert_eq!(stirling_number_second_kind(4, 2), 7);
/// assert_eq!(stirling_number_second_kind(5, 3), 25);
/// ```
#[inline]
pub fn stirling_number_second_kind<T>(n: T, k: T) -> T
where
    T: PrimInt + Copy,
{
    // Base cases
    if n.is_zero() && k.is_zero() {
        return T::one();
    }
    if n.is_zero() || k.is_zero() {
        return T::zero();
    }
    if k > n {
        return T::zero();
    }
    if k == n {
        return T::one();
    }
    if k == T::one() {
        return T::one();
    }

    // Build table using dynamic programming
    let n_usize = to_usize(n);
    let k_usize = to_usize(k);

    let mut table = vec![vec![T::zero(); k_usize + 1]; n_usize + 1];
    table[0][0] = T::one();

    for i in 1..=n_usize {
        for j in 1..=k_usize.min(i) {
            if i == j {
                table[i][j] = T::one();
            } else if j == 1 {
                table[i][j] = T::one();
            } else {
                let j_t = from_usize::<T>(j);
                table[i][j] = j_t * table[i - 1][j] + table[i - 1][j - 1];
            }
        }
    }

    table[n_usize][k_usize]
}

// Helper functions for type conversion
#[inline]
fn to_usize<T: PrimInt>(n: T) -> usize {
    n.to_usize().expect("Value too large to convert to usize")
}

#[inline]
fn from_usize<T: PrimInt>(n: usize) -> T {
    T::from(n).expect("Failed to convert from usize")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_factorial_basic() {
        assert_eq!(factorial(0_u64), 1);
        assert_eq!(factorial(1_u64), 1);
        assert_eq!(factorial(2_u64), 2);
        assert_eq!(factorial(3_u64), 6);
        assert_eq!(factorial(4_u64), 24);
        assert_eq!(factorial(5_u64), 120);
        assert_eq!(factorial(6_u64), 720);
        assert_eq!(factorial(10_u64), 3628800);
    }

    #[test]
    fn test_binomial_coefficient_basic() {
        assert_eq!(binomial_coefficient(5, 2), 10);
        assert_eq!(binomial_coefficient(5, 3), 10); // Symmetry
        assert_eq!(binomial_coefficient(10, 3), 120);
        assert_eq!(binomial_coefficient(6, 0), 1);
        assert_eq!(binomial_coefficient(6, 6), 1);
    }

    #[test]
    fn test_binomial_coefficient_symmetry() {
        // C(n, k) = C(n, n-k)
        assert_eq!(
            binomial_coefficient(10, 3),
            binomial_coefficient(10, 7)
        );
        assert_eq!(
            binomial_coefficient(8, 2),
            binomial_coefficient(8, 6)
        );
    }

    #[test]
    fn test_binomial_coefficient_pascal() {
        // Pascal's identity: C(n, k) = C(n-1, k-1) + C(n-1, k)
        for n in 2..10 {
            for k in 1..n {
                let left = binomial_coefficient(n, k);
                let right = binomial_coefficient(n - 1, k - 1) + binomial_coefficient(n - 1, k);
                assert_eq!(left, right);
            }
        }
    }

    #[test]
    #[should_panic(expected = "k must be less than or equal to n")]
    fn test_binomial_coefficient_invalid() {
        binomial_coefficient(5, 7);
    }

    #[test]
    fn test_permutations_basic() {
        assert_eq!(permutations(5, 0), 1);
        assert_eq!(permutations(5, 1), 5);
        assert_eq!(permutations(5, 2), 20);
        assert_eq!(permutations(5, 5), 120);
        assert_eq!(permutations(10, 3), 720);
    }

    #[test]
    fn test_permutations_relationship() {
        // P(n, k) = C(n, k) × k!
        for n in 5..10 {
            for k in 0..=n {
                let perm = permutations(n, k);
                let comb = binomial_coefficient(n, k);
                let fact = factorial(k);
                assert_eq!(perm, comb * fact);
            }
        }
    }

    #[test]
    #[should_panic(expected = "k must be less than or equal to n")]
    fn test_permutations_invalid() {
        permutations(5, 7);
    }

    #[test]
    fn test_stirling_first_kind() {
        // Known values
        assert_eq!(stirling_number_first_kind(3, 2), 3);
        assert_eq!(stirling_number_first_kind(4, 2), 11);
        assert_eq!(stirling_number_first_kind(5, 2), 50);
        assert_eq!(stirling_number_first_kind(4, 3), 6);
    }

    #[test]
    fn test_stirling_first_kind_edge_cases() {
        assert_eq!(stirling_number_first_kind(0, 0), 1);
        assert_eq!(stirling_number_first_kind(5, 0), 0);
        assert_eq!(stirling_number_first_kind(0, 5), 0);
        assert_eq!(stirling_number_first_kind(5, 5), 1);
    }

    #[test]
    fn test_stirling_second_kind() {
        // Known values
        assert_eq!(stirling_number_second_kind(3, 2), 3);
        assert_eq!(stirling_number_second_kind(4, 2), 7);
        assert_eq!(stirling_number_second_kind(5, 2), 15);
        assert_eq!(stirling_number_second_kind(5, 3), 25);
        assert_eq!(stirling_number_second_kind(6, 3), 90);
    }

    #[test]
    fn test_stirling_second_kind_edge_cases() {
        assert_eq!(stirling_number_second_kind(0, 0), 1);
        assert_eq!(stirling_number_second_kind(5, 0), 0);
        assert_eq!(stirling_number_second_kind(0, 5), 0);
        assert_eq!(stirling_number_second_kind(5, 5), 1);
        assert_eq!(stirling_number_second_kind(5, 1), 1);
    }
}
