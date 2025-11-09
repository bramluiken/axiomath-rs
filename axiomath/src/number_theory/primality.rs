//! Primality testing and prime number generation.
//!
//! This module provides functions for working with prime numbers:
//! - Primality testing (trial division + Miller-Rabin)
//! - Prime factorization
//! - Sieve of Eratosthenes
//! - Prime counting and nth prime

use num_traits::{One, PrimInt, ToPrimitive, Zero};

/// Tests whether a number is prime using trial division and Miller-Rabin.
///
/// # Algorithm
///
/// For small numbers (n < 1000), uses trial division.
/// For larger numbers, uses the Miller-Rabin primality test with deterministic witnesses.
///
/// # Examples
///
/// ```
/// use axiomath::number_theory::primality::is_prime;
///
/// assert!(is_prime(2_u64));
/// assert!(is_prime(17_u64));
/// assert!(is_prime(97_u64));
/// assert!(!is_prime(1_u64));
/// assert!(!is_prime(4_u64));
/// assert!(!is_prime(100_u64));
/// ```
#[inline]
pub fn is_prime<T>(n: T) -> bool
where
    T: PrimInt + Copy + ToPrimitive,
{
    let n_u64 = n.to_u64().expect("Number too large to convert to u64");

    if n_u64 < 2 {
        return false;
    }
    if n_u64 == 2 || n_u64 == 3 {
        return true;
    }
    if n_u64 % 2 == 0 {
        return false;
    }

    // For small numbers, use trial division
    if n_u64 < 1000 {
        return is_prime_trial_division(n_u64);
    }

    // For larger numbers, use Miller-Rabin
    is_prime_miller_rabin(n_u64)
}

/// Tests primality using trial division up to √n.
#[inline]
fn is_prime_trial_division(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    if n == 2 {
        return true;
    }
    if n % 2 == 0 {
        return false;
    }

    let sqrt_n = (n as f64).sqrt() as u64;
    let mut i = 3;

    while i <= sqrt_n {
        if n % i == 0 {
            return false;
        }
        i += 2;
    }

    true
}

/// Miller-Rabin primality test with deterministic witnesses.
#[inline]
fn is_prime_miller_rabin(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    if n == 2 || n == 3 {
        return true;
    }
    if n % 2 == 0 {
        return false;
    }

    // Write n-1 as 2^r × d
    let mut d = n - 1;
    let mut r = 0;
    while d % 2 == 0 {
        d /= 2;
        r += 1;
    }

    // Deterministic witnesses for numbers up to 2^64
    let witnesses = if n < 2_047 {
        vec![2]
    } else if n < 1_373_653 {
        vec![2, 3]
    } else if n < 9_080_191 {
        vec![31, 73]
    } else if n < 25_326_001 {
        vec![2, 3, 5]
    } else if n < 3_215_031_751 {
        vec![2, 3, 5, 7]
    } else if n < 4_759_123_141 {
        vec![2, 7, 61]
    } else if n < 1_122_004_669_633 {
        vec![2, 13, 23, 1662803]
    } else {
        vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    };

    for &a in &witnesses {
        if a >= n {
            continue;
        }

        let mut x = mod_pow(a, d, n);

        if x == 1 || x == n - 1 {
            continue;
        }

        let mut composite = true;
        for _ in 0..r - 1 {
            x = mod_mul(x, x, n);
            if x == n - 1 {
                composite = false;
                break;
            }
        }

        if composite {
            return false;
        }
    }

    true
}

/// Modular exponentiation: computes (base^exp) mod modulus
#[inline]
fn mod_pow(mut base: u64, mut exp: u64, modulus: u64) -> u64 {
    if modulus == 1 {
        return 0;
    }

    let mut result = 1;
    base %= modulus;

    while exp > 0 {
        if exp % 2 == 1 {
            result = mod_mul(result, base, modulus);
        }
        base = mod_mul(base, base, modulus);
        exp /= 2;
    }

    result
}

/// Modular multiplication: computes (a × b) mod modulus avoiding overflow
#[inline]
fn mod_mul(a: u64, b: u64, modulus: u64) -> u64 {
    ((a as u128 * b as u128) % modulus as u128) as u64
}

/// Computes the prime factorization of a number.
///
/// Returns a vector of (prime, exponent) pairs such that:
/// ```text
/// n = ∏ prime^exponent
/// ```
///
/// # Algorithm
///
/// Uses trial division for small factors, then checks larger primes.
///
/// # Examples
///
/// ```
/// use axiomath::number_theory::primality::prime_factorization;
///
/// assert_eq!(prime_factorization(12_u64), vec![(2, 2), (3, 1)]);
/// assert_eq!(prime_factorization(100_u64), vec![(2, 2), (5, 2)]);
/// assert_eq!(prime_factorization(17_u64), vec![(17, 1)]);
/// ```
#[inline]
pub fn prime_factorization<T>(n: T) -> Vec<(T, u32)>
where
    T: PrimInt + Copy + ToPrimitive,
{
    let mut factors = Vec::new();
    let mut num = n;

    if num.is_zero() || num.is_one() {
        return factors;
    }

    // Factor out 2s
    let two = T::one() + T::one();
    let mut count = 0;
    while num % two == T::zero() {
        num = num / two;
        count += 1;
    }
    if count > 0 {
        factors.push((two, count));
    }

    // Factor out odd numbers
    let mut divisor = two + T::one(); // Start at 3
    while divisor * divisor <= num {
        count = 0;
        while num % divisor == T::zero() {
            num = num / divisor;
            count += 1;
        }
        if count > 0 {
            factors.push((divisor, count));
        }
        divisor = divisor + two; // Check only odd divisors
    }

    // If num > 1, then it's a prime factor
    if num > T::one() {
        factors.push((num, 1));
    }

    factors
}

/// Generates all primes up to a given limit using the Sieve of Eratosthenes.
///
/// # Algorithm
///
/// The Sieve of Eratosthenes works by iteratively marking multiples of each prime as composite.
///
/// Time complexity: O(n log log n)
/// Space complexity: O(n)
///
/// # Examples
///
/// ```
/// use axiomath::number_theory::primality::sieve_of_eratosthenes;
///
/// let primes = sieve_of_eratosthenes(20);
/// assert_eq!(primes, vec![2, 3, 5, 7, 11, 13, 17, 19]);
///
/// let primes = sieve_of_eratosthenes(10);
/// assert_eq!(primes, vec![2, 3, 5, 7]);
/// ```
#[inline]
pub fn sieve_of_eratosthenes(limit: usize) -> Vec<usize> {
    if limit < 2 {
        return vec![];
    }

    let mut is_prime = vec![true; limit + 1];
    is_prime[0] = false;
    is_prime[1] = false;

    let sqrt_limit = (limit as f64).sqrt() as usize;

    for i in 2..=sqrt_limit {
        if is_prime[i] {
            let mut j = i * i;
            while j <= limit {
                is_prime[j] = false;
                j += i;
            }
        }
    }

    is_prime
        .iter()
        .enumerate()
        .filter_map(|(num, &prime)| if prime { Some(num) } else { None })
        .collect()
}

/// Returns the nth prime number (1-indexed).
///
/// # Examples
///
/// ```
/// use axiomath::number_theory::primality::nth_prime;
///
/// assert_eq!(nth_prime(1), 2);
/// assert_eq!(nth_prime(2), 3);
/// assert_eq!(nth_prime(10), 29);
/// assert_eq!(nth_prime(100), 541);
/// ```
///
/// # Panics
///
/// Panics if n is 0.
#[inline]
pub fn nth_prime(n: usize) -> u64 {
    if n == 0 {
        panic!("n must be at least 1 (1-indexed)");
    }

    if n == 1 {
        return 2;
    }

    // Estimate upper bound using prime number theorem
    // π(x) ≈ x / ln(x), so x ≈ n * ln(n) for the nth prime
    let estimate = if n < 6 {
        15
    } else {
        (n as f64 * (n as f64).ln() * 1.3) as usize
    };

    let primes = sieve_of_eratosthenes(estimate);

    if n <= primes.len() {
        primes[n - 1] as u64
    } else {
        // If estimate was too small, extend the search
        let mut candidate = estimate as u64 + 1;
        let mut count = primes.len();

        while count < n {
            if is_prime(candidate) {
                count += 1;
                if count == n {
                    return candidate;
                }
            }
            candidate += 2; // Check only odd numbers
        }

        candidate
    }
}

/// Counts the number of primes less than or equal to the limit (prime counting function π(x)).
///
/// # Examples
///
/// ```
/// use axiomath::number_theory::primality::prime_count;
///
/// assert_eq!(prime_count(10), 4); // 2, 3, 5, 7
/// assert_eq!(prime_count(20), 8); // 2, 3, 5, 7, 11, 13, 17, 19
/// assert_eq!(prime_count(100), 25);
/// ```
#[inline]
pub fn prime_count(limit: usize) -> usize {
    sieve_of_eratosthenes(limit).len()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_prime_small() {
        assert!(!is_prime(0_u64));
        assert!(!is_prime(1_u64));
        assert!(is_prime(2_u64));
        assert!(is_prime(3_u64));
        assert!(!is_prime(4_u64));
        assert!(is_prime(5_u64));
        assert!(!is_prime(6_u64));
        assert!(is_prime(7_u64));
        assert!(!is_prime(8_u64));
        assert!(!is_prime(9_u64));
        assert!(!is_prime(10_u64));
        assert!(is_prime(11_u64));
    }

    #[test]
    fn test_is_prime_larger() {
        assert!(is_prime(97_u64));
        assert!(is_prime(101_u64));
        assert!(is_prime(1009_u64));
        assert!(is_prime(10007_u64));
        assert!(!is_prime(1000_u64));
        assert!(!is_prime(10000_u64));
    }

    #[test]
    fn test_is_prime_mersenne() {
        // Mersenne primes: 2^p - 1
        assert!(is_prime(31_u64)); // 2^5 - 1
        assert!(is_prime(127_u64)); // 2^7 - 1
        assert!(is_prime(8191_u64)); // 2^13 - 1
    }

    #[test]
    fn test_prime_factorization_basic() {
        assert_eq!(prime_factorization(1_u64), vec![]);
        assert_eq!(prime_factorization(2_u64), vec![(2, 1)]);
        assert_eq!(prime_factorization(4_u64), vec![(2, 2)]);
        assert_eq!(prime_factorization(12_u64), vec![(2, 2), (3, 1)]);
        assert_eq!(prime_factorization(100_u64), vec![(2, 2), (5, 2)]);
    }

    #[test]
    fn test_prime_factorization_primes() {
        // Prime numbers should have themselves as the only factor
        assert_eq!(prime_factorization(17_u64), vec![(17, 1)]);
        assert_eq!(prime_factorization(97_u64), vec![(97, 1)]);
    }

    #[test]
    fn test_prime_factorization_product() {
        // Verify that the product of factors equals the original number
        for n in 2..100_u64 {
            let factors = prime_factorization(n);
            let product: u64 = factors
                .iter()
                .map(|&(prime, exp)| prime.pow(exp))
                .product();
            assert_eq!(product, n);
        }
    }

    #[test]
    fn test_sieve_of_eratosthenes() {
        let primes = sieve_of_eratosthenes(20);
        assert_eq!(primes, vec![2, 3, 5, 7, 11, 13, 17, 19]);

        let primes = sieve_of_eratosthenes(10);
        assert_eq!(primes, vec![2, 3, 5, 7]);

        let primes = sieve_of_eratosthenes(2);
        assert_eq!(primes, vec![2]);

        let primes = sieve_of_eratosthenes(1);
        assert_eq!(primes, vec![]);
    }

    #[test]
    fn test_nth_prime() {
        assert_eq!(nth_prime(1), 2);
        assert_eq!(nth_prime(2), 3);
        assert_eq!(nth_prime(3), 5);
        assert_eq!(nth_prime(4), 7);
        assert_eq!(nth_prime(10), 29);
        assert_eq!(nth_prime(25), 97);
    }

    #[test]
    fn test_nth_prime_larger() {
        assert_eq!(nth_prime(100), 541);
        assert_eq!(nth_prime(1000), 7919);
    }

    #[test]
    #[should_panic(expected = "n must be at least 1")]
    fn test_nth_prime_zero() {
        nth_prime(0);
    }

    #[test]
    fn test_prime_count() {
        assert_eq!(prime_count(10), 4);
        assert_eq!(prime_count(20), 8);
        assert_eq!(prime_count(100), 25);
    }

    #[test]
    fn test_prime_count_consistency() {
        // Verify that nth_prime and prime_count are consistent
        for n in 1..50 {
            let p = nth_prime(n);
            let count = prime_count(p as usize);
            assert!(count >= n);
        }
    }
}
