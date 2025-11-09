//! Number theory functions.
//!
//! This module implements number theoretic functions including
//! primality testing, factorization, and modular arithmetic.

pub mod combinatorics;
pub mod divisibility;
pub mod modular;
pub mod primality;

// Re-export all public functions for convenient access
pub use combinatorics::{
    binomial_coefficient, factorial, permutations, stirling_number_first_kind,
    stirling_number_second_kind,
};
pub use divisibility::{extended_gcd, greatest_common_divisor, least_common_multiple};
pub use modular::{modular_add, modular_inverse, modular_multiply, modular_power};
pub use primality::{
    is_prime, nth_prime, prime_count, prime_factorization, sieve_of_eratosthenes,
};
