//! Discrete probability distributions.
//!
//! This module provides implementations of common discrete probability distributions,
//! including their probability mass functions (PMF), cumulative distribution functions (CDF),
//! and statistical properties.
//!
//! # Distributions
//!
//! - **Bernoulli**: Models a single trial with two possible outcomes (success/failure)
//! - **Binomial**: Number of successes in n independent Bernoulli trials
//! - **Poisson**: Number of events occurring in a fixed interval of time/space
//! - **Geometric**: Number of trials until the first success
//!
//! # Examples
//!
//! ```
//! use axiomath::statistics::distributions::discrete::{
//!     bernoulli_pmf, binomial_pmf, poisson_pmf, geometric_pmf
//! };
//!
//! // Bernoulli: probability of success
//! let prob = bernoulli_pmf(1, 0.6); // P(X=1) with p=0.6
//! assert!((prob - 0.6).abs() < 1e-10);
//!
//! // Binomial: 3 successes in 5 trials with p=0.5
//! let prob = binomial_pmf(3, 5, 0.5);
//! assert!((prob - 0.3125).abs() < 1e-10);
//!
//! // Poisson: 2 events with rate λ=1.5
//! let prob = poisson_pmf(2, 1.5);
//! assert!((prob - 0.2510).abs() < 1e-3);
//!
//! // Geometric: first success on 3rd trial with p=0.5
//! let prob = geometric_pmf(3, 0.5);
//! assert!((prob - 0.125).abs() < 1e-10);
//! ```

use crate::core::traits::Float;
use crate::number_theory::{binomial_coefficient, factorial};
use crate::scalar::exponential;
use crate::statistics::StatisticsError;

// ============================================================================
// Bernoulli Distribution
// ============================================================================

/// Bernoulli probability mass function.
///
/// Computes P(X = k) for a Bernoulli distribution with success probability p.
///
/// # Mathematical Definition
///
/// P(X = k) = p^k * (1-p)^(1-k) for k ∈ {0, 1}
///
/// # Parameters
///
/// - `k`: The outcome (0 for failure, 1 for success)
/// - `p`: Probability of success (must be in [0, 1])
///
/// # Returns
///
/// The probability of observing k, or an error if parameters are invalid.
///
/// # Examples
///
/// ```
/// use axiomath::statistics::distributions::discrete::bernoulli_pmf;
///
/// let prob_success = bernoulli_pmf(1, 0.7).unwrap();
/// assert!((prob_success - 0.7).abs() < 1e-10);
///
/// let prob_failure = bernoulli_pmf(0, 0.7).unwrap();
/// assert!((prob_failure - 0.3).abs() < 1e-10);
/// ```
pub fn bernoulli_pmf<T: Float>(k: i32, p: T) -> Result<T, StatisticsError> {
    if p < T::zero() || p > T::one() {
        return Err(StatisticsError::InvalidProbability);
    }

    match k {
        0 => Ok(T::one() - p),
        1 => Ok(p),
        _ => Ok(T::zero()),
    }
}

/// Bernoulli cumulative distribution function.
///
/// Computes P(X ≤ k) for a Bernoulli distribution.
///
/// # Examples
///
/// ```
/// use axiomath::statistics::distributions::discrete::bernoulli_cdf;
///
/// let cdf = bernoulli_cdf(0, 0.7).unwrap();
/// assert!((cdf - 0.3).abs() < 1e-10);
/// ```
pub fn bernoulli_cdf<T: Float>(k: i32, p: T) -> Result<T, StatisticsError> {
    if p < T::zero() || p > T::one() {
        return Err(StatisticsError::InvalidProbability);
    }

    if k < 0 {
        Ok(T::zero())
    } else if k < 1 {
        Ok(T::one() - p)
    } else {
        Ok(T::one())
    }
}

/// Bernoulli distribution mean.
///
/// E[X] = p
pub fn bernoulli_mean<T: Float>(p: T) -> Result<T, StatisticsError> {
    if p < T::zero() || p > T::one() {
        return Err(StatisticsError::InvalidProbability);
    }
    Ok(p)
}

/// Bernoulli distribution variance.
///
/// Var(X) = p(1-p)
pub fn bernoulli_variance<T: Float>(p: T) -> Result<T, StatisticsError> {
    if p < T::zero() || p > T::one() {
        return Err(StatisticsError::InvalidProbability);
    }
    Ok(p * (T::one() - p))
}

// ============================================================================
// Binomial Distribution
// ============================================================================

/// Binomial probability mass function.
///
/// Computes P(X = k) for a binomial distribution with n trials and success probability p.
///
/// # Mathematical Definition
///
/// P(X = k) = C(n,k) * p^k * (1-p)^(n-k)
///
/// where C(n,k) is the binomial coefficient "n choose k".
///
/// # Parameters
///
/// - `k`: Number of successes
/// - `n`: Number of trials
/// - `p`: Probability of success per trial (must be in [0, 1])
///
/// # Complexity
///
/// O(min(k, n-k)) for computing the binomial coefficient
///
/// # Examples
///
/// ```
/// use axiomath::statistics::distributions::discrete::binomial_pmf;
///
/// // Probability of 3 heads in 5 coin flips
/// let prob = binomial_pmf(3, 5, 0.5).unwrap();
/// assert!((prob - 0.3125).abs() < 1e-10);
/// ```
pub fn binomial_pmf<T: Float>(k: u64, n: u64, p: T) -> Result<T, StatisticsError> {
    if p < T::zero() || p > T::one() {
        return Err(StatisticsError::InvalidProbability);
    }

    if k > n {
        return Ok(T::zero());
    }

    let binom_coeff = T::from(binomial_coefficient(n, k)).unwrap();
    let prob_success = p.powi(k as i32);
    let prob_failure = (T::one() - p).powi((n - k) as i32);

    Ok(binom_coeff * prob_success * prob_failure)
}

/// Binomial cumulative distribution function.
///
/// Computes P(X ≤ k) for a binomial distribution.
///
/// # Examples
///
/// ```
/// use axiomath::statistics::distributions::discrete::binomial_cdf;
///
/// let cdf = binomial_cdf(2, 5, 0.5).unwrap();
/// assert!(cdf > 0.5); // More than half the probability mass is at k ≤ 2
/// ```
pub fn binomial_cdf<T: Float>(k: u64, n: u64, p: T) -> Result<T, StatisticsError> {
    if p < T::zero() || p > T::one() {
        return Err(StatisticsError::InvalidProbability);
    }

    let mut cdf = T::zero();
    for i in 0..=k.min(n) {
        cdf = cdf + binomial_pmf(i, n, p)?;
    }

    Ok(cdf)
}

/// Binomial distribution mean.
///
/// E[X] = np
pub fn binomial_mean<T: Float>(n: u64, p: T) -> Result<T, StatisticsError> {
    if p < T::zero() || p > T::one() {
        return Err(StatisticsError::InvalidProbability);
    }
    Ok(T::from(n).unwrap() * p)
}

/// Binomial distribution variance.
///
/// Var(X) = np(1-p)
pub fn binomial_variance<T: Float>(n: u64, p: T) -> Result<T, StatisticsError> {
    if p < T::zero() || p > T::one() {
        return Err(StatisticsError::InvalidProbability);
    }
    Ok(T::from(n).unwrap() * p * (T::one() - p))
}

// ============================================================================
// Poisson Distribution
// ============================================================================

/// Poisson probability mass function.
///
/// Computes P(X = k) for a Poisson distribution with rate parameter λ.
///
/// # Mathematical Definition
///
/// P(X = k) = (λ^k * e^(-λ)) / k!
///
/// # Parameters
///
/// - `k`: Number of events
/// - `lambda`: Rate parameter (must be positive)
///
/// # Applications
///
/// - Number of arrivals in a queue
/// - Number of radioactive decays
/// - Number of phone calls received
/// - Rare events in a fixed interval
///
/// # Examples
///
/// ```
/// use axiomath::statistics::distributions::discrete::poisson_pmf;
///
/// // Probability of 3 events with rate 2.5
/// let prob = poisson_pmf(3, 2.5).unwrap();
/// assert!((prob - 0.2138).abs() < 1e-3);
/// ```
pub fn poisson_pmf<T: Float>(k: u64, lambda: T) -> Result<T, StatisticsError> {
    if lambda <= T::zero() {
        return Err(StatisticsError::InvalidParameter(
            "Lambda must be positive".to_string(),
        ));
    }

    // P(X = k) = (λ^k * e^(-λ)) / k!
    let lambda_k = lambda.powi(k as i32);
    let exp_neg_lambda = exponential(-lambda);
    let k_factorial = T::from(factorial(k)).unwrap();

    Ok(lambda_k * exp_neg_lambda / k_factorial)
}

/// Poisson cumulative distribution function.
///
/// Computes P(X ≤ k) for a Poisson distribution.
///
/// # Examples
///
/// ```
/// use axiomath::statistics::distributions::discrete::poisson_cdf;
///
/// let cdf = poisson_cdf(5, 3.0).unwrap();
/// assert!(cdf > 0.9); // Most probability mass is at k ≤ 5 when λ=3
/// ```
pub fn poisson_cdf<T: Float>(k: u64, lambda: T) -> Result<T, StatisticsError> {
    if lambda <= T::zero() {
        return Err(StatisticsError::InvalidParameter(
            "Lambda must be positive".to_string(),
        ));
    }

    let mut cdf = T::zero();
    for i in 0..=k {
        cdf = cdf + poisson_pmf(i, lambda)?;
    }

    Ok(cdf)
}

/// Poisson distribution mean.
///
/// E[X] = λ
pub fn poisson_mean<T: Float>(lambda: T) -> Result<T, StatisticsError> {
    if lambda <= T::zero() {
        return Err(StatisticsError::InvalidParameter(
            "Lambda must be positive".to_string(),
        ));
    }
    Ok(lambda)
}

/// Poisson distribution variance.
///
/// Var(X) = λ
pub fn poisson_variance<T: Float>(lambda: T) -> Result<T, StatisticsError> {
    if lambda <= T::zero() {
        return Err(StatisticsError::InvalidParameter(
            "Lambda must be positive".to_string(),
        ));
    }
    Ok(lambda)
}

// ============================================================================
// Geometric Distribution
// ============================================================================

/// Geometric probability mass function.
///
/// Computes P(X = k) for a geometric distribution with success probability p.
/// This is the probability that the first success occurs on the k-th trial.
///
/// # Mathematical Definition
///
/// P(X = k) = (1-p)^(k-1) * p for k ≥ 1
///
/// # Parameters
///
/// - `k`: Trial number of first success (must be ≥ 1)
/// - `p`: Probability of success per trial (must be in (0, 1])
///
/// # Applications
///
/// - Number of attempts until first success
/// - Waiting time until first occurrence
/// - Memoryless processes
///
/// # Examples
///
/// ```
/// use axiomath::statistics::distributions::discrete::geometric_pmf;
///
/// // Probability first success is on 3rd trial with p=0.5
/// let prob = geometric_pmf(3, 0.5).unwrap();
/// assert!((prob - 0.125).abs() < 1e-10);
/// ```
pub fn geometric_pmf<T: Float>(k: u64, p: T) -> Result<T, StatisticsError> {
    if p <= T::zero() || p > T::one() {
        return Err(StatisticsError::InvalidProbability);
    }

    if k < 1 {
        return Err(StatisticsError::InvalidParameter(
            "k must be at least 1".to_string(),
        ));
    }

    // P(X = k) = (1-p)^(k-1) * p
    let prob_failure = (T::one() - p).powi((k - 1) as i32);
    Ok(prob_failure * p)
}

/// Geometric cumulative distribution function.
///
/// Computes P(X ≤ k) for a geometric distribution.
///
/// # Examples
///
/// ```
/// use axiomath::statistics::distributions::discrete::geometric_cdf;
///
/// let cdf = geometric_cdf(3, 0.5).unwrap();
/// assert!((cdf - 0.875).abs() < 1e-10);
/// ```
pub fn geometric_cdf<T: Float>(k: u64, p: T) -> Result<T, StatisticsError> {
    if p <= T::zero() || p > T::one() {
        return Err(StatisticsError::InvalidProbability);
    }

    if k < 1 {
        return Ok(T::zero());
    }

    // P(X ≤ k) = 1 - (1-p)^k
    Ok(T::one() - (T::one() - p).powi(k as i32))
}

/// Geometric distribution mean.
///
/// E[X] = 1/p
pub fn geometric_mean<T: Float>(p: T) -> Result<T, StatisticsError> {
    if p <= T::zero() || p > T::one() {
        return Err(StatisticsError::InvalidProbability);
    }
    Ok(T::one() / p)
}

/// Geometric distribution variance.
///
/// Var(X) = (1-p)/p²
pub fn geometric_variance<T: Float>(p: T) -> Result<T, StatisticsError> {
    if p <= T::zero() || p > T::one() {
        return Err(StatisticsError::InvalidProbability);
    }
    Ok((T::one() - p) / (p * p))
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    // Bernoulli tests
    #[test]
    fn test_bernoulli_pmf() {
        let p = 0.7;
        assert_abs_diff_eq!(bernoulli_pmf(1, p).unwrap(), 0.7, epsilon = 1e-10);
        assert_abs_diff_eq!(bernoulli_pmf(0, p).unwrap(), 0.3, epsilon = 1e-10);
        assert_eq!(bernoulli_pmf(2, p).unwrap(), 0.0);
    }

    #[test]
    fn test_bernoulli_cdf() {
        let p = 0.6;
        assert_eq!(bernoulli_cdf(-1, p).unwrap(), 0.0);
        assert_abs_diff_eq!(bernoulli_cdf(0, p).unwrap(), 0.4, epsilon = 1e-10);
        assert_eq!(bernoulli_cdf(1, p).unwrap(), 1.0);
        assert_eq!(bernoulli_cdf(2, p).unwrap(), 1.0);
    }

    #[test]
    fn test_bernoulli_moments() {
        let p = 0.6;
        assert_abs_diff_eq!(bernoulli_mean(p).unwrap(), 0.6, epsilon = 1e-10);
        assert_abs_diff_eq!(bernoulli_variance(p).unwrap(), 0.24, epsilon = 1e-10);
    }

    // Binomial tests
    #[test]
    fn test_binomial_pmf() {
        // 3 successes in 5 trials with p=0.5
        let prob = binomial_pmf(3, 5, 0.5).unwrap();
        assert_abs_diff_eq!(prob, 0.3125, epsilon = 1e-10);

        // Edge case: 0 successes
        let prob = binomial_pmf(0, 5, 0.5).unwrap();
        assert_abs_diff_eq!(prob, 0.03125, epsilon = 1e-10);

        // Edge case: all successes
        let prob = binomial_pmf(5, 5, 0.5).unwrap();
        assert_abs_diff_eq!(prob, 0.03125, epsilon = 1e-10);
    }

    #[test]
    fn test_binomial_cdf() {
        let cdf = binomial_cdf(2, 5, 0.5).unwrap();
        // P(X ≤ 2) = P(0) + P(1) + P(2)
        let expected = 0.03125 + 0.15625 + 0.3125;
        assert_abs_diff_eq!(cdf, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_binomial_moments() {
        let n = 10;
        let p = 0.3;
        assert_abs_diff_eq!(binomial_mean(n, p).unwrap(), 3.0, epsilon = 1e-10);
        assert_abs_diff_eq!(binomial_variance(n, p).unwrap(), 2.1, epsilon = 1e-10);
    }

    #[test]
    fn test_binomial_k_exceeds_n() {
        let prob = binomial_pmf(10, 5, 0.5).unwrap();
        assert_eq!(prob, 0.0);
    }

    // Poisson tests
    #[test]
    fn test_poisson_pmf() {
        // λ = 2, k = 3: P(X=3) = 2³ * e^(-2) / 3!
        let prob = poisson_pmf(3, 2.0).unwrap();
        let expected = 8.0 * (-2.0_f64).exp() / 6.0;
        assert_abs_diff_eq!(prob, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_poisson_pmf_small_lambda() {
        let prob = poisson_pmf(0, 0.5).unwrap();
        let expected = (-0.5_f64).exp();
        assert_abs_diff_eq!(prob, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_poisson_moments() {
        let lambda = 3.5;
        assert_abs_diff_eq!(poisson_mean(lambda).unwrap(), 3.5, epsilon = 1e-10);
        assert_abs_diff_eq!(poisson_variance(lambda).unwrap(), 3.5, epsilon = 1e-10);
    }

    #[test]
    fn test_poisson_cdf() {
        let cdf = poisson_cdf(2, 1.5).unwrap();
        // Sum of P(0), P(1), P(2)
        let p0 = poisson_pmf(0, 1.5).unwrap();
        let p1 = poisson_pmf(1, 1.5).unwrap();
        let p2 = poisson_pmf(2, 1.5).unwrap();
        assert_abs_diff_eq!(cdf, p0 + p1 + p2, epsilon = 1e-10);
    }

    // Geometric tests
    #[test]
    fn test_geometric_pmf() {
        let p = 0.5;
        // First success on trial 1: p
        assert_abs_diff_eq!(geometric_pmf(1, p).unwrap(), 0.5, epsilon = 1e-10);
        // First success on trial 2: (1-p)*p
        assert_abs_diff_eq!(geometric_pmf(2, p).unwrap(), 0.25, epsilon = 1e-10);
        // First success on trial 3: (1-p)²*p
        assert_abs_diff_eq!(geometric_pmf(3, p).unwrap(), 0.125, epsilon = 1e-10);
    }

    #[test]
    fn test_geometric_cdf() {
        let p = 0.5;
        // P(X ≤ 3) = 1 - (1-p)³
        let cdf = geometric_cdf(3, p).unwrap();
        assert_abs_diff_eq!(cdf, 0.875, epsilon = 1e-10);
    }

    #[test]
    fn test_geometric_moments() {
        let p = 0.25;
        // E[X] = 1/p = 4
        assert_abs_diff_eq!(geometric_mean(p).unwrap(), 4.0, epsilon = 1e-10);
        // Var(X) = (1-p)/p² = 0.75/0.0625 = 12
        assert_abs_diff_eq!(geometric_variance(p).unwrap(), 12.0, epsilon = 1e-10);
    }

    #[test]
    fn test_geometric_invalid_k() {
        assert!(geometric_pmf(0, 0.5).is_err());
    }

    // Error handling tests
    #[test]
    fn test_invalid_probability() {
        assert!(bernoulli_pmf(1, -0.1).is_err());
        assert!(bernoulli_pmf(1, 1.5).is_err());
        assert!(binomial_pmf(3, 5, -0.1).is_err());
        assert!(binomial_pmf(3, 5, 1.5).is_err());
        assert!(geometric_pmf(1, 0.0).is_err());
        assert!(geometric_pmf(1, 1.5).is_err());
    }

    #[test]
    fn test_invalid_lambda() {
        assert!(poisson_pmf(3, 0.0).is_err());
        assert!(poisson_pmf(3, -1.5).is_err());
    }

    // Type tests
    #[test]
    fn test_f32_support() {
        let prob: f32 = binomial_pmf(3, 5, 0.5_f32).unwrap();
        assert_abs_diff_eq!(prob, 0.3125_f32, epsilon = 1e-6);
    }
}
