//! Probability distributions.
//!
//! This module provides implementations of both discrete and continuous probability
//! distributions, including their probability mass/density functions, cumulative
//! distribution functions, and statistical properties.

pub mod continuous;
pub mod discrete;

// Re-export discrete distributions
pub use discrete::{
    bernoulli_cdf, bernoulli_mean, bernoulli_pmf, bernoulli_variance, binomial_cdf, binomial_mean,
    binomial_pmf, binomial_variance, geometric_cdf, geometric_mean, geometric_pmf,
    geometric_variance, poisson_cdf, poisson_mean, poisson_pmf, poisson_variance,
};

// Re-export continuous distributions
pub use continuous::{
    beta_cdf, beta_mean, beta_pdf, beta_variance, chi_squared_cdf, chi_squared_mean,
    chi_squared_pdf, chi_squared_variance, exponential_cdf, exponential_mean, exponential_pdf,
    exponential_variance, f_cdf, f_mean, f_pdf, f_variance, gamma_cdf, gamma_mean, gamma_pdf,
    gamma_variance, normal_cdf, normal_mean, normal_pdf, normal_variance, students_t_cdf,
    students_t_mean, students_t_pdf, students_t_variance, uniform_cdf, uniform_mean, uniform_pdf,
    uniform_variance,
};
