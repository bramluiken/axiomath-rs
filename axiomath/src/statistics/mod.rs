//! Statistical analysis and probability distributions.
//!
//! This module provides comprehensive statistical functionality including descriptive
//! statistics, probability distributions (discrete and continuous), and statistical
//! inference tools.
//!
//! # Module Organization
//!
//! - [`descriptive`]: Descriptive statistics (mean, variance, skewness, correlation, etc.)
//! - [`distributions`]: Probability distributions (Bernoulli, Binomial, Normal, etc.)
//!
//! # Examples
//!
//! ## Descriptive Statistics
//!
//! ```
//! use axiomath::statistics::descriptive::{mean, standard_deviation, correlation};
//!
//! let data = [1.0, 2.0, 3.0, 4.0, 5.0];
//! let avg = mean(&data).unwrap();
//! let std = standard_deviation(&data, true).unwrap();
//!
//! assert!((avg - 3.0).abs() < 1e-10);
//! ```
//!
//! ## Probability Distributions
//!
//! ```
//! use axiomath::statistics::distributions::{
//!     binomial_pmf, normal_pdf, exponential_pdf
//! };
//!
//! // Discrete: Binomial distribution
//! let prob = binomial_pmf(3, 5, 0.5).unwrap();
//!
//! // Continuous: Normal distribution
//! let density = normal_pdf(0.0, 0.0, 1.0).unwrap();
//!
//! // Continuous: Exponential distribution
//! let density = exponential_pdf(1.0, 2.0).unwrap();
//! ```
//!
//! # Statistical Theory
//!
//! This module implements statistical methods from first principles:
//!
//! - **Descriptive Statistics**: Sample and population statistics
//! - **Discrete Distributions**: Bernoulli, Binomial, Poisson, Geometric
//! - **Continuous Distributions**: Normal, Uniform, Exponential, Gamma, Beta,
//!   Chi-squared, Student's t, F-distribution
//!
//! All distributions provide:
//! - Probability mass/density functions (PMF/PDF)
//! - Cumulative distribution functions (CDF)
//! - Mean and variance
//!
//! # Dependencies
//!
//! The statistics module builds on:
//! - [`crate::special`]: Special functions (gamma, beta, error function)
//! - [`crate::number_theory`]: Combinatorics (factorial, binomial coefficient)
//! - [`crate::scalar`]: Elementary functions (exp, ln, sqrt, power)

pub mod descriptive;
pub mod distributions;

// Re-export error type
pub use descriptive::StatisticsError;
