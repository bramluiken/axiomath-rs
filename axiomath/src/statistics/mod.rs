//! Statistical analysis and probability distributions.
//!
//! This module provides comprehensive statistical functionality including descriptive
//! statistics, probability distributions (discrete and continuous), statistical
//! inference, and information theory.
//!
//! # Module Organization
//!
//! - [`descriptive`]: Descriptive statistics (mean, variance, skewness, correlation, etc.)
//! - [`distributions`]: Probability distributions (Bernoulli, Binomial, Normal, etc.)
//! - [`inference`]: Statistical inference (hypothesis testing, confidence intervals)
//! - [`information`]: Information theory (entropy, mutual information, divergence)
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
//! ## Statistical Inference
//!
//! ```
//! use axiomath::statistics::inference::{t_test_one_sample, mean_confidence_interval};
//!
//! let data = [2.1, 2.3, 1.9, 2.0, 2.2];
//!
//! // Hypothesis test
//! let result = t_test_one_sample(&data, 2.0).unwrap();
//! // Check if p_value < 0.05 for significance
//!
//! // Confidence interval
//! let (lower, upper) = mean_confidence_interval(&data, 0.95).unwrap();
//! ```
//!
//! ## Information Theory
//!
//! ```
//! use axiomath::statistics::information::{shannon_entropy, kl_divergence};
//!
//! // Shannon entropy
//! let probs = [0.5, 0.5];
//! let h = shannon_entropy(&probs).unwrap();
//!
//! // KL divergence
//! let p = [0.5, 0.5];
//! let q = [0.6, 0.4];
//! let div = kl_divergence(&p, &q).unwrap();
//! ```
//!
//! # Statistical Theory
//!
//! This module implements statistical methods from first principles:
//!
//! - **Descriptive Statistics**: Sample and population statistics, correlations, rank-based methods
//! - **Discrete Distributions**: Bernoulli, Binomial, Poisson, Geometric
//! - **Continuous Distributions**: Normal, Uniform, Exponential, Gamma, Beta,
//!   Chi-squared, Student's t, F-distribution
//! - **Statistical Inference**: t-tests, z-tests, chi-squared tests, ANOVA, confidence intervals
//! - **Information Theory**: Entropy measures, mutual information, KL divergence, cross-entropy
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
pub mod inference;
pub mod information;

// Re-export error type
pub use descriptive::StatisticsError;
