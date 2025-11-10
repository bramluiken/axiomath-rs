//! Descriptive statistics for data analysis.
//!
//! This module provides functions for computing summary statistics that describe
//! the central tendency, dispersion, and shape of data distributions.
//!
//! # Categories
//!
//! - **Central Tendency**: mean, median, mode, geometric mean, harmonic mean
//! - **Dispersion**: variance, standard deviation, range, MAD, IQR
//! - **Shape**: skewness, kurtosis
//! - **Order Statistics**: percentiles, quantiles, quartiles
//! - **Multivariate**: covariance, correlation
//!
//! # Examples
//!
//! ```
//! use axiomath::statistics::descriptive::{mean, standard_deviation, median};
//!
//! let data = [1.0, 2.0, 3.0, 4.0, 5.0];
//!
//! let m = mean(&data).unwrap();
//! assert!((m - 3.0).abs() < 1e-10);
//!
//! let sd = standard_deviation(&data, true).unwrap();
//! assert!(sd > 0.0);
//! ```
//!
//! # References
//!
//! - Hogg, R. V., & Tanis, E. A. (2009). "Probability and Statistical Inference"
//! - DeGroot, M. H., & Schervish, M. J. (2012). "Probability and Statistics"

use crate::core::Float;
use num_traits::Float as NumFloat;
use std::collections::HashMap;
use std::hash::Hash;

/// Defines errors that can occur in statistical computations.
#[derive(Debug, Clone, PartialEq)]
pub enum StatisticsError {
    /// The input data is empty.
    EmptyData,

    /// Insufficient data for the requested operation.
    InsufficientData {
        /// Number of data points required.
        required: usize,
        /// Number of data points provided.
        provided: usize,
    },

    /// A probability value is outside the valid range [0, 1].
    InvalidProbability,

    /// A parameter has an invalid value.
    InvalidParameter(String),

    /// Data contains invalid values (NaN or infinite).
    InvalidData,
}

impl std::fmt::Display for StatisticsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            StatisticsError::EmptyData => write!(f, "Input data is empty"),
            StatisticsError::InsufficientData { required, provided } => {
                write!(f, "Insufficient data: need {} points, got {}", required, provided)
            }
            StatisticsError::InvalidProbability => {
                write!(f, "Probability must be in range [0, 1]")
            }
            StatisticsError::InvalidParameter(msg) => write!(f, "Invalid parameter: {}", msg),
            StatisticsError::InvalidData => write!(f, "Data contains NaN or infinite values"),
        }
    }
}

impl std::error::Error for StatisticsError {}

// ============================================================================
// CENTRAL TENDENCY
// ============================================================================

/// Computes the arithmetic mean (average) of a dataset.
///
/// The mean is defined as:
/// ```text
/// μ = (1/n) ∑ᵢ xᵢ
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::mean;
///
/// let data = [1.0, 2.0, 3.0, 4.0, 5.0];
/// let m = mean(&data).unwrap();
/// assert!((m - 3.0).abs() < 1e-10);
/// ```
///
/// # Errors
///
/// Returns `StatisticsError::EmptyData` if the input is empty.
pub fn mean<T: Float>(data: &[T]) -> Result<T, StatisticsError> {
    if data.is_empty() {
        return Err(StatisticsError::EmptyData);
    }

    let sum: T = data.iter().copied().fold(T::zero(), |acc, x| acc + x);
    let n = T::from(data.len()).unwrap();
    Ok(sum / n)
}

/// Computes the weighted mean of a dataset.
///
/// The weighted mean is defined as:
/// ```text
/// μ_w = (∑ᵢ wᵢxᵢ) / (∑ᵢ wᵢ)
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::weighted_mean;
///
/// let data = [1.0, 2.0, 3.0];
/// let weights = [1.0, 2.0, 1.0];
/// let m = weighted_mean(&data, &weights).unwrap();
/// assert!((m - 2.0).abs() < 1e-10);
/// ```
///
/// # Errors
///
/// - `EmptyData` if data is empty
/// - `InvalidParameter` if data and weights have different lengths
pub fn weighted_mean<T: Float>(data: &[T], weights: &[T]) -> Result<T, StatisticsError> {
    if data.is_empty() {
        return Err(StatisticsError::EmptyData);
    }
    if data.len() != weights.len() {
        return Err(StatisticsError::InvalidParameter(
            "Data and weights must have the same length".to_string(),
        ));
    }

    let weighted_sum: T = data.iter().zip(weights.iter())
        .map(|(&x, &w)| x * w)
        .fold(T::zero(), |acc, x| acc + x);

    let weight_sum: T = weights.iter().copied().fold(T::zero(), |acc, w| acc + w);

    if weight_sum.is_zero() {
        return Err(StatisticsError::InvalidParameter("Sum of weights cannot be zero".to_string()));
    }

    Ok(weighted_sum / weight_sum)
}

/// Computes the median (50th percentile) of a dataset.
///
/// The median is the middle value when data is sorted. For even-sized datasets,
/// it's the average of the two middle values.
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::median;
///
/// let data = vec![1.0, 3.0, 2.0, 5.0, 4.0];
/// let m = median(&data).unwrap();
/// assert!((m - 3.0).abs() < 1e-10);
/// ```
///
/// # Note
///
/// This function sorts the data internally, so the input is consumed and modified.
pub fn median<T: Float>(data: &mut [T]) -> Result<T, StatisticsError> {
    if data.is_empty() {
        return Err(StatisticsError::EmptyData);
    }

    // Sort the data
    data.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let n = data.len();
    if n % 2 == 1 {
        // Odd number of elements
        Ok(data[n / 2])
    } else {
        // Even number of elements - average the two middle values
        let mid1 = data[n / 2 - 1];
        let mid2 = data[n / 2];
        Ok((mid1 + mid2) / T::from(2.0).unwrap())
    }
}

/// Computes the mode(s) of a dataset.
///
/// The mode is the most frequently occurring value. Multiple modes can exist.
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::mode;
///
/// let data = vec![1, 2, 2, 3, 3, 3, 4];
/// let modes = mode(&data).unwrap();
/// assert_eq!(modes, vec![3]);
/// ```
pub fn mode<T: Eq + Hash + Copy>(data: &[T]) -> Result<Vec<T>, StatisticsError> {
    if data.is_empty() {
        return Err(StatisticsError::EmptyData);
    }

    let mut counts = HashMap::new();
    for &value in data {
        *counts.entry(value).or_insert(0) += 1;
    }

    let max_count = counts.values().max().unwrap();
    let modes: Vec<T> = counts.iter()
        .filter(|(_, &count)| count == *max_count)
        .map(|(&value, _)| value)
        .collect();

    Ok(modes)
}

/// Computes the geometric mean of a dataset.
///
/// The geometric mean is defined as:
/// ```text
/// GM = (∏ᵢ xᵢ)^(1/n)
/// ```
///
/// Computed as `exp(mean(log(x)))` for numerical stability.
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::geometric_mean;
///
/// let data = [1.0, 2.0, 4.0, 8.0];
/// let gm = geometric_mean(&data).unwrap();
/// assert!((gm - (1.0 * 2.0 * 4.0 * 8.0_f64).powf(0.25)).abs() < 1e-10);
/// ```
///
/// # Errors
///
/// Returns error if any value is non-positive.
pub fn geometric_mean<T: Float>(data: &[T]) -> Result<T, StatisticsError> {
    if data.is_empty() {
        return Err(StatisticsError::EmptyData);
    }

    // Check for non-positive values
    for &x in data {
        if x <= T::zero() {
            return Err(StatisticsError::InvalidParameter(
                "Geometric mean requires positive values".to_string(),
            ));
        }
    }

    let log_sum: T = data.iter().map(|&x| x.ln()).fold(T::zero(), |acc, x| acc + x);
    let n = T::from(data.len()).unwrap();
    Ok((log_sum / n).exp())
}

// ============================================================================
// DISPERSION
// ============================================================================

/// Computes the variance of a dataset.
///
/// Variance measures how spread out the data is:
/// ```text
/// σ² = (1/n) ∑ᵢ (xᵢ - μ)²     (population variance, sample=false)
/// s² = (1/(n-1)) ∑ᵢ (xᵢ - μ)²  (sample variance, sample=true)
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::variance;
///
/// let data = [1.0, 2.0, 3.0, 4.0, 5.0];
/// let var = variance(&data, true).unwrap();  // Sample variance
/// assert!(var > 0.0);
/// ```
///
/// # Parameters
///
/// - `sample`: If true, computes sample variance (divides by n-1). If false, population variance (divides by n).
pub fn variance<T: Float>(data: &[T], sample: bool) -> Result<T, StatisticsError> {
    if data.is_empty() {
        return Err(StatisticsError::EmptyData);
    }
    if sample && data.len() < 2 {
        return Err(StatisticsError::InsufficientData {
            required: 2,
            provided: data.len(),
        });
    }

    let m = mean(data)?;
    let sum_sq_diff: T = data.iter()
        .map(|&x| {
            let diff = x - m;
            diff * diff
        })
        .fold(T::zero(), |acc, x| acc + x);

    let divisor = if sample {
        T::from(data.len() - 1).unwrap()
    } else {
        T::from(data.len()).unwrap()
    };

    Ok(sum_sq_diff / divisor)
}

/// Computes the standard deviation of a dataset.
///
/// Standard deviation is the square root of variance:
/// ```text
/// σ = √(variance)
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::standard_deviation;
///
/// let data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
/// let sd = standard_deviation(&data, true).unwrap();
/// assert!((sd - 2.0).abs() < 0.1);
/// ```
pub fn standard_deviation<T: Float>(data: &[T], sample: bool) -> Result<T, StatisticsError> {
    Ok(variance(data, sample)?.sqrt())
}

/// Computes the mean absolute deviation (MAD) from the mean.
///
/// ```text
/// MAD = (1/n) ∑ᵢ |xᵢ - μ|
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::mean_absolute_deviation;
///
/// let data = [1.0, 2.0, 3.0, 4.0, 5.0];
/// let mad = mean_absolute_deviation(&data).unwrap();
/// assert!(mad > 0.0);
/// ```
pub fn mean_absolute_deviation<T: Float>(data: &[T]) -> Result<T, StatisticsError> {
    if data.is_empty() {
        return Err(StatisticsError::EmptyData);
    }

    let m = mean(data)?;
    let sum_abs_diff: T = data.iter()
        .map(|&x| (x - m).abs())
        .fold(T::zero(), |acc, x| acc + x);

    let n = T::from(data.len()).unwrap();
    Ok(sum_abs_diff / n)
}

/// Computes the range (max - min) of a dataset.
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::range;
///
/// let data = [1.0, 2.0, 5.0, 3.0];
/// let r = range(&data).unwrap();
/// assert!((r - 4.0).abs() < 1e-10);
/// ```
pub fn range<T: Float>(data: &[T]) -> Result<T, StatisticsError> {
    if data.is_empty() {
        return Err(StatisticsError::EmptyData);
    }

    let min = data.iter().copied().fold(T::infinity(), T::min);
    let max = data.iter().copied().fold(T::neg_infinity(), T::max);

    Ok(max - min)
}

/// Computes the coefficient of variation (CV).
///
/// CV is the ratio of standard deviation to mean:
/// ```text
/// CV = σ / μ
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::coefficient_of_variation;
///
/// let data = [2.0, 4.0, 6.0, 8.0];
/// let cv = coefficient_of_variation(&data, true).unwrap();
/// assert!(cv > 0.0 && cv < 1.0);
/// ```
pub fn coefficient_of_variation<T: Float>(data: &[T], sample: bool) -> Result<T, StatisticsError> {
    let m = mean(data)?;
    let sd = standard_deviation(data, sample)?;

    if m.is_zero() {
        return Err(StatisticsError::InvalidParameter(
            "Mean is zero, CV is undefined".to_string(),
        ));
    }

    Ok(sd / m)
}

// ============================================================================
// SHAPE
// ============================================================================

/// Computes the skewness of a dataset.
///
/// Skewness measures the asymmetry of the distribution:
/// ```text
/// γ₁ = E[(X - μ)³] / σ³
/// ```
///
/// - γ₁ > 0: right-skewed (tail on right)
/// - γ₁ < 0: left-skewed (tail on left)
/// - γ₁ ≈ 0: symmetric
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::skewness;
///
/// let data = [1.0, 2.0, 3.0, 4.0, 10.0];  // Right-skewed
/// let skew = skewness(&data).unwrap();
/// assert!(skew > 0.0);
/// ```
pub fn skewness<T: Float>(data: &[T]) -> Result<T, StatisticsError> {
    if data.len() < 3 {
        return Err(StatisticsError::InsufficientData {
            required: 3,
            provided: data.len(),
        });
    }

    let m = mean(data)?;
    let sd = standard_deviation(data, false)?;

    if sd.is_zero() {
        return Err(StatisticsError::InvalidParameter(
            "Standard deviation is zero".to_string(),
        ));
    }

    let n = T::from(data.len()).unwrap();
    let sum_cubed: T = data.iter()
        .map(|&x| {
            let z = (x - m) / sd;
            z * z * z
        })
        .fold(T::zero(), |acc, x| acc + x);

    Ok(sum_cubed / n)
}

/// Computes the kurtosis of a dataset.
///
/// Kurtosis measures the "tailedness" of the distribution:
/// ```text
/// γ₂ = E[(X - μ)⁴] / σ⁴
/// ```
///
/// Returns excess kurtosis (kurtosis - 3), where:
/// - Normal distribution has excess kurtosis = 0
/// - Heavy-tailed distributions have excess kurtosis > 0
/// - Light-tailed distributions have excess kurtosis < 0
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::kurtosis;
///
/// let data = [1.0, 2.0, 3.0, 4.0, 5.0];
/// let kurt = kurtosis(&data).unwrap();
/// // Uniform-like distribution has negative excess kurtosis
/// assert!(kurt < 0.0);
/// ```
pub fn kurtosis<T: Float>(data: &[T]) -> Result<T, StatisticsError> {
    if data.len() < 4 {
        return Err(StatisticsError::InsufficientData {
            required: 4,
            provided: data.len(),
        });
    }

    let m = mean(data)?;
    let sd = standard_deviation(data, false)?;

    if sd.is_zero() {
        return Err(StatisticsError::InvalidParameter(
            "Standard deviation is zero".to_string(),
        ));
    }

    let n = T::from(data.len()).unwrap();
    let sum_fourth: T = data.iter()
        .map(|&x| {
            let z = (x - m) / sd;
            let z_sq = z * z;
            z_sq * z_sq
        })
        .fold(T::zero(), |acc, x| acc + x);

    // Return excess kurtosis (kurtosis - 3)
    let kurtosis = sum_fourth / n;
    Ok(kurtosis - T::from(3.0).unwrap())
}

// ============================================================================
// ORDER STATISTICS
// ============================================================================

/// Computes a percentile of a dataset.
///
/// The pth percentile is a value below which p% of the data falls.
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::percentile;
///
/// let mut data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
/// let p50 = percentile(&mut data, 50.0).unwrap();
/// assert!((p50 - 3.0).abs() < 1e-10);
/// ```
///
/// # Parameters
///
/// - `p`: Percentile to compute (0 to 100)
pub fn percentile<T: Float>(data: &mut [T], p: T) -> Result<T, StatisticsError> {
    let zero = T::zero();
    let hundred = T::from(100.0).unwrap();

    if p < zero || p > hundred {
        return Err(StatisticsError::InvalidParameter(
            "Percentile must be between 0 and 100".to_string(),
        ));
    }

    quantile(data, p / hundred)
}

/// Computes a quantile of a dataset.
///
/// The qth quantile is a value below which a fraction q of the data falls.
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::quantile;
///
/// let mut data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
/// let q50 = quantile(&mut data, 0.5).unwrap();  // Median
/// assert!((q50 - 3.0).abs() < 1e-10);
/// ```
///
/// # Parameters
///
/// - `q`: Quantile to compute (0.0 to 1.0)
pub fn quantile<T: Float>(data: &mut [T], q: T) -> Result<T, StatisticsError> {
    if data.is_empty() {
        return Err(StatisticsError::EmptyData);
    }

    let zero = T::zero();
    let one = T::one();

    if q < zero || q > one {
        return Err(StatisticsError::InvalidParameter(
            "Quantile must be between 0 and 1".to_string(),
        ));
    }

    // Sort the data
    data.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let n = T::from(data.len()).unwrap();
    let index = q * (n - one);
    let lower_index = index.floor().to_usize().unwrap();
    let fraction = index - index.floor();

    if lower_index >= data.len() - 1 {
        Ok(data[data.len() - 1])
    } else {
        // Linear interpolation between two adjacent values
        let lower_value = data[lower_index];
        let upper_value = data[lower_index + 1];
        Ok(lower_value + fraction * (upper_value - lower_value))
    }
}

/// Computes the quartiles (Q1, Q2, Q3) of a dataset.
///
/// - Q1: 25th percentile
/// - Q2: 50th percentile (median)
/// - Q3: 75th percentile
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::quartiles;
///
/// let mut data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
/// let (q1, q2, q3) = quartiles(&mut data).unwrap();
/// assert!(q1 < q2 && q2 < q3);
/// ```
pub fn quartiles<T: Float>(data: &mut [T]) -> Result<(T, T, T), StatisticsError> {
    let q1 = quantile(data, T::from(0.25).unwrap())?;
    let q2 = quantile(data, T::from(0.50).unwrap())?;
    let q3 = quantile(data, T::from(0.75).unwrap())?;
    Ok((q1, q2, q3))
}

/// Computes the interquartile range (IQR).
///
/// IQR = Q3 - Q1, the range of the middle 50% of data.
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::interquartile_range;
///
/// let mut data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
/// let iqr = interquartile_range(&mut data).unwrap();
/// assert!(iqr > 0.0);
/// ```
pub fn interquartile_range<T: Float>(data: &mut [T]) -> Result<T, StatisticsError> {
    let (q1, _, q3) = quartiles(data)?;
    Ok(q3 - q1)
}

/// Computes the five-number summary of a dataset.
///
/// Returns (min, Q1, median, Q3, max).
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::five_number_summary;
///
/// let mut data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
/// let (min, q1, med, q3, max) = five_number_summary(&mut data).unwrap();
/// assert!(min <= q1 && q1 <= med && med <= q3 && q3 <= max);
/// ```
pub fn five_number_summary<T: Float>(data: &mut [T]) -> Result<(T, T, T, T, T), StatisticsError> {
    if data.is_empty() {
        return Err(StatisticsError::EmptyData);
    }

    data.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let min = data[0];
    let max = data[data.len() - 1];
    let (q1, q2, q3) = quartiles(data)?;

    Ok((min, q1, q2, q3, max))
}

// ============================================================================
// MULTIVARIATE
// ============================================================================

/// Computes the covariance between two datasets.
///
/// Covariance measures how two variables vary together:
/// ```text
/// cov(X,Y) = E[(X - μₓ)(Y - μᵧ)]
/// ```
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::covariance;
///
/// let x = [1.0, 2.0, 3.0, 4.0];
/// let y = [2.0, 4.0, 6.0, 8.0];
/// let cov = covariance(&x, &y, true).unwrap();
/// assert!(cov > 0.0);  // Positive covariance
/// ```
///
/// # Parameters
///
/// - `sample`: If true, divides by n-1 (sample covariance). If false, divides by n (population covariance).
pub fn covariance<T: Float>(x: &[T], y: &[T], sample: bool) -> Result<T, StatisticsError> {
    if x.is_empty() || y.is_empty() {
        return Err(StatisticsError::EmptyData);
    }
    if x.len() != y.len() {
        return Err(StatisticsError::InvalidParameter(
            "X and Y must have the same length".to_string(),
        ));
    }
    if sample && x.len() < 2 {
        return Err(StatisticsError::InsufficientData {
            required: 2,
            provided: x.len(),
        });
    }

    let mean_x = mean(x)?;
    let mean_y = mean(y)?;

    let sum: T = x.iter().zip(y.iter())
        .map(|(&xi, &yi)| (xi - mean_x) * (yi - mean_y))
        .fold(T::zero(), |acc, val| acc + val);

    let divisor = if sample {
        T::from(x.len() - 1).unwrap()
    } else {
        T::from(x.len()).unwrap()
    };

    Ok(sum / divisor)
}

/// Computes the Pearson correlation coefficient between two datasets.
///
/// Correlation measures the linear relationship between two variables:
/// ```text
/// ρ = cov(X,Y) / (σₓ σᵧ)
/// ```
///
/// Returns a value in [-1, 1]:
/// - ρ = 1: perfect positive linear relationship
/// - ρ = 0: no linear relationship
/// - ρ = -1: perfect negative linear relationship
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::correlation;
///
/// let x = [1.0, 2.0, 3.0, 4.0];
/// let y = [2.0, 4.0, 6.0, 8.0];
/// let r = correlation(&x, &y).unwrap();
/// assert!((r - 1.0).abs() < 1e-10);  // Perfect correlation
/// ```
pub fn correlation<T: Float>(x: &[T], y: &[T]) -> Result<T, StatisticsError> {
    if x.len() < 2 {
        return Err(StatisticsError::InsufficientData {
            required: 2,
            provided: x.len(),
        });
    }

    let cov = covariance(x, y, true)?;
    let sd_x = standard_deviation(x, true)?;
    let sd_y = standard_deviation(y, true)?;

    if sd_x.is_zero() || sd_y.is_zero() {
        return Err(StatisticsError::InvalidParameter(
            "Standard deviation is zero".to_string(),
        ));
    }

    Ok(cov / (sd_x * sd_y))
}

// ============================================================================
// MATRICES (Deferred from Phase 8)
// ============================================================================

/// Computes the covariance matrix for multivariate data.
///
/// For data with n observations and p variables, returns a p×p covariance matrix
/// where element (i,j) is the covariance between variables i and j.
///
/// # Mathematical Background
///
/// ```text
/// Cov(X,Y) = E[(X - μₓ)(Y - μᵧ)]
/// ```
///
/// The covariance matrix Σ is symmetric and positive semi-definite.
///
/// # Parameters
///
/// - `data`: 2D data where each row is an observation and each column is a variable.
///   Input is flattened: data[i*cols + j] = observation i, variable j
/// - `rows`: Number of observations
/// - `cols`: Number of variables
/// - `sample`: If true, uses sample covariance (divides by n-1). Otherwise population (divides by n).
///
/// # Returns
///
/// A flattened p×p covariance matrix where element at (i,j) represents cov(var_i, var_j).
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::covariance_matrix;
///
/// // 3 observations, 2 variables
/// let data = vec![
///     1.0, 2.0,  // observation 1
///     2.0, 4.0,  // observation 2
///     3.0, 6.0,  // observation 3
/// ];
/// let cov_mat = covariance_matrix(&data, 3, 2, true).unwrap();
/// // cov_mat is a 2×2 matrix (flattened to length 4)
/// ```
pub fn covariance_matrix<T: Float>(
    data: &[T],
    rows: usize,
    cols: usize,
    sample: bool,
) -> Result<Vec<T>, StatisticsError> {
    if rows == 0 || cols == 0 {
        return Err(StatisticsError::EmptyData);
    }
    if data.len() != rows * cols {
        return Err(StatisticsError::InvalidParameter(
            "Data length must equal rows × cols".to_string(),
        ));
    }
    if sample && rows < 2 {
        return Err(StatisticsError::InsufficientData {
            required: 2,
            provided: rows,
        });
    }

    // Compute means for each variable (column)
    let means: Vec<T> = (0..cols)
        .map(|col| {
            let col_sum: T = (0..rows)
                .map(|row| data[row * cols + col])
                .fold(T::zero(), |acc, x| acc + x);
            col_sum / T::from(rows).unwrap()
        })
        .collect();

    // Compute covariance matrix
    let mut cov_matrix = vec![T::zero(); cols * cols];
    let divisor = if sample {
        T::from(rows - 1).unwrap()
    } else {
        T::from(rows).unwrap()
    };

    for i in 0..cols {
        for j in 0..cols {
            let mut sum = T::zero();
            for row in 0..rows {
                let xi = data[row * cols + i] - means[i];
                let xj = data[row * cols + j] - means[j];
                sum = sum + xi * xj;
            }
            cov_matrix[i * cols + j] = sum / divisor;
        }
    }

    Ok(cov_matrix)
}

/// Computes the correlation matrix (Pearson) for multivariate data.
///
/// The correlation matrix is the normalized covariance matrix where each element
/// is divided by the product of the standard deviations.
///
/// # Mathematical Background
///
/// ```text
/// Corr(X,Y) = Cov(X,Y) / (σₓ σᵧ)
/// ```
///
/// All diagonal elements are 1, and off-diagonal elements are in [-1, 1].
///
/// # Parameters
///
/// - `data`: 2D data where each row is an observation, each column is a variable (flattened)
/// - `rows`: Number of observations
/// - `cols`: Number of variables
///
/// # Returns
///
/// A flattened p×p correlation matrix.
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::correlation_matrix;
///
/// let data = vec![
///     1.0, 2.0,
///     2.0, 4.0,
///     3.0, 6.0,
/// ];
/// let corr_mat = correlation_matrix(&data, 3, 2).unwrap();
/// // Diagonal elements should be 1.0
/// ```
pub fn correlation_matrix<T: Float>(
    data: &[T],
    rows: usize,
    cols: usize,
) -> Result<Vec<T>, StatisticsError> {
    // Get covariance matrix
    let cov_mat = covariance_matrix(data, rows, cols, true)?;

    // Compute standard deviations for each variable
    let std_devs: Vec<T> = (0..cols)
        .map(|col| {
            let col_data: Vec<T> = (0..rows)
                .map(|row| data[row * cols + col])
                .collect();
            standard_deviation(&col_data, true).unwrap()
        })
        .collect();

    // Check for zero standard deviations
    for &sd in &std_devs {
        if sd.is_zero() {
            return Err(StatisticsError::InvalidParameter(
                "Cannot compute correlation with zero standard deviation".to_string(),
            ));
        }
    }

    // Normalize covariance matrix to get correlation matrix
    let mut corr_matrix = vec![T::zero(); cols * cols];
    for i in 0..cols {
        for j in 0..cols {
            corr_matrix[i * cols + j] = cov_mat[i * cols + j] / (std_devs[i] * std_devs[j]);
        }
    }

    Ok(corr_matrix)
}

/// Computes Spearman's rank correlation coefficient.
///
/// Spearman's ρ is a non-parametric measure of rank correlation. It assesses
/// how well the relationship between two variables can be described using a
/// monotonic function.
///
/// # Mathematical Background
///
/// Spearman's ρ is the Pearson correlation applied to the ranks of the data:
/// ```text
/// ρ = 1 - (6 Σdᵢ²) / (n(n²-1))
/// ```
/// where dᵢ is the difference between the ranks of corresponding values.
///
/// For data without ties, this simplifies to the formula above. With ties,
/// we compute Pearson correlation on the ranks.
///
/// # Properties
///
/// - -1 ≤ ρ ≤ 1
/// - ρ = 1: perfect monotonic increasing relationship
/// - ρ = -1: perfect monotonic decreasing relationship
/// - ρ = 0: no monotonic relationship
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::spearman_correlation;
///
/// let x = [1.0, 2.0, 3.0, 4.0, 5.0];
/// let y = [2.0, 4.0, 6.0, 8.0, 10.0];
/// let rho = spearman_correlation(&x, &y).unwrap();
/// assert!((rho - 1.0).abs() < 1e-10); // Perfect monotonic relationship
/// ```
///
/// # References
///
/// - Spearman, C. (1904). "The proof and measurement of association between two things"
pub fn spearman_correlation<T: Float>(x: &[T], y: &[T]) -> Result<T, StatisticsError> {
    if x.is_empty() || y.is_empty() {
        return Err(StatisticsError::EmptyData);
    }
    if x.len() != y.len() {
        return Err(StatisticsError::InvalidParameter(
            "Arrays must have the same length".to_string(),
        ));
    }
    if x.len() < 2 {
        return Err(StatisticsError::InsufficientData {
            required: 2,
            provided: x.len(),
        });
    }

    // Convert to ranks
    let ranks_x = compute_ranks(x);
    let ranks_y = compute_ranks(y);

    // Compute Pearson correlation on ranks
    correlation(&ranks_x, &ranks_y)
}

/// Computes Kendall's tau correlation coefficient.
///
/// Kendall's τ is a non-parametric measure of ordinal association based on
/// the number of concordant and discordant pairs.
///
/// # Mathematical Background
///
/// ```text
/// τ = (nc - nd) / (n(n-1)/2)
/// ```
/// where:
/// - nc = number of concordant pairs
/// - nd = number of discordant pairs
/// - n = sample size
///
/// A pair (xᵢ, yᵢ) and (xⱼ, yⱼ) is:
/// - Concordant if (xᵢ - xⱼ)(yᵢ - yⱼ) > 0
/// - Discordant if (xᵢ - xⱼ)(yᵢ - yⱼ) < 0
///
/// # Properties
///
/// - -1 ≤ τ ≤ 1
/// - τ = 1: all pairs are concordant
/// - τ = -1: all pairs are discordant
/// - τ = 0: equal numbers of concordant and discordant pairs
///
/// # Complexity
///
/// - Time: O(n²) for naive implementation
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::statistics::descriptive::kendall_tau;
///
/// let x = [1.0, 2.0, 3.0, 4.0, 5.0];
/// let y = [1.0, 3.0, 2.0, 5.0, 4.0];
/// let tau = kendall_tau(&x, &y).unwrap();
/// assert!(tau > 0.0 && tau < 1.0);
/// ```
///
/// # References
///
/// - Kendall, M. G. (1938). "A new measure of rank correlation"
pub fn kendall_tau<T: Float>(x: &[T], y: &[T]) -> Result<T, StatisticsError> {
    if x.is_empty() || y.is_empty() {
        return Err(StatisticsError::EmptyData);
    }
    if x.len() != y.len() {
        return Err(StatisticsError::InvalidParameter(
            "Arrays must have the same length".to_string(),
        ));
    }
    if x.len() < 2 {
        return Err(StatisticsError::InsufficientData {
            required: 2,
            provided: x.len(),
        });
    }

    let n = x.len();
    let mut concordant = 0;
    let mut discordant = 0;

    // Count concordant and discordant pairs
    for i in 0..n {
        for j in (i + 1)..n {
            let x_diff = x[i] - x[j];
            let y_diff = y[i] - y[j];
            let product = x_diff * y_diff;

            if product > T::zero() {
                concordant += 1;
            } else if product < T::zero() {
                discordant += 1;
            }
            // If product == 0, it's a tie, not counted
        }
    }

    // τ = (nc - nd) / (n(n-1)/2)
    let total_pairs = T::from(n * (n - 1) / 2).unwrap();
    let tau = (T::from(concordant).unwrap() - T::from(discordant).unwrap()) / total_pairs;

    Ok(tau)
}

/// Helper function to compute ranks of data.
///
/// Returns a vector of ranks where the smallest value gets rank 1.
/// Handles ties by assigning average ranks.
fn compute_ranks<T: Float>(data: &[T]) -> Vec<T> {
    let n = data.len();

    // Create index-value pairs and sort by value
    let mut indexed: Vec<(usize, T)> = data.iter().copied().enumerate().collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    // Assign ranks (1-indexed)
    let mut ranks = vec![T::zero(); n];
    let mut i = 0;
    while i < n {
        // Find range of tied values
        let mut j = i + 1;
        while j < n && (indexed[j].1 - indexed[i].1).abs() < T::from(1e-10).unwrap() {
            j += 1;
        }

        // Assign average rank to all tied values
        let avg_rank = T::from(i + 1 + j).unwrap() / T::from(2.0).unwrap();
        for k in i..j {
            ranks[indexed[k].0] = avg_rank;
        }

        i = j;
    }

    ranks
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    // Central tendency tests
    #[test]
    fn test_mean() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0];
        assert_abs_diff_eq!(mean(&data).unwrap(), 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_mean_empty() {
        let data: [f64; 0] = [];
        assert!(mean(&data).is_err());
    }

    #[test]
    fn test_weighted_mean() {
        let data = [1.0, 2.0, 3.0];
        let weights = [1.0, 2.0, 1.0];
        assert_abs_diff_eq!(weighted_mean(&data, &weights).unwrap(), 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_median_odd() {
        let mut data = vec![1.0, 3.0, 2.0, 5.0, 4.0];
        assert_abs_diff_eq!(median(&mut data).unwrap(), 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_median_even() {
        let mut data = vec![1.0, 2.0, 3.0, 4.0];
        assert_abs_diff_eq!(median(&mut data).unwrap(), 2.5, epsilon = 1e-10);
    }

    #[test]
    fn test_mode() {
        let data = vec![1, 2, 2, 3, 3, 3, 4];
        assert_eq!(mode(&data).unwrap(), vec![3]);
    }

    #[test]
    fn test_geometric_mean() {
        let data = [1.0, 2.0, 4.0, 8.0];
        let gm = geometric_mean(&data).unwrap();
        let expected = (1.0 * 2.0 * 4.0 * 8.0_f64).powf(0.25);
        assert_abs_diff_eq!(gm, expected, epsilon = 1e-10);
    }

    // Dispersion tests
    #[test]
    fn test_variance() {
        let data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let var = variance(&data, false).unwrap();
        assert_abs_diff_eq!(var, 4.0, epsilon = 0.01);
    }

    #[test]
    fn test_standard_deviation() {
        let data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let sd = standard_deviation(&data, false).unwrap();
        assert_abs_diff_eq!(sd, 2.0, epsilon = 0.01);
    }

    #[test]
    fn test_range() {
        let data = [1.0, 2.0, 5.0, 3.0];
        assert_abs_diff_eq!(range(&data).unwrap(), 4.0, epsilon = 1e-10);
    }

    #[test]
    fn test_coefficient_of_variation() {
        let data = [2.0, 4.0, 6.0, 8.0];
        let cv = coefficient_of_variation(&data, false).unwrap();
        assert!(cv > 0.0 && cv < 1.0);
    }

    // Shape tests
    #[test]
    fn test_skewness() {
        let data = [1.0, 2.0, 3.0, 4.0, 10.0]; // Right-skewed
        let skew = skewness(&data).unwrap();
        assert!(skew > 0.0);
    }

    #[test]
    fn test_kurtosis() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0];
        let kurt = kurtosis(&data).unwrap();
        // Uniform-like has negative excess kurtosis
        assert!(kurt < 0.0);
    }

    // Order statistics tests
    #[test]
    fn test_percentile() {
        let mut data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert_abs_diff_eq!(percentile(&mut data, 50.0).unwrap(), 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_quantile() {
        let mut data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert_abs_diff_eq!(quantile(&mut data, 0.5).unwrap(), 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_quartiles() {
        let mut data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let (q1, q2, q3) = quartiles(&mut data).unwrap();
        assert!(q1 < q2 && q2 < q3);
    }

    #[test]
    fn test_interquartile_range() {
        let mut data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        let iqr = interquartile_range(&mut data).unwrap();
        assert!(iqr > 0.0);
    }

    #[test]
    fn test_five_number_summary() {
        let mut data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let (min, q1, med, q3, max) = five_number_summary(&mut data).unwrap();
        assert!(min <= q1 && q1 <= med && med <= q3 && q3 <= max);
        assert_eq!(min, 1.0);
        assert_eq!(max, 5.0);
    }

    // Multivariate tests
    #[test]
    fn test_covariance() {
        let x = [1.0, 2.0, 3.0, 4.0];
        let y = [2.0, 4.0, 6.0, 8.0];
        let cov = covariance(&x, &y, false).unwrap();
        assert!(cov > 0.0);
    }

    #[test]
    fn test_correlation_perfect() {
        let x = [1.0, 2.0, 3.0, 4.0];
        let y = [2.0, 4.0, 6.0, 8.0];
        let r = correlation(&x, &y).unwrap();
        assert_abs_diff_eq!(r, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_correlation_negative() {
        let x = [1.0, 2.0, 3.0, 4.0];
        let y = [8.0, 6.0, 4.0, 2.0];
        let r = correlation(&x, &y).unwrap();
        assert_abs_diff_eq!(r, -1.0, epsilon = 1e-10);
    }

    // Matrix tests (Phase 9 deferred functions)
    #[test]
    fn test_covariance_matrix() {
        // 3 observations, 2 variables
        let data = vec![
            1.0, 2.0,  // obs 1
            2.0, 4.0,  // obs 2
            3.0, 6.0,  // obs 3
        ];
        let cov_mat = covariance_matrix(&data, 3, 2, true).unwrap();

        // Should be 2×2 matrix (4 elements)
        assert_eq!(cov_mat.len(), 4);

        // Diagonal elements (variances) should be positive
        assert!(cov_mat[0] > 0.0); // var(X)
        assert!(cov_mat[3] > 0.0); // var(Y)

        // Matrix should be symmetric
        assert_abs_diff_eq!(cov_mat[1], cov_mat[2], epsilon = 1e-10);
    }

    #[test]
    fn test_correlation_matrix() {
        let data = vec![
            1.0, 2.0,
            2.0, 4.0,
            3.0, 6.0,
        ];
        let corr_mat = correlation_matrix(&data, 3, 2).unwrap();

        // Diagonal should be 1
        assert_abs_diff_eq!(corr_mat[0], 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(corr_mat[3], 1.0, epsilon = 1e-10);

        // Off-diagonal should be in [-1, 1]
        assert!(corr_mat[1] >= -1.0 && corr_mat[1] <= 1.0);
        assert!(corr_mat[2] >= -1.0 && corr_mat[2] <= 1.0);

        // Symmetric
        assert_abs_diff_eq!(corr_mat[1], corr_mat[2], epsilon = 1e-10);
    }

    #[test]
    fn test_spearman_correlation_perfect() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [2.0, 4.0, 6.0, 8.0, 10.0];
        let rho = spearman_correlation(&x, &y).unwrap();
        assert_abs_diff_eq!(rho, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_spearman_correlation_monotonic() {
        // Monotonic but not linear
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [1.0, 4.0, 9.0, 16.0, 25.0]; // y = x²
        let rho = spearman_correlation(&x, &y).unwrap();
        assert_abs_diff_eq!(rho, 1.0, epsilon = 1e-10); // Perfect monotonic
    }

    #[test]
    fn test_kendall_tau_perfect() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [1.0, 2.0, 3.0, 4.0, 5.0];
        let tau = kendall_tau(&x, &y).unwrap();
        assert_abs_diff_eq!(tau, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_kendall_tau_partial() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [1.0, 3.0, 2.0, 5.0, 4.0]; // Some inversions
        let tau = kendall_tau(&x, &y).unwrap();
        assert!(tau > 0.0 && tau < 1.0); // Positive but not perfect
    }

    #[test]
    fn test_kendall_tau_negative() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [5.0, 4.0, 3.0, 2.0, 1.0]; // Perfect negative
        let tau = kendall_tau(&x, &y).unwrap();
        assert_abs_diff_eq!(tau, -1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_compute_ranks() {
        let data = [3.0, 1.0, 4.0, 1.0, 5.0];
        let ranks = compute_ranks(&data);

        // Check that smallest values get smallest ranks
        // 1.0 appears twice at positions 1,3 - should get average rank 1.5
        // 3.0 at position 0 - should get rank 3
        // 4.0 at position 2 - should get rank 4
        // 5.0 at position 4 - should get rank 5
        assert_abs_diff_eq!(ranks[1], 1.5, epsilon = 1e-10);
        assert_abs_diff_eq!(ranks[3], 1.5, epsilon = 1e-10);
        assert_abs_diff_eq!(ranks[0], 3.0, epsilon = 1e-10);
        assert_abs_diff_eq!(ranks[2], 4.0, epsilon = 1e-10);
        assert_abs_diff_eq!(ranks[4], 5.0, epsilon = 1e-10);
    }
}
