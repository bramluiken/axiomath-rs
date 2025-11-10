//! Information theory: entropy, mutual information, and divergence measures.
//!
//! This module implements fundamental concepts from information theory, measuring
//! the information content, uncertainty, and similarity of probability distributions.
//!
//! # Theory
//!
//! Information theory, founded by Claude Shannon in 1948, provides a mathematical
//! framework for quantifying information. Key concepts include:
//!
//! - **Entropy**: Measures the average uncertainty in a random variable
//! - **Mutual Information**: Measures the dependence between two variables
//! - **Divergence**: Measures the difference between probability distributions
//!
//! # Module Organization
//!
//! - **Entropy Measures**: Shannon entropy, joint entropy, conditional entropy
//! - **Information Measures**: Mutual information, normalized mutual information
//! - **Divergence Measures**: KL divergence, JS divergence, cross-entropy
//!
//! # Examples
//!
//! ## Shannon Entropy
//!
//! ```
//! use axiomath::statistics::information::shannon_entropy;
//!
//! // Fair coin: maximum entropy for binary variable
//! let probs = [0.5, 0.5];
//! let h = shannon_entropy(&probs).unwrap();
//! assert!((h - 1.0).abs() < 1e-10); // H = 1 bit
//! ```
//!
//! ## Kullback-Leibler Divergence
//!
//! ```
//! use axiomath::statistics::information::kl_divergence;
//!
//! let p = [0.5, 0.5];
//! let q = [0.6, 0.4];
//! let div = kl_divergence(&p, &q).unwrap();
//! assert!(div > 0.0); // KL divergence is always non-negative
//! ```
//!
//! # References
//!
//! - Shannon, C. E. (1948). "A Mathematical Theory of Communication"
//! - Cover, T. M., & Thomas, J. A. (2006). "Elements of Information Theory"
//! - MacKay, D. J. (2003). "Information Theory, Inference, and Learning Algorithms"

use crate::core::Float;
use crate::statistics::StatisticsError;
use num_traits::Float as NumFloat;

// ============================================================================
// ENTROPY MEASURES
// ============================================================================

/// Computes the Shannon entropy of a discrete probability distribution.
///
/// Shannon entropy measures the average information content (uncertainty) in a
/// random variable. It's maximal when the distribution is uniform.
///
/// # Mathematical Background
///
/// For a discrete random variable X with probability mass function p(x):
/// ```text
/// H(X) = -Σ p(x) log₂ p(x)
/// ```
///
/// By convention, 0 log 0 = 0.
///
/// # Properties
///
/// - H(X) ≥ 0 (non-negative)
/// - H(X) = 0 iff X is deterministic (one probability is 1, rest are 0)
/// - H(X) ≤ log₂(n) where n is the number of outcomes (maximum for uniform distribution)
///
/// # Units
///
/// Uses base-2 logarithm, so entropy is measured in **bits**.
///
/// # Complexity
///
/// - Time: O(n) where n is the number of probabilities
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::statistics::information::shannon_entropy;
///
/// // Fair coin: H = 1 bit
/// let fair_coin = [0.5, 0.5];
/// let h1 = shannon_entropy(&fair_coin).unwrap();
/// assert!((h1 - 1.0).abs() < 1e-10);
///
/// // Biased coin: H < 1 bit
/// let biased_coin = [0.9, 0.1];
/// let h2 = shannon_entropy(&biased_coin).unwrap();
/// assert!(h2 < h1);
///
/// // Deterministic: H = 0
/// let deterministic = [1.0, 0.0];
/// let h3 = shannon_entropy(&deterministic).unwrap();
/// assert!(h3.abs() < 1e-10);
/// ```
///
/// # Errors
///
/// Returns error if:
/// - Probabilities array is empty
/// - Any probability is negative
/// - Probabilities don't sum to 1 (within tolerance)
///
/// # References
///
/// - Shannon, C. E. (1948). "A Mathematical Theory of Communication". Bell System Technical Journal.
pub fn shannon_entropy<T: Float>(probabilities: &[T]) -> Result<T, StatisticsError> {
    if probabilities.is_empty() {
        return Err(StatisticsError::EmptyData);
    }

    // Check that all probabilities are non-negative
    for &p in probabilities {
        if p < T::zero() {
            return Err(StatisticsError::InvalidProbability);
        }
    }

    // Check that probabilities sum to 1 (within tolerance)
    let sum: T = probabilities.iter().copied().fold(T::zero(), |acc, p| acc + p);
    let tolerance = T::from(1e-6).unwrap();
    if (sum - T::one()).abs() > tolerance {
        return Err(StatisticsError::InvalidProbability);
    }

    // Compute entropy: H = -Σ p log₂ p
    // Use convention that 0 log 0 = 0
    let log2 = T::from(2.0_f64.ln()).unwrap();
    let entropy: T = probabilities.iter()
        .filter(|&&p| p > T::zero())
        .map(|&p| {
            let log_p = p.ln() / log2; // Convert to log base 2
            -p * log_p
        })
        .fold(T::zero(), |acc, x| acc + x);

    Ok(entropy)
}

/// Computes the joint entropy of a two-dimensional probability distribution.
///
/// Joint entropy measures the total uncertainty in two random variables together.
///
/// # Mathematical Background
///
/// For random variables X and Y with joint probability distribution p(x,y):
/// ```text
/// H(X,Y) = -ΣΣ p(x,y) log₂ p(x,y)
/// ```
///
/// # Properties
///
/// - H(X,Y) ≤ H(X) + H(Y) (equality iff X and Y are independent)
/// - H(X,Y) ≥ max(H(X), H(Y))
///
/// # Parameters
///
/// - `joint_probabilities`: Flattened 2D array of probabilities p(x,y)
///
/// # Examples
///
/// ```
/// use axiomath::statistics::information::joint_entropy;
///
/// // Independent events: p(x,y) = p(x)p(y)
/// // X: [0.5, 0.5], Y: [0.5, 0.5]
/// let joint_probs = [0.25, 0.25, 0.25, 0.25];
/// let h_xy = joint_entropy(&joint_probs).unwrap();
/// assert!((h_xy - 2.0).abs() < 1e-10); // H(X,Y) = H(X) + H(Y) = 2 bits
/// ```
pub fn joint_entropy<T: Float>(joint_probabilities: &[T]) -> Result<T, StatisticsError> {
    // Joint entropy is just Shannon entropy of the joint distribution
    shannon_entropy(joint_probabilities)
}

/// Computes the conditional entropy H(Y|X).
///
/// Conditional entropy measures the average uncertainty in Y given knowledge of X.
///
/// # Mathematical Background
///
/// ```text
/// H(Y|X) = H(X,Y) - H(X)
///        = -ΣΣ p(x,y) log₂ p(y|x)
/// ```
///
/// # Properties
///
/// - H(Y|X) ≥ 0
/// - H(Y|X) ≤ H(Y) (knowledge can only reduce uncertainty)
/// - H(Y|X) = 0 iff Y is a deterministic function of X
/// - H(Y|X) = H(Y) iff X and Y are independent
///
/// # Parameters
///
/// - `joint_probs`: Joint probability distribution p(x,y) as flattened 2D array
/// - `rows`: Number of rows (values of X)
/// - `cols`: Number of columns (values of Y)
///
/// # Examples
///
/// ```
/// use axiomath::statistics::information::conditional_entropy;
///
/// // Perfect dependence: Y = X
/// // For 2×2 case: p(0,0)=0.5, p(1,1)=0.5, rest are 0
/// let joint = [0.5, 0.0, 0.0, 0.5];
/// let h_y_given_x = conditional_entropy(&joint, 2, 2).unwrap();
/// assert!(h_y_given_x.abs() < 1e-10); // H(Y|X) = 0
/// ```
pub fn conditional_entropy<T: Float>(
    joint_probs: &[T],
    rows: usize,
    cols: usize,
) -> Result<T, StatisticsError> {
    if joint_probs.len() != rows * cols {
        return Err(StatisticsError::InvalidParameter(
            "Joint probabilities size must equal rows × cols".to_string(),
        ));
    }

    // Compute marginal distribution p(x) by summing over y
    let marginal_x: Vec<T> = (0..rows)
        .map(|i| {
            (0..cols)
                .map(|j| joint_probs[i * cols + j])
                .fold(T::zero(), |acc, p| acc + p)
        })
        .collect();

    // Compute H(X,Y)
    let h_xy = joint_entropy(joint_probs)?;

    // Compute H(X)
    let h_x = shannon_entropy(&marginal_x)?;

    // H(Y|X) = H(X,Y) - H(X)
    Ok(h_xy - h_x)
}

// ============================================================================
// INFORMATION MEASURES
// ============================================================================

/// Computes the mutual information between two random variables.
///
/// Mutual information measures how much information one variable provides about another.
/// It quantifies the reduction in uncertainty about one variable given knowledge of the other.
///
/// # Mathematical Background
///
/// ```text
/// I(X;Y) = H(X) + H(Y) - H(X,Y)
///        = H(X) - H(X|Y)
///        = H(Y) - H(Y|X)
///        = ΣΣ p(x,y) log₂(p(x,y) / (p(x)p(y)))
/// ```
///
/// # Properties
///
/// - I(X;Y) ≥ 0 (non-negative)
/// - I(X;Y) = 0 iff X and Y are independent
/// - I(X;Y) = I(Y;X) (symmetric)
/// - I(X;Y) ≤ min(H(X), H(Y))
///
/// # Parameters
///
/// - `joint_probs`: Joint probability distribution p(x,y) as flattened 2D array
/// - `rows`: Number of rows (values of X)
/// - `cols`: Number of columns (values of Y)
///
/// # Complexity
///
/// - Time: O(n·m) where n = rows, m = cols
/// - Space: O(n + m) for marginal distributions
///
/// # Examples
///
/// ```
/// use axiomath::statistics::information::mutual_information;
///
/// // Independent variables: I(X;Y) = 0
/// let independent = [0.25, 0.25, 0.25, 0.25];
/// let mi1 = mutual_information(&independent, 2, 2).unwrap();
/// assert!(mi1.abs() < 1e-10);
///
/// // Perfect dependence: I(X;Y) = H(X) = H(Y)
/// let dependent = [0.5, 0.0, 0.0, 0.5];
/// let mi2 = mutual_information(&dependent, 2, 2).unwrap();
/// assert!((mi2 - 1.0).abs() < 1e-10); // I(X;Y) = 1 bit
/// ```
///
/// # References
///
/// - Cover & Thomas (2006). "Elements of Information Theory", Chapter 2.
pub fn mutual_information<T: Float>(
    joint_probs: &[T],
    rows: usize,
    cols: usize,
) -> Result<T, StatisticsError> {
    if joint_probs.len() != rows * cols {
        return Err(StatisticsError::InvalidParameter(
            "Joint probabilities size must equal rows × cols".to_string(),
        ));
    }

    // Compute marginal distributions
    let marginal_x: Vec<T> = (0..rows)
        .map(|i| {
            (0..cols)
                .map(|j| joint_probs[i * cols + j])
                .fold(T::zero(), |acc, p| acc + p)
        })
        .collect();

    let marginal_y: Vec<T> = (0..cols)
        .map(|j| {
            (0..rows)
                .map(|i| joint_probs[i * cols + j])
                .fold(T::zero(), |acc, p| acc + p)
        })
        .collect();

    // Compute H(X) and H(Y)
    let h_x = shannon_entropy(&marginal_x)?;
    let h_y = shannon_entropy(&marginal_y)?;

    // Compute H(X,Y)
    let h_xy = joint_entropy(joint_probs)?;

    // I(X;Y) = H(X) + H(Y) - H(X,Y)
    Ok(h_x + h_y - h_xy)
}

/// Computes the normalized mutual information.
///
/// Normalized mutual information scales MI to the range [0, 1] by dividing
/// by the average of the individual entropies.
///
/// # Mathematical Background
///
/// ```text
/// NMI(X;Y) = I(X;Y) / ((H(X) + H(Y)) / 2)
/// ```
///
/// Alternative normalizations exist (e.g., divide by max or min of H(X), H(Y)).
///
/// # Properties
///
/// - 0 ≤ NMI ≤ 1
/// - NMI = 0 iff X and Y are independent
/// - NMI = 1 iff X and Y are perfectly dependent
///
/// # Examples
///
/// ```
/// use axiomath::statistics::information::normalized_mutual_information;
///
/// // Perfect dependence: NMI = 1
/// let dependent = [0.5, 0.0, 0.0, 0.5];
/// let nmi = normalized_mutual_information(&dependent, 2, 2).unwrap();
/// assert!((nmi - 1.0).abs() < 1e-10);
/// ```
///
/// # References
///
/// - Strehl & Ghosh (2002). "Cluster ensembles: A knowledge reuse framework"
pub fn normalized_mutual_information<T: Float>(
    joint_probs: &[T],
    rows: usize,
    cols: usize,
) -> Result<T, StatisticsError> {
    // Compute marginal distributions
    let marginal_x: Vec<T> = (0..rows)
        .map(|i| {
            (0..cols)
                .map(|j| joint_probs[i * cols + j])
                .fold(T::zero(), |acc, p| acc + p)
        })
        .collect();

    let marginal_y: Vec<T> = (0..cols)
        .map(|j| {
            (0..rows)
                .map(|i| joint_probs[i * cols + j])
                .fold(T::zero(), |acc, p| acc + p)
        })
        .collect();

    // Compute entropies
    let h_x = shannon_entropy(&marginal_x)?;
    let h_y = shannon_entropy(&marginal_y)?;

    // Compute mutual information
    let mi = mutual_information(joint_probs, rows, cols)?;

    // Normalize by average entropy
    let avg_entropy = (h_x + h_y) / T::from(2.0).unwrap();

    if avg_entropy.is_zero() {
        // Both distributions are deterministic
        return Ok(T::zero());
    }

    Ok(mi / avg_entropy)
}

// ============================================================================
// DIVERGENCE MEASURES
// ============================================================================

/// Computes the Kullback-Leibler (KL) divergence from Q to P.
///
/// KL divergence measures how one probability distribution differs from another.
/// It quantifies the extra bits needed to encode data from P using a code optimized for Q.
///
/// # Mathematical Background
///
/// ```text
/// D_KL(P || Q) = Σ p(x) log₂(p(x) / q(x))
/// ```
///
/// # Properties
///
/// - D_KL(P || Q) ≥ 0 (non-negative, Gibbs' inequality)
/// - D_KL(P || Q) = 0 iff P = Q
/// - **NOT symmetric**: D_KL(P || Q) ≠ D_KL(Q || P) in general
/// - **NOT a metric**: doesn't satisfy triangle inequality
///
/// # Interpretation
///
/// - D_KL(P || Q): "Surprise" when using distribution Q to model data from P
/// - Measures the information lost when Q is used to approximate P
///
/// # Complexity
///
/// - Time: O(n) where n is the number of probability values
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::statistics::information::kl_divergence;
///
/// // Identical distributions: KL = 0
/// let p = [0.5, 0.5];
/// let q = [0.5, 0.5];
/// let div1 = kl_divergence(&p, &q).unwrap();
/// assert!(div1.abs() < 1e-10);
///
/// // Different distributions: KL > 0
/// let p = [0.5, 0.5];
/// let q = [0.6, 0.4];
/// let div2 = kl_divergence(&p, &q).unwrap();
/// assert!(div2 > 0.0);
///
/// // Asymmetry: D(P||Q) ≠ D(Q||P)
/// let div3 = kl_divergence(&q, &p).unwrap();
/// assert!((div2 - div3).abs() > 1e-10);
/// ```
///
/// # Errors
///
/// Returns error if:
/// - Arrays have different lengths
/// - Any probability is negative
/// - Q has a zero where P is non-zero (would cause division by zero)
/// - Probabilities don't sum to 1
///
/// # References
///
/// - Kullback, S., & Leibler, R. A. (1951). "On Information and Sufficiency"
/// - Cover & Thomas (2006). "Elements of Information Theory", Section 2.3
pub fn kl_divergence<T: Float>(p: &[T], q: &[T]) -> Result<T, StatisticsError> {
    if p.is_empty() || q.is_empty() {
        return Err(StatisticsError::EmptyData);
    }
    if p.len() != q.len() {
        return Err(StatisticsError::InvalidParameter(
            "Distributions must have the same length".to_string(),
        ));
    }

    // Check validity of probabilities
    for (&pi, &qi) in p.iter().zip(q.iter()) {
        if pi < T::zero() || qi < T::zero() {
            return Err(StatisticsError::InvalidProbability);
        }
        // Check for division by zero: if p(x) > 0 but q(x) = 0, KL divergence is infinite
        if pi > T::zero() && qi.is_zero() {
            return Err(StatisticsError::InvalidParameter(
                "Q cannot have zero probability where P is non-zero".to_string(),
            ));
        }
    }

    // Check that probabilities sum to 1
    let sum_p: T = p.iter().copied().fold(T::zero(), |acc, x| acc + x);
    let sum_q: T = q.iter().copied().fold(T::zero(), |acc, x| acc + x);
    let tolerance = T::from(1e-6).unwrap();
    if (sum_p - T::one()).abs() > tolerance || (sum_q - T::one()).abs() > tolerance {
        return Err(StatisticsError::InvalidProbability);
    }

    // Compute KL divergence: Σ p(x) log₂(p(x) / q(x))
    let log2 = T::from(2.0_f64.ln()).unwrap();
    let divergence: T = p.iter()
        .zip(q.iter())
        .filter(|(&pi, &qi)| pi > T::zero() && qi > T::zero())
        .map(|(&pi, &qi)| {
            let ratio = pi / qi;
            let log_ratio = ratio.ln() / log2;
            pi * log_ratio
        })
        .fold(T::zero(), |acc, x| acc + x);

    Ok(divergence)
}

/// Computes the Jensen-Shannon (JS) divergence between two distributions.
///
/// JS divergence is a symmetrized and smoothed version of KL divergence.
///
/// # Mathematical Background
///
/// ```text
/// JS(P || Q) = (1/2) D_KL(P || M) + (1/2) D_KL(Q || M)
///
/// where M = (P + Q) / 2
/// ```
///
/// # Properties
///
/// - 0 ≤ JS(P || Q) ≤ 1 (bounded)
/// - JS(P || Q) = 0 iff P = Q
/// - **Symmetric**: JS(P || Q) = JS(Q || P)
/// - √JS is a valid metric (satisfies triangle inequality)
///
/// # Complexity
///
/// - Time: O(n) where n is the number of probability values
/// - Space: O(n) for storing the midpoint distribution
///
/// # Examples
///
/// ```
/// use axiomath::statistics::information::js_divergence;
///
/// // Identical distributions: JS = 0
/// let p = [0.5, 0.5];
/// let q = [0.5, 0.5];
/// let div1 = js_divergence(&p, &q).unwrap();
/// assert!(div1.abs() < 1e-10);
///
/// // Different distributions: 0 < JS ≤ 1
/// let p = [0.5, 0.5];
/// let q = [0.9, 0.1];
/// let div2 = js_divergence(&p, &q).unwrap();
/// assert!(div2 > 0.0 && div2 <= 1.0);
///
/// // Symmetry: JS(P||Q) = JS(Q||P)
/// let div3 = js_divergence(&q, &p).unwrap();
/// assert!((div2 - div3).abs() < 1e-10);
/// ```
///
/// # References
///
/// - Lin, J. (1991). "Divergence measures based on the Shannon entropy"
/// - Endres & Schindelin (2003). "A new metric for probability distributions"
pub fn js_divergence<T: Float>(p: &[T], q: &[T]) -> Result<T, StatisticsError> {
    if p.is_empty() || q.is_empty() {
        return Err(StatisticsError::EmptyData);
    }
    if p.len() != q.len() {
        return Err(StatisticsError::InvalidParameter(
            "Distributions must have the same length".to_string(),
        ));
    }

    // Compute midpoint distribution M = (P + Q) / 2
    let two = T::from(2.0).unwrap();
    let m: Vec<T> = p.iter()
        .zip(q.iter())
        .map(|(&pi, &qi)| (pi + qi) / two)
        .collect();

    // Compute JS(P || Q) = (1/2) D_KL(P || M) + (1/2) D_KL(Q || M)
    let kl_pm = kl_divergence(p, &m)?;
    let kl_qm = kl_divergence(q, &m)?;

    Ok((kl_pm + kl_qm) / two)
}

/// Computes the cross-entropy from P to Q.
///
/// Cross-entropy measures the average number of bits needed to encode data from P
/// using a code optimized for Q.
///
/// # Mathematical Background
///
/// ```text
/// H(P, Q) = -Σ p(x) log₂ q(x)
///         = H(P) + D_KL(P || Q)
/// ```
///
/// # Relationship to Other Measures
///
/// - Cross-entropy ≥ Entropy (equality iff P = Q)
/// - Cross-entropy = Entropy + KL divergence
///
/// # Applications
///
/// - Machine learning loss function (classification)
/// - Measures coding inefficiency when using wrong distribution
///
/// # Complexity
///
/// - Time: O(n) where n is the number of probability values
/// - Space: O(1)
///
/// # Examples
///
/// ```
/// use axiomath::statistics::information::cross_entropy;
///
/// // Same distribution: H(P,Q) = H(P)
/// let p = [0.5, 0.5];
/// let q = [0.5, 0.5];
/// let ce1 = cross_entropy(&p, &q).unwrap();
/// assert!((ce1 - 1.0).abs() < 1e-10); // Equal to H(P)
///
/// // Different distributions: H(P,Q) > H(P)
/// let p = [0.5, 0.5];
/// let q = [0.9, 0.1];
/// let ce2 = cross_entropy(&p, &q).unwrap();
/// assert!(ce2 > 1.0);
/// ```
///
/// # References
///
/// - Cover & Thomas (2006). "Elements of Information Theory", Section 2.8
/// - Goodfellow et al. (2016). "Deep Learning", Section 3.13
pub fn cross_entropy<T: Float>(p: &[T], q: &[T]) -> Result<T, StatisticsError> {
    if p.is_empty() || q.is_empty() {
        return Err(StatisticsError::EmptyData);
    }
    if p.len() != q.len() {
        return Err(StatisticsError::InvalidParameter(
            "Distributions must have the same length".to_string(),
        ));
    }

    // Check validity of probabilities
    for (&pi, &qi) in p.iter().zip(q.iter()) {
        if pi < T::zero() || qi < T::zero() {
            return Err(StatisticsError::InvalidProbability);
        }
        // Check for log(0): if p(x) > 0 but q(x) = 0, cross-entropy is infinite
        if pi > T::zero() && qi.is_zero() {
            return Err(StatisticsError::InvalidParameter(
                "Q cannot have zero probability where P is non-zero".to_string(),
            ));
        }
    }

    // Compute cross-entropy: H(P,Q) = -Σ p(x) log₂ q(x)
    let log2 = T::from(2.0_f64.ln()).unwrap();
    let ce: T = p.iter()
        .zip(q.iter())
        .filter(|(&pi, &qi)| pi > T::zero() && qi > T::zero())
        .map(|(&pi, &qi)| {
            let log_qi = qi.ln() / log2;
            -pi * log_qi
        })
        .fold(T::zero(), |acc, x| acc + x);

    Ok(ce)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    // Entropy tests
    #[test]
    fn test_shannon_entropy_fair_coin() {
        let probs = [0.5, 0.5];
        let h = shannon_entropy(&probs).unwrap();
        assert_abs_diff_eq!(h, 1.0, epsilon = 1e-10); // 1 bit
    }

    #[test]
    fn test_shannon_entropy_deterministic() {
        let probs = [1.0, 0.0, 0.0];
        let h = shannon_entropy(&probs).unwrap();
        assert_abs_diff_eq!(h, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_shannon_entropy_uniform() {
        let probs = [0.25, 0.25, 0.25, 0.25];
        let h = shannon_entropy(&probs).unwrap();
        assert_abs_diff_eq!(h, 2.0, epsilon = 1e-10); // log₂(4) = 2 bits
    }

    #[test]
    fn test_shannon_entropy_invalid_probabilities() {
        let probs = [0.5, 0.6]; // Don't sum to 1
        assert!(shannon_entropy(&probs).is_err());

        let probs = [0.5, -0.5]; // Negative
        assert!(shannon_entropy(&probs).is_err());
    }

    #[test]
    fn test_joint_entropy_independent() {
        // Independent: H(X,Y) = H(X) + H(Y)
        let joint = [0.25, 0.25, 0.25, 0.25];
        let h_xy = joint_entropy(&joint).unwrap();
        assert_abs_diff_eq!(h_xy, 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_conditional_entropy_perfect_dependence() {
        // Perfect dependence: Y = X, so H(Y|X) = 0
        let joint = [0.5, 0.0, 0.0, 0.5];
        let h_y_given_x = conditional_entropy(&joint, 2, 2).unwrap();
        assert_abs_diff_eq!(h_y_given_x, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_conditional_entropy_independent() {
        // Independent: H(Y|X) = H(Y)
        let joint = [0.25, 0.25, 0.25, 0.25];
        let h_y_given_x = conditional_entropy(&joint, 2, 2).unwrap();
        assert_abs_diff_eq!(h_y_given_x, 1.0, epsilon = 1e-10); // H(Y) = 1
    }

    // Mutual information tests
    #[test]
    fn test_mutual_information_independent() {
        // Independent variables: I(X;Y) = 0
        let joint = [0.25, 0.25, 0.25, 0.25];
        let mi = mutual_information(&joint, 2, 2).unwrap();
        assert_abs_diff_eq!(mi, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_mutual_information_perfect_dependence() {
        // Perfect dependence: I(X;Y) = H(X) = H(Y) = 1
        let joint = [0.5, 0.0, 0.0, 0.5];
        let mi = mutual_information(&joint, 2, 2).unwrap();
        assert_abs_diff_eq!(mi, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_normalized_mutual_information() {
        // Perfect dependence: NMI = 1
        let joint = [0.5, 0.0, 0.0, 0.5];
        let nmi = normalized_mutual_information(&joint, 2, 2).unwrap();
        assert_abs_diff_eq!(nmi, 1.0, epsilon = 1e-10);

        // Independent: NMI = 0
        let joint = [0.25, 0.25, 0.25, 0.25];
        let nmi = normalized_mutual_information(&joint, 2, 2).unwrap();
        assert_abs_diff_eq!(nmi, 0.0, epsilon = 1e-10);
    }

    // KL divergence tests
    #[test]
    fn test_kl_divergence_identical() {
        let p = [0.5, 0.5];
        let q = [0.5, 0.5];
        let div = kl_divergence(&p, &q).unwrap();
        assert_abs_diff_eq!(div, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_kl_divergence_non_negative() {
        let p = [0.5, 0.5];
        let q = [0.6, 0.4];
        let div = kl_divergence(&p, &q).unwrap();
        assert!(div >= 0.0);
    }

    #[test]
    fn test_kl_divergence_asymmetry() {
        let p = [0.5, 0.5];
        let q = [0.9, 0.1];
        let div_pq = kl_divergence(&p, &q).unwrap();
        let div_qp = kl_divergence(&q, &p).unwrap();
        assert!((div_pq - div_qp).abs() > 1e-6); // Not equal
    }

    #[test]
    fn test_kl_divergence_zero_in_q() {
        let p = [0.5, 0.5];
        let q = [1.0, 0.0]; // q has 0 where p is non-zero
        assert!(kl_divergence(&p, &q).is_err());
    }

    // JS divergence tests
    #[test]
    fn test_js_divergence_identical() {
        let p = [0.5, 0.5];
        let q = [0.5, 0.5];
        let div = js_divergence(&p, &q).unwrap();
        assert_abs_diff_eq!(div, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_js_divergence_symmetry() {
        let p = [0.5, 0.5];
        let q = [0.9, 0.1];
        let div_pq = js_divergence(&p, &q).unwrap();
        let div_qp = js_divergence(&q, &p).unwrap();
        assert_abs_diff_eq!(div_pq, div_qp, epsilon = 1e-10);
    }

    #[test]
    fn test_js_divergence_bounded() {
        let p = [0.5, 0.5];
        let q = [0.9, 0.1];
        let div = js_divergence(&p, &q).unwrap();
        assert!(div >= 0.0 && div <= 1.0);
    }

    // Cross-entropy tests
    #[test]
    fn test_cross_entropy_same_distribution() {
        let p = [0.5, 0.5];
        let q = [0.5, 0.5];
        let ce = cross_entropy(&p, &q).unwrap();
        let h = shannon_entropy(&p).unwrap();
        assert_abs_diff_eq!(ce, h, epsilon = 1e-10);
    }

    #[test]
    fn test_cross_entropy_greater_than_entropy() {
        let p = [0.5, 0.5];
        let q = [0.9, 0.1];
        let ce = cross_entropy(&p, &q).unwrap();
        let h = shannon_entropy(&p).unwrap();
        assert!(ce >= h);
    }

    #[test]
    fn test_cross_entropy_relationship_to_kl() {
        // H(P,Q) = H(P) + D_KL(P||Q)
        let p = [0.5, 0.5];
        let q = [0.6, 0.4];
        let ce = cross_entropy(&p, &q).unwrap();
        let h = shannon_entropy(&p).unwrap();
        let kl = kl_divergence(&p, &q).unwrap();
        assert_abs_diff_eq!(ce, h + kl, epsilon = 1e-10);
    }

    #[test]
    fn test_cross_entropy_zero_in_q() {
        let p = [0.5, 0.5];
        let q = [1.0, 0.0]; // q has 0 where p is non-zero
        assert!(cross_entropy(&p, &q).is_err());
    }
}
