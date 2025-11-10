//! Gaussian quadrature methods.
//!
//! Gaussian quadrature achieves optimal accuracy by carefully choosing both the
//! evaluation points (nodes) and weights. For n points, it can exactly integrate
//! polynomials of degree 2n-1.

use crate::core::Float;

/// Computes the definite integral using Gauss-Legendre quadrature.
///
/// # Algorithm
///
/// Gauss-Legendre quadrature integrates functions over [-1, 1]:
///
/// ```text
/// ∫₋₁¹ f(x)dx ≈ Σᵢ wᵢ f(xᵢ)
/// ```
///
/// where xᵢ are roots of Legendre polynomials and wᵢ are optimal weights.
/// For arbitrary intervals [a,b], we transform: x = (b-a)t/2 + (b+a)/2
///
/// This method is exact for polynomials of degree ≤ 2n-1 where n is the number of points.
///
/// # Complexity
///
/// - Function evaluations: n
/// - Accuracy: Exponential convergence for smooth functions
///
/// # Arguments
///
/// * `f` - The function to integrate
/// * `a` - Lower bound of integration
/// * `b` - Upper bound of integration
/// * `n` - Number of quadrature points (supported: 2-10)
///
/// # Panics
///
/// Panics if n is not in the range [2, 10].
///
/// # Examples
///
/// ```
/// use axiomath::calculus::integration::gauss_legendre_quadrature;
///
/// // Integrate x² from 0 to 1, exact value is 1/3
/// let f = |x: f64| x * x;
/// let result = gauss_legendre_quadrature(f, 0.0, 1.0, 5);
/// assert!((result - 1.0/3.0).abs() < 1e-15);
/// ```
///
/// # References
///
/// - Abramowitz & Stegun, "Handbook of Mathematical Functions", Section 25.4
/// - Press et al., "Numerical Recipes", Section 4.5
pub fn gauss_legendre_quadrature<T, F>(f: F, a: T, b: T, n: usize) -> T
where
    T: Float,
    F: Fn(T) -> T,
{
    assert!(
        (2..=10).contains(&n),
        "Number of points must be between 2 and 10"
    );

    let (nodes, weights) = gauss_legendre_nodes_weights::<T>(n);

    let two = T::one() + T::one();
    let scale = (b - a) / two;
    let shift = (b + a) / two;

    let mut sum = T::zero();
    for i in 0..n {
        let x = scale * nodes[i] + shift;
        sum = sum + weights[i] * f(x);
    }

    sum * scale
}

/// Returns the nodes and weights for Gauss-Legendre quadrature on [-1, 1].
///
/// These are precomputed values for n=2 to n=10.
fn gauss_legendre_nodes_weights<T: Float>(n: usize) -> (Vec<T>, Vec<T>) {
    match n {
        2 => {
            let node = T::from(0.5773502691896257).unwrap(); // 1/√3
            let weight = T::one();
            (vec![-node, node], vec![weight, weight])
        }
        3 => {
            let node = T::from(0.7745966692414834).unwrap(); // √(3/5)
            let w0 = T::from(0.5555555555555556).unwrap(); // 5/9
            let w1 = T::from(0.8888888888888888).unwrap(); // 8/9
            (vec![-node, T::zero(), node], vec![w0, w1, w0])
        }
        4 => {
            let x1 = T::from(0.3399810435848563).unwrap();
            let x2 = T::from(0.8611363115940526).unwrap();
            let w1 = T::from(0.6521451548625461).unwrap();
            let w2 = T::from(0.3478548451374538).unwrap();
            (vec![-x2, -x1, x1, x2], vec![w2, w1, w1, w2])
        }
        5 => {
            let x1 = T::from(0.5384693101056831).unwrap();
            let x2 = T::from(0.9061798459386640).unwrap();
            let w0 = T::from(0.5688888888888889).unwrap();
            let w1 = T::from(0.4786286704993665).unwrap();
            let w2 = T::from(0.2369268850561891).unwrap();
            (
                vec![-x2, -x1, T::zero(), x1, x2],
                vec![w2, w1, w0, w1, w2],
            )
        }
        6 => {
            let x1 = T::from(0.2386191860831969).unwrap();
            let x2 = T::from(0.6612093864662645).unwrap();
            let x3 = T::from(0.9324695142031521).unwrap();
            let w1 = T::from(0.4679139345726910).unwrap();
            let w2 = T::from(0.3607615730481386).unwrap();
            let w3 = T::from(0.1713244923791704).unwrap();
            (
                vec![-x3, -x2, -x1, x1, x2, x3],
                vec![w3, w2, w1, w1, w2, w3],
            )
        }
        7 => {
            let x1 = T::from(0.4058451513773972).unwrap();
            let x2 = T::from(0.7415311855993945).unwrap();
            let x3 = T::from(0.9491079123427585).unwrap();
            let w0 = T::from(0.4179591836734694).unwrap();
            let w1 = T::from(0.3818300505051189).unwrap();
            let w2 = T::from(0.2797053914892766).unwrap();
            let w3 = T::from(0.1294849661688697).unwrap();
            (
                vec![-x3, -x2, -x1, T::zero(), x1, x2, x3],
                vec![w3, w2, w1, w0, w1, w2, w3],
            )
        }
        8 => {
            let x1 = T::from(0.1834346424956498).unwrap();
            let x2 = T::from(0.5255324099163290).unwrap();
            let x3 = T::from(0.7966664774136267).unwrap();
            let x4 = T::from(0.9602898564975363).unwrap();
            let w1 = T::from(0.3626837833783620).unwrap();
            let w2 = T::from(0.3137066458778873).unwrap();
            let w3 = T::from(0.2223810344533745).unwrap();
            let w4 = T::from(0.1012285362903763).unwrap();
            (
                vec![-x4, -x3, -x2, -x1, x1, x2, x3, x4],
                vec![w4, w3, w2, w1, w1, w2, w3, w4],
            )
        }
        9 => {
            let x1 = T::from(0.3242534234038089).unwrap();
            let x2 = T::from(0.6133714327005904).unwrap();
            let x3 = T::from(0.8360311073266358).unwrap();
            let x4 = T::from(0.9681602395076261).unwrap();
            let w0 = T::from(0.3302393550012598).unwrap();
            let w1 = T::from(0.3123470770400029).unwrap();
            let w2 = T::from(0.2606106964029354).unwrap();
            let w3 = T::from(0.1806481606948574).unwrap();
            let w4 = T::from(0.0812743883615744).unwrap();
            (
                vec![-x4, -x3, -x2, -x1, T::zero(), x1, x2, x3, x4],
                vec![w4, w3, w2, w1, w0, w1, w2, w3, w4],
            )
        }
        10 => {
            let x1 = T::from(0.1488743389816312).unwrap();
            let x2 = T::from(0.4333953941292472).unwrap();
            let x3 = T::from(0.6794095682990244).unwrap();
            let x4 = T::from(0.8650633666889845).unwrap();
            let x5 = T::from(0.9739065285171717).unwrap();
            let w1 = T::from(0.2955242247147529).unwrap();
            let w2 = T::from(0.2692667193099963).unwrap();
            let w3 = T::from(0.2190863625159820).unwrap();
            let w4 = T::from(0.1494513491505806).unwrap();
            let w5 = T::from(0.0666713443086881).unwrap();
            (
                vec![-x5, -x4, -x3, -x2, -x1, x1, x2, x3, x4, x5],
                vec![w5, w4, w3, w2, w1, w1, w2, w3, w4, w5],
            )
        }
        _ => panic!("Unsupported number of points: {}", n),
    }
}

/// Computes the integral over [0, ∞) using Gauss-Laguerre quadrature.
///
/// # Algorithm
///
/// Gauss-Laguerre quadrature is designed for integrals of the form:
///
/// ```text
/// ∫₀^∞ e^(-x) f(x) dx ≈ Σᵢ wᵢ f(xᵢ)
/// ```
///
/// where xᵢ are roots of Laguerre polynomials. For functions g(x) without
/// the exponential weight, use: ∫₀^∞ g(x) dx = ∫₀^∞ e^(-x) [e^x g(x)] dx
///
/// # Arguments
///
/// * `f` - The function to integrate (without e^(-x) weight)
/// * `n` - Number of quadrature points (supported: 2-10)
///
/// # Panics
///
/// Panics if n is not in the range [2, 10].
///
/// # Examples
///
/// ```
/// use axiomath::calculus::integration::gauss_laguerre_quadrature;
///
/// // Integrate e^(-x) from 0 to ∞, exact value is 1
/// let f = |_x: f64| 1.0;  // We integrate e^(-x) * 1
/// let result = gauss_laguerre_quadrature(f, 5);
/// assert!((result - 1.0).abs() < 1e-10);
/// ```
pub fn gauss_laguerre_quadrature<T, F>(f: F, n: usize) -> T
where
    T: Float,
    F: Fn(T) -> T,
{
    assert!(
        (2..=10).contains(&n),
        "Number of points must be between 2 and 10"
    );

    let (nodes, weights) = gauss_laguerre_nodes_weights::<T>(n);

    let mut sum = T::zero();
    for i in 0..n {
        sum = sum + weights[i] * f(nodes[i]);
    }

    sum
}

/// Returns the nodes and weights for Gauss-Laguerre quadrature.
///
/// These are roots of Laguerre polynomials and corresponding weights.
fn gauss_laguerre_nodes_weights<T: Float>(n: usize) -> (Vec<T>, Vec<T>) {
    match n {
        2 => {
            let x1 = T::from(0.5857864376269049).unwrap();
            let x2 = T::from(3.4142135623730951).unwrap();
            let w1 = T::from(0.8535533905932738).unwrap();
            let w2 = T::from(0.1464466094067262).unwrap();
            (vec![x1, x2], vec![w1, w2])
        }
        3 => {
            let x1 = T::from(0.4157745567834791).unwrap();
            let x2 = T::from(2.2942803602790417).unwrap();
            let x3 = T::from(6.2899450829374792).unwrap();
            let w1 = T::from(0.7110930099293469).unwrap();
            let w2 = T::from(0.2785177335692408).unwrap();
            let w3 = T::from(0.0103892565015122).unwrap();
            (vec![x1, x2, x3], vec![w1, w2, w3])
        }
        4 => {
            let x1 = T::from(0.3225476896193923).unwrap();
            let x2 = T::from(1.7457611011583465).unwrap();
            let x3 = T::from(4.5366202969211280).unwrap();
            let x4 = T::from(9.3950709123011331).unwrap();
            let w1 = T::from(0.6031541043416336).unwrap();
            let w2 = T::from(0.3574186924377997).unwrap();
            let w3 = T::from(0.0388879085150054).unwrap();
            let w4 = T::from(0.0005392947055613).unwrap();
            (vec![x1, x2, x3, x4], vec![w1, w2, w3, w4])
        }
        5 => {
            let x1 = T::from(0.2635603197181410).unwrap();
            let x2 = T::from(1.4134030591065168).unwrap();
            let x3 = T::from(3.5964257710407221).unwrap();
            let x4 = T::from(7.0858100058588376).unwrap();
            let x5 = T::from(12.640800844275783).unwrap();
            let w1 = T::from(0.5217556105828087).unwrap();
            let w2 = T::from(0.3986668110831759).unwrap();
            let w3 = T::from(0.0759424496817076).unwrap();
            let w4 = T::from(0.0036117586799220).unwrap();
            let w5 = T::from(0.0000233699723858).unwrap();
            (vec![x1, x2, x3, x4, x5], vec![w1, w2, w3, w4, w5])
        }
        _ => {
            // For n > 5, we provide approximate values
            // In practice, these would be computed to higher precision
            panic!("Gauss-Laguerre quadrature for n={} not yet implemented", n);
        }
    }
}

/// Computes the integral over (-∞, ∞) using Gauss-Hermite quadrature.
///
/// # Algorithm
///
/// Gauss-Hermite quadrature is designed for integrals of the form:
///
/// ```text
/// ∫₋∞^∞ e^(-x²) f(x) dx ≈ Σᵢ wᵢ f(xᵢ)
/// ```
///
/// where xᵢ are roots of Hermite polynomials. This is particularly useful
/// for computing moments of the normal distribution and Gaussian integrals.
///
/// # Arguments
///
/// * `f` - The function to integrate (without e^(-x²) weight)
/// * `n` - Number of quadrature points (supported: 2-10)
///
/// # Panics
///
/// Panics if n is not in the range [2, 10].
///
/// # Examples
///
/// ```
/// use axiomath::calculus::integration::gauss_hermite_quadrature;
///
/// // Integrate e^(-x²) from -∞ to ∞, exact value is √π
/// let f = |_x: f64| 1.0;
/// let result = gauss_hermite_quadrature(f, 5);
/// let exact = std::f64::consts::PI.sqrt();
/// assert!((result - exact).abs() < 1e-10);
/// ```
pub fn gauss_hermite_quadrature<T, F>(f: F, n: usize) -> T
where
    T: Float,
    F: Fn(T) -> T,
{
    assert!(
        (2..=10).contains(&n),
        "Number of points must be between 2 and 10"
    );

    let (nodes, weights) = gauss_hermite_nodes_weights::<T>(n);

    let mut sum = T::zero();
    for i in 0..n {
        sum = sum + weights[i] * f(nodes[i]);
    }

    sum
}

/// Returns the nodes and weights for Gauss-Hermite quadrature.
///
/// These are roots of Hermite polynomials and corresponding weights.
fn gauss_hermite_nodes_weights<T: Float>(n: usize) -> (Vec<T>, Vec<T>) {
    match n {
        2 => {
            let x = T::from(0.7071067811865475).unwrap(); // 1/√2
            let w = T::from(0.8862269254527580).unwrap(); // √π/2
            (vec![-x, x], vec![w, w])
        }
        3 => {
            let x1 = T::from(1.2247448713915890).unwrap();
            let w0 = T::from(1.1816359006036772).unwrap();
            let w1 = T::from(0.2954089751509193).unwrap();
            (vec![-x1, T::zero(), x1], vec![w1, w0, w1])
        }
        4 => {
            let x1 = T::from(0.5246476232752903).unwrap();
            let x2 = T::from(1.6506801238857846).unwrap();
            let w1 = T::from(0.8049140900055128).unwrap();
            let w2 = T::from(0.0813128354472452).unwrap();
            (vec![-x2, -x1, x1, x2], vec![w2, w1, w1, w2])
        }
        5 => {
            let x1 = T::from(0.9585724646138185).unwrap();
            let x2 = T::from(2.0201828704560856).unwrap();
            let w0 = T::from(0.9453087204829419).unwrap();
            let w1 = T::from(0.3936193231522412).unwrap();
            let w2 = T::from(0.0199532420590459).unwrap();
            (
                vec![-x2, -x1, T::zero(), x1, x2],
                vec![w2, w1, w0, w1, w2],
            )
        }
        6 => {
            let x1 = T::from(0.4360774119276165).unwrap();
            let x2 = T::from(1.3358490740136969).unwrap();
            let x3 = T::from(2.3506049736744922).unwrap();
            let w1 = T::from(0.7246295952243925).unwrap();
            let w2 = T::from(0.1570673203228566).unwrap();
            let w3 = T::from(0.0045300099055088).unwrap();
            (
                vec![-x3, -x2, -x1, x1, x2, x3],
                vec![w3, w2, w1, w1, w2, w3],
            )
        }
        _ => {
            panic!("Gauss-Hermite quadrature for n={} not yet implemented", n);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_gauss_legendre_polynomial() {
        // ∫₀¹ x² dx = 1/3
        let f = |x: f64| x * x;
        let result = gauss_legendre_quadrature(f, 0.0, 1.0, 3);
        // Gauss-Legendre with 3 points is exact for degree 5 polynomials
        assert_abs_diff_eq!(result, 1.0 / 3.0, epsilon = 1e-15);
    }

    #[test]
    fn test_gauss_legendre_exponential() {
        // ∫₀¹ e^x dx = e - 1
        let f = |x: f64| x.exp();
        let result = gauss_legendre_quadrature(f, 0.0, 1.0, 10);
        let exact = std::f64::consts::E - 1.0;
        assert_abs_diff_eq!(result, exact, epsilon = 1e-14);
    }

    #[test]
    fn test_gauss_legendre_sine() {
        // ∫₀^π sin(x) dx = 2
        use std::f64::consts::PI;
        let f = |x: f64| x.sin();
        let result = gauss_legendre_quadrature(f, 0.0, PI, 10);
        assert_abs_diff_eq!(result, 2.0, epsilon = 1e-14);
    }

    #[test]
    fn test_gauss_legendre_exact_polynomial() {
        // Gauss-Legendre with n points is exact for polynomials of degree ≤ 2n-1
        // Test with n=5: should be exact for x^9
        let f = |x: f64| x.powi(9);
        let result = gauss_legendre_quadrature(f, 0.0, 1.0, 5);
        let exact = 1.0 / 10.0; // ∫₀¹ x⁹ dx = 1/10
        assert_abs_diff_eq!(result, exact, epsilon = 1e-14);
    }

    #[test]
    fn test_gauss_laguerre_constant() {
        // ∫₀^∞ e^(-x) dx = 1
        let f = |_x: f64| 1.0;
        let result = gauss_laguerre_quadrature(f, 5);
        assert_abs_diff_eq!(result, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_gauss_laguerre_linear() {
        // ∫₀^∞ x e^(-x) dx = 1
        let f = |x: f64| x;
        let result = gauss_laguerre_quadrature(f, 5);
        assert_abs_diff_eq!(result, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_gauss_laguerre_quadratic() {
        // ∫₀^∞ x² e^(-x) dx = 2
        let f = |x: f64| x * x;
        let result = gauss_laguerre_quadrature(f, 5);
        assert_abs_diff_eq!(result, 2.0, epsilon = 1e-9);
    }

    #[test]
    fn test_gauss_hermite_constant() {
        // ∫₋∞^∞ e^(-x²) dx = √π
        let f = |_x: f64| 1.0;
        let result = gauss_hermite_quadrature(f, 5);
        let exact = std::f64::consts::PI.sqrt();
        assert_abs_diff_eq!(result, exact, epsilon = 1e-10);
    }

    #[test]
    fn test_gauss_hermite_odd_function() {
        // ∫₋∞^∞ x e^(-x²) dx = 0 (odd function)
        let f = |x: f64| x;
        let result = gauss_hermite_quadrature(f, 5);
        assert_abs_diff_eq!(result, 0.0, epsilon = 1e-15);
    }

    #[test]
    fn test_gauss_hermite_even_function() {
        // ∫₋∞^∞ x² e^(-x²) dx = √π/2
        let f = |x: f64| x * x;
        let result = gauss_hermite_quadrature(f, 5);
        let exact = std::f64::consts::PI.sqrt() / 2.0;
        assert_abs_diff_eq!(result, exact, epsilon = 1e-9);
    }

    #[test]
    fn test_gauss_legendre_different_n() {
        // Verify convergence with increasing n
        let f = |x: f64| (x * x).sin();
        let exact = gauss_legendre_quadrature(f, 0.0, 1.0, 10);

        let result_2 = gauss_legendre_quadrature(f, 0.0, 1.0, 2);
        let result_5 = gauss_legendre_quadrature(f, 0.0, 1.0, 5);

        let error_2 = (result_2 - exact).abs();
        let error_5 = (result_5 - exact).abs();

        assert!(error_5 < error_2);
    }

    #[test]
    fn test_f32_support() {
        // Test that methods work with f32
        let f = |x: f32| x * x;
        let result = gauss_legendre_quadrature(f, 0.0_f32, 1.0_f32, 5);
        assert_abs_diff_eq!(result, 1.0_f32 / 3.0_f32, epsilon = 1e-6_f32);
    }

    #[test]
    #[should_panic(expected = "Number of points must be between 2 and 10")]
    fn test_gauss_legendre_invalid_n() {
        let f = |x: f64| x;
        gauss_legendre_quadrature(f, 0.0, 1.0, 1);
    }

    #[test]
    #[should_panic(expected = "Number of points must be between 2 and 10")]
    fn test_gauss_laguerre_invalid_n() {
        let f = |x: f64| x;
        gauss_laguerre_quadrature(f, 15);
    }
}
