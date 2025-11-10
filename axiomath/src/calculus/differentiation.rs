//! Numerical differentiation methods.
//!
//! This module implements various finite difference approximations for computing
//! derivatives of functions. The methods differ in accuracy order and the number
//! of function evaluations required.
//!
//! # Theory
//!
//! Numerical differentiation approximates the derivative using finite differences.
//! The choice of method depends on the trade-off between accuracy and computational cost.
//!
//! ## Accuracy Orders
//!
//! - Forward/Backward difference: O(h) - First order
//! - Central difference: O(h²) - Second order
//! - Five-point stencil: O(h⁴) - Fourth order
//!
//! Smaller step sizes (h) generally improve accuracy until round-off errors dominate.
//!
//! # Examples
//!
//! ```
//! use axiomath::calculus::differentiation::{central_difference, forward_difference};
//!
//! // Derivative of x² at x=3 is 2x = 6
//! let f = |x: f64| x * x;
//! let df = central_difference(f, 3.0, 1e-5);
//! assert!((df - 6.0).abs() < 1e-8);
//!
//! // Derivative of e^x is e^x
//! let exp_f = |x: f64| x.exp();
//! let df_exp = forward_difference(exp_f, 1.0, 1e-5);
//! assert!((df_exp - std::f64::consts::E).abs() < 1e-3);
//! ```

use crate::core::Float;

/// Computes the forward difference approximation of the derivative.
///
/// # Algorithm
///
/// ```text
/// f'(x) ≈ [f(x+h) - f(x)] / h
/// ```
///
/// This is a first-order accurate method with error O(h).
///
/// # Complexity
///
/// - Function evaluations: 2
/// - Accuracy order: O(h)
///
/// # Arguments
///
/// * `f` - The function to differentiate
/// * `x` - The point at which to compute the derivative
/// * `h` - The step size (should be small, typically 1e-5 to 1e-8)
///
/// # Examples
///
/// ```
/// use axiomath::calculus::differentiation::forward_difference;
///
/// let f = |x: f64| x.powi(3);  // f(x) = x³
/// let df = forward_difference(f, 2.0, 1e-5);  // f'(2) = 3x² = 12
/// assert!((df - 12.0).abs() < 1e-3);
/// ```
#[inline]
pub fn forward_difference<T, F>(f: F, x: T, h: T) -> T
where
    T: Float,
    F: Fn(T) -> T,
{
    (f(x + h) - f(x)) / h
}

/// Computes the backward difference approximation of the derivative.
///
/// # Algorithm
///
/// ```text
/// f'(x) ≈ [f(x) - f(x-h)] / h
/// ```
///
/// This is a first-order accurate method with error O(h).
///
/// # Complexity
///
/// - Function evaluations: 2
/// - Accuracy order: O(h)
///
/// # Arguments
///
/// * `f` - The function to differentiate
/// * `x` - The point at which to compute the derivative
/// * `h` - The step size (should be small, typically 1e-5 to 1e-8)
///
/// # Examples
///
/// ```
/// use axiomath::calculus::differentiation::backward_difference;
///
/// let f = |x: f64| x.sin();
/// let df = backward_difference(f, 0.0, 1e-5);  // cos(0) = 1
/// assert!((df - 1.0).abs() < 1e-5);
/// ```
#[inline]
pub fn backward_difference<T, F>(f: F, x: T, h: T) -> T
where
    T: Float,
    F: Fn(T) -> T,
{
    (f(x) - f(x - h)) / h
}

/// Computes the central difference approximation of the derivative.
///
/// # Algorithm
///
/// ```text
/// f'(x) ≈ [f(x+h) - f(x-h)] / (2h)
/// ```
///
/// This is a second-order accurate method with error O(h²), making it
/// significantly more accurate than forward or backward differences for
/// the same step size.
///
/// # Complexity
///
/// - Function evaluations: 2
/// - Accuracy order: O(h²)
///
/// # Arguments
///
/// * `f` - The function to differentiate
/// * `x` - The point at which to compute the derivative
/// * `h` - The step size (should be small, typically 1e-5 to 1e-8)
///
/// # Examples
///
/// ```
/// use axiomath::calculus::differentiation::central_difference;
///
/// let f = |x: f64| x.exp();  // f(x) = e^x
/// let df = central_difference(f, 1.0, 1e-5);  // f'(1) = e^1
/// assert!((df - std::f64::consts::E).abs() < 1e-8);
/// ```
#[inline]
pub fn central_difference<T, F>(f: F, x: T, h: T) -> T
where
    T: Float,
    F: Fn(T) -> T,
{
    let two = T::one() + T::one();
    (f(x + h) - f(x - h)) / (two * h)
}

/// Computes a five-point stencil approximation of the derivative.
///
/// # Algorithm
///
/// ```text
/// f'(x) ≈ [-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)] / (12h)
/// ```
///
/// This is a fourth-order accurate method with error O(h⁴), providing
/// excellent accuracy at the cost of four function evaluations.
///
/// # Complexity
///
/// - Function evaluations: 4
/// - Accuracy order: O(h⁴)
///
/// # Arguments
///
/// * `f` - The function to differentiate
/// * `x` - The point at which to compute the derivative
/// * `h` - The step size (should be small, typically 1e-3 to 1e-5)
///
/// # Examples
///
/// ```
/// use axiomath::calculus::differentiation::five_point_stencil;
///
/// let f = |x: f64| x.sin();  // f(x) = sin(x)
/// let df = five_point_stencil(f, 0.0, 1e-3);  // f'(0) = cos(0) = 1
/// assert!((df - 1.0).abs() < 1e-10);
/// ```
#[inline]
pub fn five_point_stencil<T, F>(f: F, x: T, h: T) -> T
where
    T: Float,
    F: Fn(T) -> T,
{
    let two = T::one() + T::one();
    let eight = two * two * two;
    let twelve = eight + two + two;

    let f_plus_2h = f(x + two * h);
    let f_plus_h = f(x + h);
    let f_minus_h = f(x - h);
    let f_minus_2h = f(x - two * h);

    (-f_plus_2h + eight * f_plus_h - eight * f_minus_h + f_minus_2h) / (twelve * h)
}

/// Computes the derivative at a point using automatic step size selection.
///
/// This function uses the central difference method with a step size automatically
/// chosen based on the magnitude of x and machine epsilon. For most purposes,
/// this provides a good balance between accuracy and numerical stability.
///
/// # Algorithm
///
/// Step size is chosen as:
/// ```text
/// h = ε^(1/3) * max(|x|, 1)
/// ```
///
/// where ε is machine epsilon. This balances truncation error and round-off error.
///
/// # Arguments
///
/// * `f` - The function to differentiate
/// * `x` - The point at which to compute the derivative
///
/// # Examples
///
/// ```
/// use axiomath::calculus::differentiation::derivative_at_point;
///
/// let f = |x: f64| x.powi(4);  // f(x) = x⁴
/// let df = derivative_at_point(f, 2.0);  // f'(2) = 4x³ = 32
/// assert!((df - 32.0).abs() < 1e-6);
/// ```
pub fn derivative_at_point<T, F>(f: F, x: T) -> T
where
    T: Float,
    F: Fn(T) -> T,
{
    let h = optimal_step_size(x);
    central_difference(f, x, h)
}

/// Computes the second derivative using central differences.
///
/// # Algorithm
///
/// ```text
/// f''(x) ≈ [f(x+h) - 2f(x) + f(x-h)] / h²
/// ```
///
/// This is a second-order accurate method with error O(h²).
///
/// # Complexity
///
/// - Function evaluations: 3
/// - Accuracy order: O(h²)
///
/// # Arguments
///
/// * `f` - The function to differentiate
/// * `x` - The point at which to compute the second derivative
/// * `h` - The step size (should be small, typically 1e-4 to 1e-6)
///
/// # Examples
///
/// ```
/// use axiomath::calculus::differentiation::second_derivative;
///
/// let f = |x: f64| x.powi(3);  // f(x) = x³
/// let d2f = second_derivative(f, 2.0, 1e-4);  // f''(2) = 6x = 12
/// assert!((d2f - 12.0).abs() < 1e-4);
/// ```
#[inline]
pub fn second_derivative<T, F>(f: F, x: T, h: T) -> T
where
    T: Float,
    F: Fn(T) -> T,
{
    let two = T::one() + T::one();
    let f_plus = f(x + h);
    let f_center = f(x);
    let f_minus = f(x - h);

    (f_plus - two * f_center + f_minus) / (h * h)
}

/// Computes the nth derivative using repeated central differences.
///
/// # Algorithm
///
/// This function iteratively applies numerical differentiation. It uses
/// a direct formula for the nth derivative based on finite differences.
///
/// **Warning:** Accuracy degrades significantly for high orders (n > 3)
/// due to accumulation of numerical errors. For n > 2, consider using
/// dedicated higher-order methods or symbolic differentiation.
///
/// # Arguments
///
/// * `f` - The function to differentiate
/// * `x` - The point at which to compute the nth derivative
/// * `n` - The order of the derivative (n ≥ 1)
/// * `h` - The step size (should be small, typically 1e-3 for n=2, 1e-2 for n=3)
///
/// # Examples
///
/// ```
/// use axiomath::calculus::differentiation::nth_derivative;
///
/// let f = |x: f64| x.powi(4);  // f(x) = x⁴
/// let d3f = nth_derivative(f, 1.0, 3, 1e-2);  // f'''(1) = 24x = 24
/// assert!((d3f - 24.0).abs() < 1e-1);
/// ```
pub fn nth_derivative<T, F>(f: F, x: T, n: usize, h: T) -> T
where
    T: Float,
    F: Fn(T) -> T,
{
    if n == 0 {
        return f(x);
    }
    if n == 1 {
        return central_difference(&f, x, h);
    }
    if n == 2 {
        return second_derivative(&f, x, h);
    }

    // For n > 2, use finite difference approximation
    // f'''(x) ≈ [f(x+2h) - 2f(x+h) + 2f(x-h) - f(x-2h)] / (2h³)
    // f''''(x) ≈ [f(x+2h) - 4f(x+h) + 6f(x) - 4f(x-h) + f(x-2h)] / h⁴

    // For general n, use symmetric finite difference
    // This is a simplified implementation that works for small n
    let two = T::one() + T::one();

    match n {
        3 => {
            // Third derivative: central difference formula
            let f_plus_2h = f(x + two * h);
            let f_plus_h = f(x + h);
            let f_minus_h = f(x - h);
            let f_minus_2h = f(x - two * h);

            (f_plus_2h - two * f_plus_h + two * f_minus_h - f_minus_2h) / (two * h * h * h)
        }
        4 => {
            // Fourth derivative: central difference formula
            let f_plus_2h = f(x + two * h);
            let f_plus_h = f(x + h);
            let f_x = f(x);
            let f_minus_h = f(x - h);
            let f_minus_2h = f(x - two * h);

            let four = two * two;
            let six = four + two;
            (f_plus_2h - four * f_plus_h + six * f_x - four * f_minus_h + f_minus_2h) / (h * h * h * h)
        }
        _ => {
            // For n > 4, use a general finite difference approximation
            // This is approximate and accuracy decreases for higher orders
            panic!("nth_derivative for n > 4 is not supported due to numerical instability. Consider using symbolic differentiation for high-order derivatives.");
        }
    }
}

/// Computes an optimal step size for numerical differentiation.
///
/// The optimal step size balances truncation error and round-off error.
/// For central differences, this is approximately ε^(1/3) where ε is
/// machine epsilon.
///
/// # Algorithm
///
/// ```text
/// h = ε^(1/3) * max(|x|, 1)
/// ```
///
/// # Arguments
///
/// * `x` - The point at which to compute the derivative
///
/// # Returns
///
/// An optimal step size for numerical differentiation at x.
#[inline]
fn optimal_step_size<T>(x: T) -> T
where
    T: Float,
{
    let epsilon = T::epsilon();
    let one = T::one();
    let three = one + one + one;

    // h = ε^(1/3) * max(|x|, 1)
    let scale = x.abs().max(one);
    epsilon.powf(one / three) * scale
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_forward_difference_polynomial() {
        // f(x) = x², f'(x) = 2x
        let f = |x: f64| x * x;
        let df = forward_difference(f, 3.0, 1e-5);
        assert_abs_diff_eq!(df, 6.0, epsilon = 1e-3);
    }

    #[test]
    fn test_backward_difference_polynomial() {
        // f(x) = x³, f'(x) = 3x²
        let f = |x: f64| x.powi(3);
        let df = backward_difference(f, 2.0, 1e-5);
        assert_abs_diff_eq!(df, 12.0, epsilon = 1e-3);
    }

    #[test]
    fn test_central_difference_exponential() {
        // f(x) = e^x, f'(x) = e^x
        let f = |x: f64| x.exp();
        let df = central_difference(f, 1.0, 1e-5);
        assert_abs_diff_eq!(df, std::f64::consts::E, epsilon = 1e-8);
    }

    #[test]
    fn test_central_difference_sine() {
        // f(x) = sin(x), f'(x) = cos(x)
        let f = |x: f64| x.sin();
        let df = central_difference(f, 0.0, 1e-5);
        assert_abs_diff_eq!(df, 1.0, epsilon = 1e-8);

        let df_pi_2 = central_difference(f, std::f64::consts::FRAC_PI_2, 1e-5);
        assert_abs_diff_eq!(df_pi_2, 0.0, epsilon = 1e-8);
    }

    #[test]
    fn test_five_point_stencil_accuracy() {
        // f(x) = sin(x), f'(0) = cos(0) = 1
        let f = |x: f64| x.sin();
        let df = five_point_stencil(f, 0.0, 1e-3);
        // Five-point should be very accurate even with larger h
        assert_abs_diff_eq!(df, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_five_point_vs_central() {
        // Compare accuracy: five-point should be better
        let f = |x: f64| x.powi(5);  // f'(x) = 5x⁴
        let x: f64 = 1.5;
        let exact = 5.0 * x.powi(4);

        let central = central_difference(f, x, 1e-2);
        let five_point = five_point_stencil(f, x, 1e-2);

        let error_central = (central - exact).abs();
        let error_five_point = (five_point - exact).abs();

        assert!(error_five_point < error_central);
    }

    #[test]
    fn test_derivative_at_point() {
        // f(x) = x⁴, f'(x) = 4x³
        let f = |x: f64| x.powi(4);
        let df = derivative_at_point(f, 2.0);
        assert_abs_diff_eq!(df, 32.0, epsilon = 1e-6);
    }

    #[test]
    fn test_second_derivative() {
        // f(x) = x³, f''(x) = 6x
        let f = |x: f64| x.powi(3);
        let d2f = second_derivative(f, 2.0, 1e-4);
        assert_abs_diff_eq!(d2f, 12.0, epsilon = 1e-4);

        // f(x) = sin(x), f''(x) = -sin(x)
        let f_sin = |x: f64| x.sin();
        let d2f_sin = second_derivative(f_sin, std::f64::consts::FRAC_PI_2, 1e-4);
        assert_abs_diff_eq!(d2f_sin, -1.0, epsilon = 1e-4);
    }

    #[test]
    fn test_nth_derivative_first_order() {
        // First derivative should match central difference
        let f = |x: f64| x.powi(3);
        let df1 = nth_derivative(f, 2.0, 1, 1e-5);
        let df_central = central_difference(f, 2.0, 1e-5);
        assert_abs_diff_eq!(df1, df_central, epsilon = 1e-10);
    }

    #[test]
    fn test_nth_derivative_second_order() {
        // f(x) = x⁴, f''(x) = 12x²
        let f = |x: f64| x.powi(4);
        let d2f = nth_derivative(f, 2.0, 2, 1e-4);
        assert_abs_diff_eq!(d2f, 48.0, epsilon = 1e-2);
    }

    #[test]
    fn test_nth_derivative_third_order() {
        // f(x) = x⁴, f'''(x) = 24x
        let f = |x: f64| x.powi(4);
        let d3f = nth_derivative(f, 1.0, 3, 1e-3);
        assert_abs_diff_eq!(d3f, 24.0, epsilon = 1e-1);
    }

    #[test]
    fn test_zero_derivative() {
        // f(x) = 5 (constant), f'(x) = 0
        let f = |_x: f64| 5.0;
        let df = central_difference(f, 10.0, 1e-5);
        assert_abs_diff_eq!(df, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_linear_function() {
        // f(x) = 3x + 2, f'(x) = 3
        let f = |x: f64| 3.0 * x + 2.0;
        let df = central_difference(f, 5.0, 1e-5);
        assert_abs_diff_eq!(df, 3.0, epsilon = 1e-8);
    }

    #[test]
    fn test_f32_support() {
        // Test that methods work with f32
        let f = |x: f32| x * x;
        let df = central_difference(f, 3.0_f32, 1e-4_f32);
        assert_abs_diff_eq!(df, 6.0_f32, epsilon = 1e-2_f32);
    }

    #[test]
    fn test_mathematical_identity() {
        // Product rule: (fg)' = f'g + fg'
        // f(x) = x², g(x) = sin(x)
        let f = |x: f64| x * x;
        let g = |x: f64| x.sin();
        let fg = |x: f64| f(x) * g(x);

        let x = 1.0;
        let h = 1e-5;

        let df = central_difference(f, x, h);  // 2x = 2
        let dg = central_difference(g, x, h);  // cos(x)
        let d_fg = central_difference(fg, x, h);

        let expected = df * g(x) + f(x) * dg;
        assert_abs_diff_eq!(d_fg, expected, epsilon = 1e-6);
    }
}
