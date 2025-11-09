//! Mathematical and physical constants.
//!
//! This module provides high-precision constants used throughout the library.
//! All constants are provided in both f32 and f64 precision.

use crate::core::Float;

/// Mathematical constant π (pi).
///
/// The ratio of a circle's circumference to its diameter.
///
/// Value (f64): 3.14159265358979323846...
pub mod pi {
    /// π in f32 precision
    pub const F32: f32 = std::f32::consts::PI;
    /// π in f64 precision
    pub const F64: f64 = std::f64::consts::PI;
}

/// Mathematical constant τ (tau) = 2π.
///
/// The ratio of a circle's circumference to its radius.
/// Some mathematicians argue this is more natural than π.
///
/// Value (f64): 6.28318530717958647692...
pub mod tau {
    /// τ in f32 precision
    pub const F32: f32 = std::f32::consts::TAU;
    /// τ in f64 precision
    pub const F64: f64 = std::f64::consts::TAU;
}

/// Mathematical constant e (Euler's number).
///
/// The base of the natural logarithm.
///
/// Value (f64): 2.71828182845904523536...
pub mod e {
    /// e in f32 precision
    pub const F32: f32 = std::f32::consts::E;
    /// e in f64 precision
    pub const F64: f64 = std::f64::consts::E;
}

/// Golden ratio φ (phi).
///
/// φ = (1 + √5) / 2
///
/// Appears in many natural phenomena and is the most irrational number
/// (hardest to approximate with rational fractions).
///
/// Value (f64): 1.61803398874989484820...
pub mod golden_ratio {
    /// φ in f32 precision
    pub const F32: f32 = 1.618_033_988_749_894_8_f32;
    /// φ in f64 precision
    pub const F64: f64 = 1.618_033_988_749_894_848_204_586_834_365_638_117_72_f64;
}

/// Square root of 2.
///
/// √2 = 1.41421356237309504880...
///
/// The length of the diagonal of a unit square.
pub mod sqrt_2 {
    /// √2 in f32 precision
    pub const F32: f32 = std::f32::consts::SQRT_2;
    /// √2 in f64 precision
    pub const F64: f64 = std::f64::consts::SQRT_2;
}

/// Square root of 3.
///
/// √3 = 1.73205080756887729352...
pub mod sqrt_3 {
    /// √3 in f32 precision
    pub const F32: f32 = 1.732_050_807_568_877_3_f32;
    /// √3 in f64 precision
    pub const F64: f64 = 1.732_050_807_568_877_293_527_446_341_505_872_366_94_f64;
}

/// Natural logarithm of 2.
///
/// ln(2) = 0.69314718055994530941...
pub mod ln_2 {
    /// ln(2) in f32 precision
    pub const F32: f32 = std::f32::consts::LN_2;
    /// ln(2) in f64 precision
    pub const F64: f64 = std::f64::consts::LN_2;
}

/// Natural logarithm of 10.
///
/// ln(10) = 2.30258509299404568401...
pub mod ln_10 {
    /// ln(10) in f32 precision
    pub const F32: f32 = std::f32::consts::LN_10;
    /// ln(10) in f64 precision
    pub const F64: f64 = std::f64::consts::LN_10;
}

/// Euler-Mascheroni constant γ (gamma).
///
/// γ = lim(n→∞) (Σ(k=1 to n) 1/k - ln(n))
///
/// Value (f64): 0.57721566490153286060...
pub mod euler_mascheroni {
    /// γ in f32 precision
    pub const F32: f32 = 0.577_215_664_901_532_9_f32;
    /// γ in f64 precision
    pub const F64: f64 = 0.577_215_664_901_532_860_606_512_090_082_402_431_04_f64;
}

/// Physical constant: Boltzmann constant (k_B).
///
/// Relates the average kinetic energy of particles in a gas to the temperature.
///
/// Value: 1.380649 × 10⁻²³ J/K (exact since 2019 SI redefinition)
pub mod boltzmann {
    /// k_B in f32 precision (J/K)
    pub const F32: f32 = 1.380_649e-23_f32;
    /// k_B in f64 precision (J/K)
    pub const F64: f64 = 1.380_649e-23_f64;
}

/// Physical constant: Planck constant (h).
///
/// Fundamental constant in quantum mechanics relating energy to frequency.
///
/// Value: 6.62607015 × 10⁻³⁴ J⋅s (exact since 2019 SI redefinition)
pub mod planck {
    /// h in f32 precision (J⋅s)
    pub const F32: f32 = 6.626_070_15e-34_f32;
    /// h in f64 precision (J⋅s)
    pub const F64: f64 = 6.626_070_15e-34_f64;
}

/// Physical constant: Reduced Planck constant (ℏ = h / 2π).
///
/// Also called "h-bar", commonly used in quantum mechanics.
///
/// Value: 1.054571817... × 10⁻³⁴ J⋅s
pub mod planck_reduced {
    /// ℏ in f32 precision (J⋅s)
    pub const F32: f32 = 1.054_571_817e-34_f32;
    /// ℏ in f64 precision (J⋅s)
    pub const F64: f64 = 1.054_571_817e-34_f64;
}

/// Physical constant: Avogadro constant (N_A).
///
/// Number of constituent particles (atoms or molecules) per mole.
///
/// Value: 6.02214076 × 10²³ mol⁻¹ (exact since 2019 SI redefinition)
pub mod avogadro {
    /// N_A in f32 precision (mol⁻¹)
    pub const F32: f32 = 6.022_140_76e23_f32;
    /// N_A in f64 precision (mol⁻¹)
    pub const F64: f64 = 6.022_140_76e23_f64;
}

/// Physical constant: Gas constant (R = N_A × k_B).
///
/// Relates pressure, volume, temperature in ideal gas law.
///
/// Value: 8.314462618... J/(mol⋅K)
pub mod gas_constant {
    /// R in f32 precision (J/(mol⋅K))
    pub const F32: f32 = 8.314_462_618_f32;
    /// R in f64 precision (J/(mol⋅K))
    pub const F64: f64 = 8.314_462_618_153_24_f64;
}

/// Physical constant: Speed of light in vacuum (c).
///
/// Value: 299792458 m/s (exact, defines the meter)
pub mod speed_of_light {
    /// c in f32 precision (m/s)
    pub const F32: f32 = 299_792_458.0_f32;
    /// c in f64 precision (m/s)
    pub const F64: f64 = 299_792_458.0_f64;
}

/// Physical constant: Elementary charge (e).
///
/// The electric charge carried by a single proton.
///
/// Value: 1.602176634 × 10⁻¹⁹ C (exact since 2019 SI redefinition)
pub mod elementary_charge {
    /// e in f32 precision (C)
    pub const F32: f32 = 1.602_176_634e-19_f32;
    /// e in f64 precision (C)
    pub const F64: f64 = 1.602_176_634e-19_f64;
}

/// Physical constant: Electron mass (m_e).
///
/// Value: 9.1093837015 × 10⁻³¹ kg
pub mod electron_mass {
    /// m_e in f32 precision (kg)
    pub const F32: f32 = 9.109_383_701_5e-31_f32;
    /// m_e in f64 precision (kg)
    pub const F64: f64 = 9.109_383_701_5e-31_f64;
}

/// Physical constant: Proton mass (m_p).
///
/// Value: 1.67262192369 × 10⁻²⁷ kg
pub mod proton_mass {
    /// m_p in f32 precision (kg)
    pub const F32: f32 = 1.672_621_923_69e-27_f32;
    /// m_p in f64 precision (kg)
    pub const F64: f64 = 1.672_621_923_69e-27_f64;
}

/// Physical constant: Stefan-Boltzmann constant (σ).
///
/// Relates the power radiated by a black body to its temperature.
///
/// σ = π²k⁴/(60ℏ³c²)
///
/// Value: 5.670374419... × 10⁻⁸ W/(m²⋅K⁴)
pub mod stefan_boltzmann {
    /// σ in f32 precision (W/(m²⋅K⁴))
    pub const F32: f32 = 5.670_374_419e-8_f32;
    /// σ in f64 precision (W/(m²⋅K⁴))
    pub const F64: f64 = 5.670_374_419e-8_f64;
}

/// Physical constant: Gravitational constant (G).
///
/// Fundamental constant in Newton's law of universal gravitation.
///
/// Value: 6.67430 × 10⁻¹¹ m³/(kg⋅s²)
pub mod gravitational {
    /// G in f32 precision (m³/(kg⋅s²))
    pub const F32: f32 = 6.67430e-11_f32;
    /// G in f64 precision (m³/(kg⋅s²))
    pub const F64: f64 = 6.67430e-11_f64;
}

/// Generic constant retrieval for a given float type.
///
/// # Examples
///
/// ```
/// use axiomath::core::{Float, constants};
///
/// let pi_f32: f32 = constants::pi::<f32>();
/// let pi_f64: f64 = constants::pi::<f64>();
///
/// assert_eq!(pi_f32, std::f32::consts::PI);
/// assert_eq!(pi_f64, std::f64::consts::PI);
/// ```

/// Returns π for the given float type.
#[inline]
pub fn pi<T: Float>() -> T {
    T::from(std::f64::consts::PI).unwrap()
}

/// Returns τ (2π) for the given float type.
#[inline]
pub fn tau<T: Float>() -> T {
    T::from(std::f64::consts::TAU).unwrap()
}

/// Returns e (Euler's number) for the given float type.
#[inline]
pub fn e<T: Float>() -> T {
    T::from(std::f64::consts::E).unwrap()
}

/// Returns the golden ratio for the given float type.
#[inline]
pub fn golden_ratio<T: Float>() -> T {
    T::from(golden_ratio::F64).unwrap()
}

/// Returns √2 for the given float type.
#[inline]
pub fn sqrt_2<T: Float>() -> T {
    T::from(std::f64::consts::SQRT_2).unwrap()
}

/// Returns the Boltzmann constant for the given float type.
#[inline]
pub fn boltzmann<T: Float>() -> T {
    T::from(boltzmann::F64).unwrap()
}

/// Returns the Planck constant for the given float type.
#[inline]
pub fn planck<T: Float>() -> T {
    T::from(planck::F64).unwrap()
}

/// Returns the Avogadro constant for the given float type.
#[inline]
pub fn avogadro<T: Float>() -> T {
    T::from(avogadro::F64).unwrap()
}

/// Returns the gas constant for the given float type.
#[inline]
pub fn gas_constant<T: Float>() -> T {
    T::from(gas_constant::F64).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mathematical_constants() {
        assert!((pi::F64 - std::f64::consts::PI).abs() < 1e-15);
        assert!((e::F64 - std::f64::consts::E).abs() < 1e-15);
        assert_eq!(tau::F64, 2.0 * pi::F64);
    }

    #[test]
    fn test_generic_constant_retrieval() {
        let pi_f32: f32 = pi();
        let pi_f64: f64 = pi();

        assert_eq!(pi_f32, std::f32::consts::PI);
        assert_eq!(pi_f64, std::f64::consts::PI);
    }

    #[test]
    fn test_physical_constants_magnitude() {
        // Basic sanity checks for physical constants
        assert!(boltzmann::F64 > 0.0);
        assert!(boltzmann::F64 < 1e-20);

        assert!(planck::F64 > 0.0);
        assert!(planck::F64 < 1e-30);

        assert!(avogadro::F64 > 6e23);
        assert!(avogadro::F64 < 7e23);

        assert!(gas_constant::F64 > 8.0);
        assert!(gas_constant::F64 < 9.0);
    }

    #[test]
    fn test_derived_constants() {
        // Check that R ≈ N_A × k_B
        let r_calculated = avogadro::F64 * boltzmann::F64;
        let relative_error = (r_calculated - gas_constant::F64).abs() / gas_constant::F64;
        assert!(relative_error < 1e-6);
    }
}
