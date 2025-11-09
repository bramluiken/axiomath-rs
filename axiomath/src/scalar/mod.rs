//! Scalar mathematical operations.
//!
//! This module implements fundamental mathematical operations on scalar values,
//! including trigonometric functions, exponential/logarithmic functions, and
//! other elementary functions, all from first principles.

pub mod exponential;
pub mod hyperbolic;
pub mod logarithmic;
pub mod power;
pub mod rounding;
pub mod trigonometric;

// Re-export all public functions for convenient access
pub use exponential::{exp10, exp2, exponential};
pub use hyperbolic::{
    hyperbolic_cosine, hyperbolic_sine, hyperbolic_tangent, inverse_hyperbolic_cosine,
    inverse_hyperbolic_sine, inverse_hyperbolic_tangent,
};
pub use logarithmic::{logarithm_10, logarithm_2, logarithm_base, natural_logarithm};
pub use power::{
    absolute_value, cube_root, nth_root, power, sign, square_root,
};
pub use rounding::{ceiling, clamp, floor, maximum, minimum, round, truncate};
pub use trigonometric::{
    arccosine, arcsine, arctangent, arctangent2, cosine, sine, tangent,
};
