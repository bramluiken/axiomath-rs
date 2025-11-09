//! Prelude module providing commonly used types and traits.
//!
//! This module re-exports the most commonly used items from axiomath,
//! allowing users to quickly import everything they need with a single use statement.
//!
//! # Examples
//!
//! ```
//! use axiomath::prelude::*;
//!
//! // Core traits and constants
//! let x: f64 = 2.0;
//! let y = x.sqrt();
//!
//! // Vector operations
//! let v1 = vec3(1.0, 2.0, 3.0);
//! let v2 = vec3(4.0, 5.0, 6.0);
//! let dot = dot_product(&v1, &v2);
//! let mag = magnitude(&v1);
//! ```

// Core traits
pub use crate::core::{Float, Real, Scalar};
pub use crate::core::constants::*;

// Vector types and constructors
pub use crate::vector::{Vector, Vec2, Vec3, Vec4, vec2, vec3, vec4};

// Common vector operations
pub use crate::vector::{
    // Arithmetic
    add, subtract, multiply, divide, negate, reciprocal,
    // Scalar operations
    scale, divide_by_scalar,
    // Norms
    magnitude, magnitude_squared, normalize, distance, euclidean_norm,
    // Products
    dot_product, cross_product,
    // Angles
    angle_between,
};

// Re-export num-traits for convenience
pub use num_traits::{One, Zero};
