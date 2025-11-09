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
//! // Now you have access to common traits and types
//! let x: f64 = 2.0;
//! let y = x.sqrt();
//! ```

pub use crate::core::{Float, Real, Scalar};
pub use crate::core::constants::*;

// Re-export num-traits for convenience
pub use num_traits::{One, Zero};
