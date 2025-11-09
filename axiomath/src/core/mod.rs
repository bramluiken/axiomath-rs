//! Core traits and foundational types for axiomath.
//!
//! This module defines the fundamental traits that govern mathematical behavior
//! throughout the library. These traits extend and compose the traits from
//! `num-traits` to provide a coherent type system for mathematical operations.
//!
//! # Trait Hierarchy
//!
//! ```text
//! Scalar
//!   ├── Real
//!   │   └── Float
//!   └── Integer (future)
//! ```
//!
//! # Design Principles
//!
//! - **Modularity**: Traits are small and composable
//! - **Genericity**: Algorithms work over trait bounds, not concrete types
//! - **Zero-cost abstractions**: All trait methods inline to concrete implementations
//! - **Mathematical correctness**: Traits encode mathematical properties

pub mod traits;
pub mod constants;

pub use traits::*;
pub use constants::*;
