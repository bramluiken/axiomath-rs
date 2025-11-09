//! # axiomath
//!
//! A comprehensive first-principle mathematics library for Rust.
//!
//! This library implements mathematical operations from fundamental axioms,
//! providing a complete suite of mathematical tools ranging from basic scalar
//! operations to advanced statistical mechanics and thermodynamics.
//!
//! ## Design Philosophy
//!
//! - **First Principles**: All algorithms implemented from mathematical foundations
//! - **Educational**: Clear, readable code with extensive documentation
//! - **Generic**: Works with f32, f64, and other numeric types
//! - **No-std Compatible**: Can be used in embedded environments
//! - **Type-Safe**: Leverages Rust's type system for compile-time guarantees
//!
//! ## Features
//!
//! - `std` (default): Enable standard library support
//! - `libm`: Use libm for no_std math functions
//! - `serde`: Enable serialization support
//! - `parallel`: Enable parallel computation via rayon
//! - `random`: Enable random number generation
//! - `simd`: Enable SIMD optimizations (requires nightly)
//!
//! ## Module Organization
//!
//! - [`core`]: Fundamental traits and types
//! - [`scalar`]: Scalar mathematical operations
//! - [`vector`]: Vector types and operations
//! - [`matrix`]: Matrix types and operations
//! - [`linear_algebra`]: Linear algebra algorithms
//! - [`geometry`]: Geometric operations and transformations
//! - [`calculus`]: Differentiation and integration
//! - [`statistics`]: Statistical analysis and probability
//! - [`number_theory`]: Number theoretic functions
//! - [`complex`]: Complex number arithmetic
//! - [`polynomial`]: Polynomial operations
//! - [`special`]: Special mathematical functions
//! - [`thermodynamics`]: Statistical mechanics and thermodynamics
//! - [`optimization`]: Optimization algorithms
//! - [`utils`]: Utility functions and constants

#![cfg_attr(not(feature = "std"), no_std)]
// SIMD support will be added in Phase 15
// #![cfg_attr(feature = "simd", feature(portable_simd))]
#![warn(missing_docs)]
#![warn(rustdoc::missing_crate_level_docs)]

#[cfg(not(feature = "std"))]
extern crate core as std;

// Re-export num-traits for convenience
pub use num_traits;

// Core modules
pub mod core;
pub mod prelude;
pub mod utils;

// Basic mathematics
pub mod scalar;
pub mod vector;
pub mod matrix;

// Advanced mathematics
pub mod linear_algebra;
pub mod geometry;
pub mod calculus;
pub mod statistics;
pub mod number_theory;
pub mod complex;
pub mod polynomial;
pub mod special;

// Physics and thermodynamics
pub mod thermodynamics;

// Computational tools
pub mod optimization;
