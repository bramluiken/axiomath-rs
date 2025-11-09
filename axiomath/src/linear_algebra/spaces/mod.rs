//! Vector space operations.
//!
//! This module implements operations on vector spaces and subspaces, providing
//! tools for analyzing the structure and properties of vector collections.
//!
//! # Basis Operations
//!
//! Methods for working with bases and orthogonalization:
//!
//! - **Gram-Schmidt Orthogonalization**: Convert any basis to an orthonormal basis
//! - **Linear Independence**: Test if vectors are linearly independent
//! - **Span Dimension**: Compute the dimension of the span of a set of vectors
//! - **Basis Finding**: Extract a basis from a set of vectors
//!
//! # Fundamental Subspaces
//!
//! The four fundamental subspaces of a matrix:
//!
//! - **Null Space (Kernel)**: Solutions to Ax = 0
//! - **Column Space (Range/Image)**: Span of columns of A
//! - **Row Space**: Span of rows of A (column space of A^T)
//! - **Left Null Space**: Null space of A^T
//!
//! These subspaces satisfy the fundamental theorem of linear algebra:
//! - dim(column space) + dim(null space) = n
//! - dim(row space) + dim(left null space) = m
//! - dim(column space) = dim(row space) = rank(A)
//!
//! # References
//!
//! - Strang, Gilbert. "Introduction to Linear Algebra." 5th ed., 2016.
//! - Trefethen, Lloyd N., and David Bau III. "Numerical Linear Algebra." SIAM, 1997.

pub mod basis;
pub mod subspaces;

// Re-export main functions for convenience
pub use basis::{find_basis, gram_schmidt, is_linearly_independent, span_dimension};
pub use subspaces::{column_space, left_null_space, null_space, row_space};
