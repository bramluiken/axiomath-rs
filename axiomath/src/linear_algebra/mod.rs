//! Linear algebra algorithms.
//!
//! This module implements advanced linear algebra operations from first principles,
//! including matrix decompositions, eigenvalue computations, and linear system solvers.
//!
//! # Theory
//!
//! Linear algebra provides fundamental tools for solving systems of linear equations,
//! analyzing transformations, and understanding vector spaces. This module implements
//! key algorithms that form the foundation of numerical linear algebra.
//!
//! # Decompositions
//!
//! Matrix decompositions express a matrix as a product of matrices with special
//! properties. These decompositions are fundamental to:
//!
//! - Solving linear systems efficiently
//! - Computing matrix properties (rank, determinant, inverse)
//! - Numerical stability in algorithms
//! - Dimensionality reduction and data analysis
//!
//! Available decompositions:
//!
//! - **LU Decomposition**: A = PLU where P is a permutation matrix, L is lower
//!   triangular, U is upper triangular
//! - **QR Decomposition**: A = QR where Q is orthogonal, R is upper triangular
//! - **Cholesky Decomposition**: A = LL* for symmetric positive definite matrices
//! - **SVD**: A = UΣV* singular value decomposition
//! - **Eigenvalue Decomposition**: A = QΛQ⁻¹ for diagonalizable matrices
//!
//! # Linear Systems
//!
//! Methods for solving Ax = b:
//!
//! - Direct methods (Gaussian elimination, LU decomposition)
//! - Triangular system solvers (forward/backward substitution)
//! - Matrix inversion
//! - Iterative methods (Jacobi, Gauss-Seidel, Conjugate Gradient)
//!
//! # Vector Spaces
//!
//! Operations on vector spaces and subspaces:
//!
//! - Orthogonalization (Gram-Schmidt)
//! - Linear independence testing
//! - Basis computation
//! - Fundamental subspaces (null space, column space, row space)
//!
//! # Examples
//!
//! ```
//! use axiomath::prelude::*;
//! use axiomath::linear_algebra::decomposition::lu_decompose;
//! use axiomath::linear_algebra::systems::solve_linear_system;
//!
//! // LU decomposition
//! let a = Matrix::from_rows([
//!     [2.0, 1.0, 1.0],
//!     [4.0, 3.0, 3.0],
//!     [8.0, 7.0, 9.0],
//! ]);
//!
//! let (l, u, p) = lu_decompose(&a).unwrap();
//!
//! // Solve a linear system
//! let b = Vector::from([4.0, 10.0, 24.0]);
//! let x = solve_linear_system(&a, &b).unwrap();
//! ```
//!
//! # References
//!
//! - Trefethen, Lloyd N., and David Bau III. "Numerical Linear Algebra." SIAM, 1997.
//! - Golub, Gene H., and Charles F. Van Loan. "Matrix Computations." 4th ed., 2013.
//! - Strang, Gilbert. "Introduction to Linear Algebra." 5th ed., 2016.

use core::fmt;

/// Error types for linear algebra operations.
#[derive(Debug, Clone, PartialEq)]
pub enum LinearAlgebraError {
    /// Matrix is singular (determinant is zero or numerically close to zero).
    SingularMatrix,

    /// Matrix is not positive definite (required for Cholesky decomposition).
    NotPositiveDefinite,

    /// Matrix dimensions are incompatible for the operation.
    DimensionMismatch,

    /// Iterative algorithm failed to converge within the specified number of iterations.
    ConvergenceFailure {
        /// Number of iterations performed.
        iterations: usize,
    },

    /// Algorithm did not converge within maximum iterations.
    NotConverged,

    /// Numerical instability detected (e.g., condition number too large).
    NumericalInstability,

    /// Matrix is not square (required for certain operations).
    NotSquare,

    /// Matrix is not symmetric (required for certain decompositions).
    NotSymmetric,

    /// Feature or algorithm not yet implemented.
    NotImplemented(String),
}

impl fmt::Display for LinearAlgebraError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SingularMatrix => {
                write!(f, "Matrix is singular (determinant is zero or near zero)")
            }
            Self::NotPositiveDefinite => {
                write!(f, "Matrix is not positive definite")
            }
            Self::DimensionMismatch => {
                write!(f, "Matrix dimensions are incompatible")
            }
            Self::ConvergenceFailure { iterations } => {
                write!(
                    f,
                    "Iterative algorithm failed to converge after {} iterations",
                    iterations
                )
            }
            Self::NotConverged => {
                write!(f, "Algorithm did not converge within maximum iterations")
            }
            Self::NumericalInstability => {
                write!(f, "Numerical instability detected")
            }
            Self::NotSquare => {
                write!(f, "Matrix is not square")
            }
            Self::NotSymmetric => {
                write!(f, "Matrix is not symmetric")
            }
            Self::NotImplemented(msg) => {
                write!(f, "Not implemented: {}", msg)
            }
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for LinearAlgebraError {}

pub mod decomposition;
pub mod systems;
pub mod spaces;
pub mod row_operations;
pub mod rref;
