//! Linear system solvers.
//!
//! This module implements methods for solving linear systems of equations Ax = b
//! from first principles, using both direct and iterative approaches.
//!
//! # Direct Methods
//!
//! Direct methods compute the exact solution (up to round-off error) in a finite
//! number of steps:
//!
//! - **Gaussian Elimination**: Classical method with partial pivoting
//! - **LU Decomposition**: Factor once, solve many times
//! - **Triangular Solvers**: Forward and backward substitution
//!
//! # Matrix Inversion
//!
//! Methods for computing matrix inverses:
//!
//! - **Standard Inversion**: Via LU decomposition or Gauss-Jordan
//! - **Moore-Penrose Pseudoinverse**: For non-square or singular matrices (via SVD)
//!
//! # Iterative Methods
//!
//! Iterative methods are useful for large, sparse matrices:
//!
//! - **Jacobi Iteration**: Simple but slow convergence
//! - **Gauss-Seidel**: Improved convergence over Jacobi
//! - **Conjugate Gradient**: Fast convergence for symmetric positive definite matrices
//!
//! # References
//!
//! - Trefethen, Lloyd N., and David Bau III. "Numerical Linear Algebra." SIAM, 1997.
//! - Golub, Gene H., and Charles F. Van Loan. "Matrix Computations." 4th ed., 2013.

pub mod direct;
pub mod triangular;
pub mod inverse;
pub mod iterative;

// Re-export main functions for convenience
pub use direct::{gaussian_elimination, solve_linear_system};
pub use triangular::{solve_lower_triangular, solve_upper_triangular};
pub use inverse::{invert, pseudoinverse};
pub use iterative::{conjugate_gradient, gauss_seidel, jacobi};
