//! Matrix decomposition algorithms.
//!
//! This module implements various matrix decomposition methods from first principles.
//! Matrix decompositions are fundamental tools in numerical linear algebra, providing
//! efficient and numerically stable methods for solving linear systems, computing
//! matrix properties, and analyzing data.
//!
//! # Available Decompositions
//!
//! - **LU Decomposition** ([`lu_decompose`]): Factors A = PLU where P is a permutation
//!   matrix, L is lower triangular with ones on the diagonal, and U is upper triangular.
//!   Used for solving linear systems and computing determinants.
//!
//! - **QR Decomposition** ([`qr_decompose`]): Factors A = QR where Q is orthogonal
//!   (Q^T Q = I) and R is upper triangular. Used for least squares problems and
//!   eigenvalue computation.
//!
//! - **Cholesky Decomposition** ([`cholesky_decompose`]): Factors A = LL* where L is
//!   lower triangular. Only applicable to symmetric positive definite matrices, but
//!   about twice as fast as LU decomposition.
//!
//! - **Singular Value Decomposition (SVD)**: Factors A = UΣV* where U and V are
//!   orthogonal and Σ is diagonal with non-negative entries. The most general and
//!   powerful decomposition.
//!
//! - **Eigenvalue Decomposition**: Factors A = QΛQ⁻¹ where Λ is diagonal containing
//!   eigenvalues. For symmetric matrices, this simplifies to A = QΛQ^T with Q orthogonal.
//!
//! # References
//!
//! - Trefethen, Lloyd N., and David Bau III. "Numerical Linear Algebra." SIAM, 1997.
//! - Golub, Gene H., and Charles F. Van Loan. "Matrix Computations." 4th ed., 2013.

pub mod lu;
pub mod qr;
pub mod cholesky;
pub mod svd;
pub mod eigen;

// Re-export main functions for convenience
pub use lu::lu_decompose;
pub use qr::{qr_decompose, qr_decompose_gram_schmidt, qr_decompose_householder};
pub use cholesky::cholesky_decompose;
pub use svd::svd;
pub use eigen::{eigenvalues, eigenpairs, power_iteration};
