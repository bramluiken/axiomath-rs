# Contributing to axiomath-rs

Thank you for your interest in contributing to axiomath-rs! This document provides guidelines and standards for implementing features in this first-principle mathematics library.

---

## Table of Contents

1. [Philosophy & Principles](#philosophy--principles)
2. [Code Style & Conventions](#code-style--conventions)
3. [Documentation Standards](#documentation-standards)
4. [Testing Requirements](#testing-requirements)
5. [Mathematical Implementation Guidelines](#mathematical-implementation-guidelines)
6. [Performance Considerations](#performance-considerations)
7. [API Design Guidelines](#api-design-guidelines)
8. [Pull Request Process](#pull-request-process)

---

## Philosophy & Principles

### First Principles Approach

This library implements mathematics from fundamental axioms and mathematical definitions, not by wrapping existing libraries.

**✅ DO:**
- Implement algorithms from mathematical definitions
- Use Taylor series, iterative methods, and analytical solutions
- Document the mathematical basis of each implementation
- Reference textbooks, papers, or mathematical proofs

**❌ DON'T:**
- Wrap external C libraries (except for `libm` in no_std)
- Use "black box" implementations without explanation
- Copy implementations without understanding them

### Educational Value

Code should teach mathematics while being production-ready.

**✅ DO:**
- Write clear, readable code that follows mathematical notation
- Use descriptive variable names (except for well-known mathematical symbols)
- Include extensive mathematical documentation
- Explain "why" not just "what"

**❌ DON'T:**
- Optimize prematurely at the cost of clarity
- Use obscure tricks without explanation
- Skip documentation because "it's obvious"

### Correctness First, Performance Second

**Priority Order:**
1. **Correctness** - Get the math right
2. **Numerical Stability** - Handle edge cases properly
3. **Clarity** - Keep code understandable
4. **Performance** - Optimize after profiling

---

## Code Style & Conventions

### Formatting

Use `rustfmt` with default settings:

```bash
cargo fmt
```

### Naming Conventions

#### Functions
Use descriptive, full names in snake_case:

**✅ GOOD:**
```rust
pub fn greatest_common_divisor(a: i64, b: i64) -> i64
pub fn matrix_multiply<T, const M: usize, const N: usize, const P: usize>(...)
pub fn fermi_dirac_distribution(energy: f64, chemical_potential: f64, temperature: f64) -> f64
```

**❌ BAD:**
```rust
pub fn gcd(a: i64, b: i64) -> i64  // Too abbreviated
pub fn mult(...)  // Unclear
pub fn fd(...)  // Cryptic abbreviation
```

**Exception:** Well-established mathematical abbreviations are acceptable:
```rust
pub fn erf(x: f64) -> f64  // Error function
pub fn gamma(x: f64) -> f64  // Gamma function
pub fn sin(x: f64) -> f64  // Sine function
```

#### Types
Use descriptive PascalCase names:

```rust
pub struct Vector<T, const N: usize> { ... }
pub struct Matrix<T, const ROWS: usize, const COLS: usize> { ... }
pub struct Quaternion<T> { ... }
```

#### Variables

**In mathematical contexts**, single-letter variables matching standard notation are preferred:

```rust
// Mathematical algorithm - OK to use standard notation
fn dot_product<T: Scalar, const N: usize>(a: &Vector<T, N>, b: &Vector<T, N>) -> T {
    let mut sum = T::zero();
    for i in 0..N {
        sum += a[i] * b[i];
    }
    sum
}

// Physics/thermodynamics - use standard symbols
fn boltzmann_distribution(E: f64, T: f64) -> f64 {
    let k_B = constants::boltzmann::F64;
    (-E / (k_B * T)).exp()
}
```

**In implementation code**, use descriptive names:

```rust
fn calculate_mean(data: &[f64]) -> f64 {
    let sum: f64 = data.iter().sum();
    let count = data.len() as f64;
    sum / count
}
```

### Module Organization

Each module should have:

```rust
//! Module-level documentation explaining the mathematical domain
//!
//! # Theory
//!
//! Mathematical background and theory
//!
//! # Examples
//!
//! ```
//! use axiomath::module::function;
//! // Example usage
//! ```

// Imports
use crate::core::Float;
use num_traits::Float as NumFloat;

// Public types and structs

// Public functions (sorted alphabetically or by logical grouping)

// Private helper functions

// Tests
#[cfg(test)]
mod tests {
    use super::*;
    // Unit tests
}
```

### Const Generics

Use const generics for compile-time dimension checking:

```rust
// ✅ GOOD - Compile-time safety
pub struct Matrix<T, const ROWS: usize, const COLS: usize> {
    data: [[T; COLS]; ROWS]
}

impl<T, const M: usize, const N: usize, const P: usize> Matrix<T, M, N> {
    pub fn multiply(&self, other: &Matrix<T, N, P>) -> Matrix<T, M, P> {
        // Type system guarantees dimensions match
    }
}

// ❌ BAD - Runtime dimension checking
pub struct Matrix<T> {
    data: Vec<Vec<T>>,
    rows: usize,
    cols: usize,
}
```

### Error Handling

For mathematical operations:

**✅ DO:**
- Return `Option<T>` for operations that may have no solution
- Return `Result<T, E>` for operations that can fail with specific errors
- Document panic conditions in docstrings
- Use `debug_assert!` for internal invariants

```rust
/// Computes the square root using Newton-Raphson method.
///
/// Returns `None` if the input is negative (for real numbers).
///
/// # Examples
/// ```
/// assert_eq!(square_root(4.0), Some(2.0));
/// assert_eq!(square_root(-1.0), None);
/// ```
pub fn square_root(x: f64) -> Option<f64> {
    if x < 0.0 {
        return None;
    }
    // Implementation
}
```

**❌ DON'T:**
- Use `unwrap()` in library code
- Panic without documentation
- Silently return incorrect values

---

## Documentation Standards

### Function Documentation

Every public function must have comprehensive documentation:

```rust
/// Computes the Nth Fibonacci number using matrix exponentiation.
///
/// This algorithm uses the matrix form of Fibonacci numbers to achieve
/// O(log n) time complexity through binary exponentiation.
///
/// # Mathematical Background
///
/// The Fibonacci sequence can be expressed as:
/// ```text
/// [F(n+1)]   [1 1]^n   [1]
/// [F(n)  ] = [1 0]   * [0]
/// ```
///
/// # Algorithm
///
/// 1. Express Fibonacci as matrix power
/// 2. Use binary exponentiation for efficiency
/// 3. Extract result from matrix
///
/// # Complexity
///
/// - Time: O(log n)
/// - Space: O(1)
///
/// # Arguments
///
/// * `n` - The index of the Fibonacci number to compute (0-indexed)
///
/// # Returns
///
/// The Nth Fibonacci number as a `u64`.
///
/// # Panics
///
/// Panics if the result would overflow `u64`.
///
/// # Examples
///
/// ```
/// use axiomath::number_theory::fibonacci;
///
/// assert_eq!(fibonacci(0), 0);
/// assert_eq!(fibonacci(1), 1);
/// assert_eq!(fibonacci(10), 55);
/// ```
///
/// # References
///
/// - Knuth, Donald. "The Art of Computer Programming, Vol 1", Section 1.2.8
pub fn fibonacci(n: u64) -> u64 {
    // Implementation
}
```

### Required Sections

1. **Summary** (first line) - One sentence describing what it does
2. **Mathematical Background** - Theory and definitions
3. **Algorithm** (if non-trivial) - Step-by-step explanation
4. **Complexity** - Time and space complexity
5. **Arguments** - Description of parameters
6. **Returns** - What the function returns
7. **Panics** - Any panic conditions
8. **Examples** - Working, testable examples
9. **References** - Academic sources (papers, books)

### Mathematical Notation

Use ASCII notation in documentation:

```rust
/// Computes the integral using Simpson's rule:
///
/// ```text
/// ∫[a,b] f(x)dx ≈ (b-a)/6 * [f(a) + 4f((a+b)/2) + f(b)]
/// ```
///
/// Or equivalently:
/// ```text
/// I ≈ h/3 * [f(x₀) + 4f(x₁) + 2f(x₂) + ... + 4f(xₙ₋₁) + f(xₙ)]
/// ```
/// where h = (b-a)/n
```

### Module Documentation

```rust
//! Linear algebra decomposition algorithms.
//!
//! This module implements various matrix decomposition methods from first principles,
//! including LU, QR, Cholesky, SVD, and eigenvalue decomposition.
//!
//! # Theory
//!
//! Matrix decompositions express a matrix as a product of matrices with special
//! properties. These decompositions are fundamental to:
//!
//! - Solving linear systems efficiently
//! - Computing matrix properties (rank, determinant)
//! - Numerical stability in algorithms
//!
//! # Available Decompositions
//!
//! - **LU Decomposition**: A = LU where L is lower triangular, U is upper triangular
//! - **QR Decomposition**: A = QR where Q is orthogonal, R is upper triangular
//! - **Cholesky Decomposition**: A = LL* for positive definite matrices
//! - **SVD**: A = UΣV* singular value decomposition
//! - **Eigenvalue Decomposition**: A = QΛQ⁻¹ for diagonalizable matrices
//!
//! # Examples
//!
//! ```
//! use axiomath::linear_algebra::decomposition::lu_decompose;
//! use axiomath::matrix::Mat3;
//!
//! let a = Mat3::new([
//!     [2.0, 1.0, 1.0],
//!     [4.0, 3.0, 3.0],
//!     [8.0, 7.0, 9.0],
//! ]);
//!
//! let (l, u, p) = lu_decompose(&a);
//! // Use decomposition...
//! ```
//!
//! # References
//!
//! - Trefethen, Lloyd N., and David Bau III. "Numerical Linear Algebra." SIAM, 1997.
//! - Golub, Gene H., and Charles F. Van Loan. "Matrix Computations." 4th ed., 2013.
```

---

## Testing Requirements

### Test Coverage

**Minimum Requirements:**
- ✅ **90%+ code coverage** for all modules
- ✅ **100% coverage** for core functionality
- ✅ All public functions must have tests

### Test Types

#### 1. Unit Tests

Test individual functions with known values:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::approximately_equal;

    #[test]
    fn test_factorial_base_cases() {
        assert_eq!(factorial(0), 1);
        assert_eq!(factorial(1), 1);
        assert_eq!(factorial(2), 2);
    }

    #[test]
    fn test_factorial_known_values() {
        assert_eq!(factorial(5), 120);
        assert_eq!(factorial(10), 3628800);
    }

    #[test]
    #[should_panic]
    fn test_factorial_overflow() {
        factorial(100); // Should panic for u64
    }
}
```

#### 2. Mathematical Identity Tests

Verify mathematical properties:

```rust
#[test]
fn test_trigonometric_identity() {
    use std::f64::consts::PI;

    for i in 0..100 {
        let x = (i as f64) * PI / 50.0;
        let sin_x = sine(x);
        let cos_x = cosine(x);

        // sin²(x) + cos²(x) = 1
        let result = sin_x * sin_x + cos_x * cos_x;
        assert!(approximately_equal(result, 1.0, 1e-10, 1e-10));
    }
}

#[test]
fn test_matrix_associativity() {
    // (AB)C = A(BC)
    let a = random_matrix::<3, 4>();
    let b = random_matrix::<4, 5>();
    let c = random_matrix::<5, 2>();

    let left = multiply(&multiply(&a, &b), &c);
    let right = multiply(&a, &multiply(&b, &c));

    assert_matrices_equal(&left, &right, 1e-10);
}
```

#### 3. Property-Based Tests

Use QuickCheck or Proptest:

```rust
#[cfg(test)]
mod property_tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn test_gcd_commutative(a in 1u64..1000, b in 1u64..1000) {
            // gcd(a, b) = gcd(b, a)
            assert_eq!(gcd(a, b), gcd(b, a));
        }

        #[test]
        fn test_gcd_divides_both(a in 1u64..1000, b in 1u64..1000) {
            let g = gcd(a, b);
            // gcd divides both numbers
            assert_eq!(a % g, 0);
            assert_eq!(b % g, 0);
        }

        #[test]
        fn test_vector_magnitude_positive(
            v in prop::array::uniform3(any::<f64>())
        ) {
            let vec = Vector::from(v);
            let mag = magnitude(&vec);
            // Magnitude is always non-negative
            assert!(mag >= 0.0);
        }
    }
}
```

#### 4. Numerical Accuracy Tests

Test against reference implementations or high-precision values:

```rust
#[test]
fn test_error_function_accuracy() {
    // Values from NIST Digital Library of Mathematical Functions
    let test_cases = [
        (0.0, 0.0),
        (0.5, 0.5204998778130465),
        (1.0, 0.8427007929497149),
        (2.0, 0.9953222650189527),
    ];

    for (x, expected) in test_cases {
        let result = error_function(x);
        let abs_error = (result - expected).abs();
        assert!(abs_error < 1e-14,
            "erf({}) = {}, expected {}, error = {}",
            x, result, expected, abs_error);
    }
}
```

#### 5. Edge Case Tests

Test boundary conditions:

```rust
#[test]
fn test_edge_cases() {
    // Zero
    assert_eq!(log_gamma(1.0), 0.0);

    // Very small
    let small = 1e-10;
    assert!(gamma(small).is_finite());

    // Very large
    let large = 100.0;
    assert!(log_gamma(large).is_finite());
    assert!(gamma(large).is_infinite()); // Expected overflow

    // Negative (should handle gracefully)
    assert!(gamma(-1.0).is_nan());

    // NaN propagation
    assert!(gamma(f64::NAN).is_nan());
}
```

### Benchmarks

Create benchmarks for performance-critical code:

```rust
// benches/matrix_ops.rs
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use axiomath::matrix::{Mat4, multiply};

fn benchmark_matrix_multiply(c: &mut Criterion) {
    let a = Mat4::identity();
    let b = Mat4::identity();

    c.bench_function("mat4_multiply", |bencher| {
        bencher.iter(|| {
            multiply(black_box(&a), black_box(&b))
        });
    });
}

criterion_group!(benches, benchmark_matrix_multiply);
criterion_main!(benches);
```

---

## Mathematical Implementation Guidelines

### Algorithm Selection

Choose algorithms based on:

1. **Numerical Stability**
2. **Accuracy Requirements**
3. **Performance Characteristics**
4. **Simplicity and Clarity**

#### Example: Square Root

```rust
/// Computes square root using Newton-Raphson iteration.
///
/// # Algorithm
///
/// Newton-Raphson for f(x) = x² - a:
/// ```text
/// xₙ₊₁ = (xₙ + a/xₙ) / 2
/// ```
///
/// Converges quadratically for good initial guesses.
///
/// # Complexity
///
/// - Iterations: O(log log(1/ε)) for accuracy ε
/// - Per iteration: O(1)
///
pub fn square_root(a: f64) -> Option<f64> {
    if a < 0.0 {
        return None;
    }
    if a == 0.0 {
        return Some(0.0);
    }

    // Initial guess: use floating-point bit manipulation for speed
    let mut x = a;

    // Newton-Raphson iterations
    let epsilon = 1e-15;
    let max_iterations = 50;

    for _ in 0..max_iterations {
        let x_new = 0.5 * (x + a / x);

        // Check convergence
        if (x_new - x).abs() < epsilon * x {
            return Some(x_new);
        }

        x = x_new;
    }

    Some(x) // Return best approximation
}
```

### Handling Special Values

Always handle IEEE 754 special values:

```rust
pub fn some_function(x: f64) -> f64 {
    // Handle special cases first
    if x.is_nan() {
        return f64::NAN;
    }
    if x.is_infinite() {
        return if x > 0.0 {
            f64::INFINITY
        } else {
            f64::NEG_INFINITY
        };
    }
    if x == 0.0 {
        return 0.0; // or appropriate value
    }

    // Main computation
    // ...
}
```

### Numerical Stability

Prefer numerically stable algorithms:

**✅ GOOD - Kahan summation:**
```rust
pub fn sum_kahan(values: &[f64]) -> f64 {
    let mut sum = 0.0;
    let mut compensation = 0.0;

    for &value in values {
        let y = value - compensation;
        let t = sum + y;
        compensation = (t - sum) - y;
        sum = t;
    }

    sum
}
```

**❌ BAD - Naive summation (accumulates error):**
```rust
pub fn sum_naive(values: &[f64]) -> f64 {
    values.iter().sum() // Can lose precision
}
```

### Range Reduction

Use range reduction for periodic functions:

```rust
pub fn sine(x: f64) -> f64 {
    use std::f64::consts::PI;

    // Range reduction to [-π, π]
    let x = x % (2.0 * PI);
    let x = if x > PI {
        x - 2.0 * PI
    } else if x < -PI {
        x + 2.0 * PI
    } else {
        x
    };

    // Now compute using Taylor series or other method
    // on reduced range for better accuracy
    sine_reduced_range(x)
}
```

---

## Performance Considerations

### Optimization Strategy

1. **Write correct code first**
2. **Write tests**
3. **Profile to find bottlenecks**
4. **Optimize hot paths only**
5. **Benchmark to verify improvement**

### When to Use Unsafe

Use `unsafe` only when:
- You can prove it's correct
- It provides significant performance benefit
- It's well-documented
- It's wrapped in safe API

**✅ GOOD:**
```rust
/// Accesses element without bounds checking.
///
/// # Safety
///
/// Caller must ensure `i < N`.
#[inline]
pub unsafe fn get_unchecked(&self, i: usize) -> &T {
    debug_assert!(i < N, "Index out of bounds");
    self.data.get_unchecked(i)
}
```

### Inline Annotations

```rust
// Use #[inline] for small functions called frequently
#[inline]
pub fn dot_product<T: Scalar>(a: &[T], b: &[T]) -> T {
    // Small, hot function
}

// Use #[inline(always)] sparingly
#[inline(always)]
pub fn multiply_scalar<T: Scalar>(a: T, b: T) -> T {
    a * b // Trivial operation
}

// Don't inline large functions
pub fn matrix_svd<T, const M: usize, const N: usize>(
    matrix: &Matrix<T, M, N>
) -> (Matrix<T, M, M>, Vector<T, N>, Matrix<T, N, N>) {
    // Complex algorithm - let compiler decide
}
```

### SIMD Considerations

```rust
#[cfg(feature = "simd")]
pub fn dot_product_simd(a: &[f64], b: &[f64]) -> f64 {
    use std::simd::*;
    // SIMD implementation
}

#[cfg(not(feature = "simd"))]
pub fn dot_product_simd(a: &[f64], b: &[f64]) -> f64 {
    // Fallback to scalar
    dot_product_scalar(a, b)
}
```

---

## API Design Guidelines

### Generic vs Concrete

Provide both:

**Concrete types for ergonomics:**
```rust
pub type Vec2<T> = Vector<T, 2>;
pub type Vec3<T> = Vector<T, 3>;
pub type Vec4<T> = Vector<T, 4>;

pub type Mat2<T> = Matrix<T, 2, 2>;
pub type Mat3<T> = Matrix<T, 3, 3>;
pub type Mat4<T> = Matrix<T, 4, 4>;
```

**Generic types for flexibility:**
```rust
pub struct Vector<T, const N: usize> {
    data: [T; N]
}
```

### Builder Pattern for Complex Types

```rust
pub struct SimulationBuilder {
    temperature: Option<f64>,
    pressure: Option<f64>,
    particle_count: Option<usize>,
    // ...
}

impl SimulationBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn temperature(mut self, temp: f64) -> Self {
        self.temperature = Some(temp);
        self
    }

    pub fn build(self) -> Result<Simulation, BuildError> {
        // Validate and construct
    }
}
```

### Free Functions vs Methods

**Use free functions for:**
- Operations on multiple types
- Mathematical operations
- Composable functions

```rust
pub fn dot_product<T, const N: usize>(a: &Vector<T, N>, b: &Vector<T, N>) -> T
pub fn cross_product(a: &Vec3<f64>, b: &Vec3<f64>) -> Vec3<f64>
pub fn matrix_multiply<T, const M: usize, const N: usize, const P: usize>(
    a: &Matrix<T, M, N>,
    b: &Matrix<T, N, P>
) -> Matrix<T, M, P>
```

**Use methods for:**
- Type-specific operations
- Mutations
- Properties

```rust
impl<T, const N: usize> Vector<T, N> {
    pub fn magnitude(&self) -> T
    pub fn normalize(&mut self)
    pub fn is_zero(&self) -> bool
}
```

---

## Pull Request Process

### Before Starting Implementation

**⚠️ CRITICAL: Ensure Dependent APIs are Established**

Before implementing any feature, verify that all dependent APIs are fully defined and stable:

**✅ DO:**
- Check that all dependencies listed in the ROADMAP are marked complete
- Verify that dependent modules have stable, documented APIs
- Ensure that dependent types, functions, and traits are finalized
- Wait for API reviews to complete before starting implementation

**❌ DON'T:**
- Start implementing features that depend on undefined APIs
- Make assumptions about how dependent APIs will work
- Begin work on phases before their dependencies are complete
- Create temporary/placeholder APIs expecting them to change later

**Example:**
```rust
// ❌ BAD: Starting Phase 7 (Calculus) before Phase 3 (Vector Math) is complete
// Don't implement gradient() if vector operations aren't finalized

// ✅ GOOD: Wait for Vector API to stabilize first
// Phase 3 complete → Vector<T, N> API is stable
// Now safe to implement gradient() using Vec3, dot_product, etc.
```

This ensures:
- No breaking changes when dependencies evolve
- Consistent API design across modules
- Reduced refactoring and rework
- Smooth integration between components

### Before Submitting

1. **Run tests**: `cargo test --all-features`
2. **Run clippy**: `cargo clippy -- -D warnings`
3. **Format code**: `cargo fmt`
4. **Update documentation**: Ensure all public items documented
5. **Add examples**: Include usage examples
6. **Update CHANGELOG**: Add entry for your changes
7. **Check benchmarks**: Run if performance-related

### PR Template

```markdown
## Description
Brief description of changes

## Motivation
Why is this change necessary?

## Implementation
Explain your approach and any design decisions

## Mathematical Background
Cite papers, textbooks, or mathematical proofs

## Testing
- [ ] Unit tests added
- [ ] Property tests added (if applicable)
- [ ] Benchmarks added (if performance-critical)
- [ ] Tested against reference values
- [ ] Edge cases tested

## Performance Impact
Describe any performance implications

## Breaking Changes
List any breaking API changes

## Checklist
- [ ] Code follows style guidelines
- [ ] Documentation added/updated
- [ ] Tests pass locally
- [ ] No clippy warnings
- [ ] Code formatted with rustfmt
- [ ] CHANGELOG updated
```

### Review Process

PRs will be reviewed for:

1. **Correctness** - Is the math right?
2. **Code Quality** - Is it readable and maintainable?
3. **Documentation** - Is it comprehensive?
4. **Tests** - Is coverage adequate?
5. **Performance** - Are there obvious improvements?
6. **API Design** - Is it consistent with existing code?

---

## Additional Resources

### Recommended Reading

**Numerical Methods:**
- "Numerical Recipes" by Press et al.
- "Numerical Linear Algebra" by Trefethen & Bau
- "Matrix Computations" by Golub & Van Loan

**Statistical Mechanics:**
- "Statistical Mechanics" by Pathria & Beale
- "Introduction to Modern Statistical Mechanics" by Chandler

**Rust:**
- The Rust Book: https://doc.rust-lang.org/book/
- Rust API Guidelines: https://rust-lang.github.io/api-guidelines/

### Tools

- **rustfmt**: Code formatting
- **clippy**: Linting
- **cargo-tarpaulin**: Code coverage
- **cargo-criterion**: Benchmarking
- **cargo-doc**: Documentation generation

---

## Questions?

Feel free to:
- Open an issue for questions
- Start a discussion in GitHub Discussions
- Check existing documentation and examples

**Thank you for contributing to axiomath-rs!**
