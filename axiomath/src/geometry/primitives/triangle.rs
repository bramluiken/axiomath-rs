//! Triangle primitive and related operations.

use crate::core::traits::{Real, Scalar};
use crate::vector::{cross_product, dot_product, magnitude, subtract, vec3, Vec3};
use num_traits::{Float, One, Zero};
use std::ops::{Add, Mul, Sub};

/// A triangle in 3D space represented by three vertices.
///
/// # Mathematical Background
///
/// A triangle is a polygon with three vertices and three edges. In 3D space,
/// three non-collinear points uniquely define a triangle and a plane containing it.
///
/// Key properties:
/// - **Area**: The magnitude of half the cross product of two edge vectors
/// - **Normal**: The unit vector perpendicular to the triangle's plane
/// - **Barycentric Coordinates**: A coordinate system where any point P in the
///   triangle's plane can be expressed as a weighted sum of the vertices
///
/// # Representation
///
/// This implementation stores:
/// - `p1`, `p2`, `p3`: The three vertices of the triangle
///
/// The vertices are ordered, and the normal direction follows the right-hand rule
/// based on the vertex order (p1 → p2 → p3).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Triangle<T>
where
    T: Scalar,
{
    /// First vertex of the triangle
    pub p1: Vec3<T>,
    /// Second vertex of the triangle
    pub p2: Vec3<T>,
    /// Third vertex of the triangle
    pub p3: Vec3<T>,
}

impl<T> Triangle<T>
where
    T: Scalar + Float,
{
    /// Creates a new triangle from three vertices.
    ///
    /// # Parameters
    ///
    /// * `p1` - The first vertex
    /// * `p2` - The second vertex
    /// * `p3` - The third vertex
    ///
    /// # Example
    ///
    /// ```
    /// use axiomath::geometry::primitives::Triangle;
    /// use axiomath::vector::vec3;
    ///
    /// let tri = Triangle::new(
    ///     vec3(0.0, 0.0, 0.0),
    ///     vec3(1.0, 0.0, 0.0),
    ///     vec3(0.0, 1.0, 0.0)
    /// );
    /// ```
    pub fn new(p1: Vec3<T>, p2: Vec3<T>, p3: Vec3<T>) -> Self {
        Triangle { p1, p2, p3 }
    }
}

/// Computes the area of a triangle given its three vertices.
///
/// # Mathematical Background
///
/// The area of a triangle with vertices P₁, P₂, P₃ is computed using the
/// cross product formula:
///
/// ```text
/// Area = (1/2) ||(P₂ - P₁) × (P₃ - P₁)||
/// ```
///
/// This formula derives from the fact that the magnitude of the cross product
/// of two vectors equals the area of the parallelogram they span, and a
/// triangle is half of that parallelogram.
///
/// # Algorithm
///
/// 1. Compute edge vector e₁ = P₂ - P₁
/// 2. Compute edge vector e₂ = P₃ - P₁
/// 3. Compute cross product n = e₁ × e₂
/// 4. Compute magnitude ||n||
/// 5. Return ||n|| / 2
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Parameters
///
/// * `p1` - The first vertex
/// * `p2` - The second vertex
/// * `p3` - The third vertex
///
/// # Returns
///
/// The area of the triangle as a scalar value (always non-negative).
/// Returns 0 if the points are collinear (degenerate triangle).
///
/// # Examples
///
/// ```
/// use axiomath::geometry::primitives::triangle_area;
/// use axiomath::vector::vec3;
///
/// // Right triangle with legs of length 1
/// let area = triangle_area(
///     vec3(0.0, 0.0, 0.0),
///     vec3(1.0, 0.0, 0.0),
///     vec3(0.0, 1.0, 0.0)
/// );
/// assert!((area - 0.5).abs() < 1e-10);
/// ```
pub fn triangle_area<T>(p1: Vec3<T>, p2: Vec3<T>, p3: Vec3<T>) -> T
where
    T: Scalar + Real + Float + Zero + Add<Output = T> + Mul<Output = T> + Sub<Output = T>,
{
    // Compute edge vectors
    let edge1 = subtract(&p2, &p1);
    let edge2 = subtract(&p3, &p1);

    // Compute cross product
    let cross = cross_product(&edge1, &edge2);

    // Area is half the magnitude of the cross product
    let mag = magnitude(&cross);
    mag / (T::one() + T::one()) // Divide by 2
}

/// Computes the unit normal vector of a triangle.
///
/// # Mathematical Background
///
/// The normal vector of a triangle is perpendicular to its plane. It is
/// computed using the cross product of two edge vectors:
///
/// ```text
/// n = normalize((P₂ - P₁) × (P₃ - P₁))
/// ```
///
/// The direction follows the right-hand rule: if fingers curl from P₁→P₂ to P₁→P₃,
/// the thumb points in the normal direction.
///
/// # Algorithm
///
/// 1. Compute edge vector e₁ = P₂ - P₁
/// 2. Compute edge vector e₂ = P₃ - P₁
/// 3. Compute cross product n = e₁ × e₂
/// 4. Normalize n to unit length
/// 5. Return normalized vector (or None if degenerate)
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Parameters
///
/// * `p1` - The first vertex
/// * `p2` - The second vertex
/// * `p3` - The third vertex
///
/// # Returns
///
/// The unit normal vector as `Vec3<T>`, or `None` if the triangle is degenerate
/// (points are collinear or coincident).
///
/// # Examples
///
/// ```
/// use axiomath::geometry::primitives::triangle_normal;
/// use axiomath::vector::vec3;
///
/// // Triangle in XY plane
/// let normal = triangle_normal(
///     vec3(0.0, 0.0, 0.0),
///     vec3(1.0, 0.0, 0.0),
///     vec3(0.0, 1.0, 0.0)
/// );
///
/// assert!(normal.is_some());
/// let n = normal.unwrap();
/// // Should point in +Z direction
/// assert!(n[2] > 0.99);
/// ```
pub fn triangle_normal<T>(p1: Vec3<T>, p2: Vec3<T>, p3: Vec3<T>) -> Option<Vec3<T>>
where
    T: Scalar + Real + Float + Zero + One + Add<Output = T> + Mul<Output = T> + Sub<Output = T>,
{
    // Compute edge vectors
    let edge1 = subtract(&p2, &p1);
    let edge2 = subtract(&p3, &p1);

    // Compute cross product
    let cross = cross_product(&edge1, &edge2);

    // Normalize (returns None if magnitude is too small)
    let mag = magnitude(&cross);
    if mag <= T::epsilon() * (T::one() + T::one() + T::one() + T::one()) {
        // Degenerate triangle (collinear or coincident points)
        return None;
    }

    Some(vec3(cross[0] / mag, cross[1] / mag, cross[2] / mag))
}

/// Computes the barycentric coordinates of a point with respect to a triangle.
///
/// # Mathematical Background
///
/// Barycentric coordinates express a point P as a weighted combination of the
/// triangle's vertices:
///
/// ```text
/// P = u·A + v·B + w·C
/// where u + v + w = 1
/// ```
///
/// For a point inside the triangle, all coordinates are in [0, 1].
/// Points outside have at least one coordinate outside this range.
///
/// # Algorithm
///
/// Uses the dot product method:
///
/// ```text
/// v0 = B - A,  v1 = C - A,  v2 = P - A
/// d00 = v0 · v0
/// d01 = v0 · v1
/// d11 = v1 · v1
/// d20 = v2 · v0
/// d21 = v2 · v1
/// denom = d00 * d11 - d01 * d01
/// v = (d11 * d20 - d01 * d21) / denom
/// w = (d00 * d21 - d01 * d20) / denom
/// u = 1 - v - w
/// ```
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Parameters
///
/// * `point` - The point to compute coordinates for
/// * `tri` - The triangle
///
/// # Returns
///
/// A tuple `(u, v, w)` representing the barycentric coordinates corresponding
/// to vertices (p1, p2, p3). These coordinates always sum to 1.
///
/// # Panics
///
/// Panics if the triangle is degenerate (area is zero or points are collinear).
///
/// # Examples
///
/// ```
/// use axiomath::geometry::primitives::{Triangle, barycentric_coordinates};
/// use axiomath::vector::vec3;
///
/// let tri = Triangle::new(
///     vec3(0.0, 0.0, 0.0),
///     vec3(1.0, 0.0, 0.0),
///     vec3(0.0, 1.0, 0.0)
/// );
///
/// // Center of triangle
/// let center = vec3(1.0/3.0, 1.0/3.0, 0.0);
/// let (u, v, w) = barycentric_coordinates(center, &tri);
///
/// // All coordinates should be approximately 1/3
/// assert!((u - 1.0/3.0).abs() < 1e-10);
/// assert!((v - 1.0/3.0).abs() < 1e-10);
/// assert!((w - 1.0/3.0).abs() < 1e-10);
/// ```
pub fn barycentric_coordinates<T>(point: Vec3<T>, tri: &Triangle<T>) -> (T, T, T)
where
    T: Scalar + Float + Zero + One + Add<Output = T> + Mul<Output = T> + Sub<Output = T>,
{
    // Compute vectors relative to p1
    let v0 = subtract(&tri.p2, &tri.p1);
    let v1 = subtract(&tri.p3, &tri.p1);
    let v2 = subtract(&point, &tri.p1);

    // Compute dot products
    let d00 = dot_product(&v0, &v0);
    let d01 = dot_product(&v0, &v1);
    let d11 = dot_product(&v1, &v1);
    let d20 = dot_product(&v2, &v0);
    let d21 = dot_product(&v2, &v1);

    // Compute denominator
    let denom = d00 * d11 - d01 * d01;

    assert!(
        Float::abs(denom) > T::epsilon(),
        "Triangle is degenerate (collinear or coincident points)"
    );

    let inv_denom = T::one() / denom;

    // Compute barycentric coordinates
    let v = (d11 * d20 - d01 * d21) * inv_denom;
    let w = (d00 * d21 - d01 * d20) * inv_denom;
    let u = T::one() - v - w;

    (u, v, w)
}

/// Determines if a point lies inside a triangle.
///
/// # Mathematical Background
///
/// A point P is inside triangle ABC if its barycentric coordinates (u, v, w)
/// all satisfy 0 ≤ coordinate ≤ 1. Since u + v + w = 1 by definition, it's
/// sufficient to check:
///
/// ```text
/// u ≥ 0  AND  v ≥ 0  AND  w ≥ 0
/// ```
///
/// Equivalently (and more efficiently):
/// ```text
/// v ≥ 0  AND  w ≥ 0  AND  (v + w) ≤ 1
/// ```
///
/// # Algorithm
///
/// 1. Compute barycentric coordinates (u, v, w)
/// 2. Check if all coordinates are in [0, 1]
/// 3. Return true if inside, false otherwise
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Parameters
///
/// * `point` - The point to test
/// * `tri` - The triangle
///
/// # Returns
///
/// `true` if the point is inside or on the boundary of the triangle,
/// `false` otherwise.
///
/// # Panics
///
/// Panics if the triangle is degenerate (collinear points).
///
/// # Examples
///
/// ```
/// use axiomath::geometry::primitives::{Triangle, point_in_triangle};
/// use axiomath::vector::vec3;
///
/// let tri = Triangle::new(
///     vec3(0.0, 0.0, 0.0),
///     vec3(1.0, 0.0, 0.0),
///     vec3(0.0, 1.0, 0.0)
/// );
///
/// // Point inside triangle
/// assert!(point_in_triangle(vec3(0.25, 0.25, 0.0), &tri));
///
/// // Point outside triangle
/// assert!(!point_in_triangle(vec3(2.0, 2.0, 0.0), &tri));
/// ```
pub fn point_in_triangle<T>(point: Vec3<T>, tri: &Triangle<T>) -> bool
where
    T: Scalar + Float + Zero + One + Add<Output = T> + Mul<Output = T> + Sub<Output = T>,
{
    let (u, v, w) = barycentric_coordinates(point, tri);

    // Small epsilon for numerical tolerance
    let epsilon = T::epsilon() * (T::one() + T::one() + T::one() + T::one()); // 4ε
    let neg_epsilon = -epsilon;

    // Check if all barycentric coordinates are in [0, 1] with tolerance
    u >= neg_epsilon && v >= neg_epsilon && w >= neg_epsilon
        && u <= T::one() + epsilon
        && v <= T::one() + epsilon
        && w <= T::one() + epsilon
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::Float;

    #[test]
    fn test_triangle_creation() {
        let tri = Triangle::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));
        assert_eq!(tri.p1, vec3(0.0, 0.0, 0.0));
        assert_eq!(tri.p2, vec3(1.0, 0.0, 0.0));
        assert_eq!(tri.p3, vec3(0.0, 1.0, 0.0));
    }

    #[test]
    fn test_triangle_area_right_triangle() {
        // Right triangle with legs of length 1
        let area = triangle_area(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));
        assert!(Float::abs(area - 0.5) < 1e-10);
    }

    #[test]
    fn test_triangle_area_equilateral() {
        // Equilateral triangle with side length 2
        let p1 = vec3(0.0, 0.0, 0.0);
        let p2 = vec3(2.0, 0.0, 0.0);
        let p3 = vec3(1.0, 3.0_f64.sqrt(), 0.0);

        let area = triangle_area(p1, p2, p3);
        let expected = 3.0_f64.sqrt(); // Area = (√3/4) * side² = (√3/4) * 4 = √3

        assert!(Float::abs(area - expected) < 1e-10);
    }

    #[test]
    fn test_triangle_area_zero_collinear() {
        // Collinear points should give zero area
        let area = triangle_area(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(2.0, 0.0, 0.0));
        assert!(Float::abs(area) < 1e-10);
    }

    #[test]
    fn test_triangle_area_3d() {
        // Triangle in 3D space
        let p1 = vec3(0.0, 0.0, 0.0);
        let p2 = vec3(1.0, 0.0, 0.0);
        let p3 = vec3(0.0, 1.0, 1.0);

        let area = triangle_area(p1, p2, p3);
        // Area should be (1/2) * ||cross([1,0,0], [0,1,1])|| = (1/2) * ||[0,-1,1]|| = (1/2)*√2
        let expected = 2.0_f64.sqrt() / 2.0;

        assert!(Float::abs(area - expected) < 1e-10);
    }

    #[test]
    fn test_triangle_normal_xy_plane() {
        // Triangle in XY plane should have normal in +Z direction
        let normal = triangle_normal(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));

        assert!(normal.is_some());
        let n = normal.unwrap();

        assert!(Float::abs(n[0]) < 1e-10);
        assert!(Float::abs(n[1]) < 1e-10);
        assert!(Float::abs(n[2] - 1.0) < 1e-10);
    }

    #[test]
    fn test_triangle_normal_unit_length() {
        let normal = triangle_normal(vec3(0.0, 0.0, 0.0), vec3(2.0, 0.0, 0.0), vec3(0.0, 3.0, 0.0));

        assert!(normal.is_some());
        let n = normal.unwrap();

        // Check that it's a unit vector
        let mag = magnitude(&n);
        assert!(Float::abs(mag - 1.0) < 1e-10);
    }

    #[test]
    fn test_triangle_normal_collinear_points() {
        // Collinear points should return None
        let normal = triangle_normal(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(2.0, 0.0, 0.0));
        assert!(normal.is_none());
    }

    #[test]
    fn test_triangle_normal_right_hand_rule() {
        // Verify right-hand rule: p1→p2 × p1→p3
        let normal = triangle_normal(vec3(0.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0), vec3(1.0, 0.0, 0.0));

        assert!(normal.is_some());
        let n = normal.unwrap();

        // Should point in -Z direction
        assert!(Float::abs(n[0]) < 1e-10);
        assert!(Float::abs(n[1]) < 1e-10);
        assert!(Float::abs(n[2] + 1.0) < 1e-10);
    }

    #[test]
    fn test_barycentric_coordinates_vertices() {
        let tri = Triangle::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));

        // At p1, should be (1, 0, 0)
        let (u, v, w) = barycentric_coordinates(tri.p1, &tri);
        assert!(Float::abs(u - 1.0) < 1e-10);
        assert!(Float::abs(v) < 1e-10);
        assert!(Float::abs(w) < 1e-10);

        // At p2, should be (0, 1, 0)
        let (u, v, w) = barycentric_coordinates(tri.p2, &tri);
        assert!(Float::abs(u) < 1e-10);
        assert!(Float::abs(v - 1.0) < 1e-10);
        assert!(Float::abs(w) < 1e-10);

        // At p3, should be (0, 0, 1)
        let (u, v, w) = barycentric_coordinates(tri.p3, &tri);
        assert!(Float::abs(u) < 1e-10);
        assert!(Float::abs(v) < 1e-10);
        assert!(Float::abs(w - 1.0) < 1e-10);
    }

    #[test]
    fn test_barycentric_coordinates_center() {
        let tri = Triangle::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));

        // Center of triangle
        let center = vec3(1.0 / 3.0, 1.0 / 3.0, 0.0);
        let (u, v, w) = barycentric_coordinates(center, &tri);

        // All should be 1/3
        assert!(Float::abs(u - 1.0 / 3.0) < 1e-10);
        assert!(Float::abs(v - 1.0 / 3.0) < 1e-10);
        assert!(Float::abs(w - 1.0 / 3.0) < 1e-10);
    }

    #[test]
    fn test_barycentric_coordinates_sum_to_one() {
        let tri = Triangle::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));

        let point = vec3(0.3, 0.4, 0.0);
        let (u, v, w) = barycentric_coordinates(point, &tri);

        // Sum should always be 1
        assert!(Float::abs(u + v + w - 1.0) < 1e-10);
    }

    #[test]
    fn test_barycentric_coordinates_edge_midpoint() {
        let tri = Triangle::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));

        // Midpoint of edge p1-p2
        let mid = vec3(0.5, 0.0, 0.0);
        let (u, v, w) = barycentric_coordinates(mid, &tri);

        assert!(Float::abs(u - 0.5) < 1e-10);
        assert!(Float::abs(v - 0.5) < 1e-10);
        assert!(Float::abs(w) < 1e-10);
    }

    #[test]
    fn test_point_in_triangle_inside() {
        let tri = Triangle::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));

        // Point clearly inside
        assert!(point_in_triangle(vec3(0.25, 0.25, 0.0), &tri));

        // Center
        assert!(point_in_triangle(vec3(1.0 / 3.0, 1.0 / 3.0, 0.0), &tri));
    }

    #[test]
    fn test_point_in_triangle_on_vertex() {
        let tri = Triangle::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));

        // On vertices
        assert!(point_in_triangle(tri.p1, &tri));
        assert!(point_in_triangle(tri.p2, &tri));
        assert!(point_in_triangle(tri.p3, &tri));
    }

    #[test]
    fn test_point_in_triangle_on_edge() {
        let tri = Triangle::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));

        // On edges
        assert!(point_in_triangle(vec3(0.5, 0.0, 0.0), &tri));
        assert!(point_in_triangle(vec3(0.0, 0.5, 0.0), &tri));
        assert!(point_in_triangle(vec3(0.5, 0.5, 0.0), &tri));
    }

    #[test]
    fn test_point_in_triangle_outside() {
        let tri = Triangle::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));

        // Points clearly outside
        assert!(!point_in_triangle(vec3(2.0, 2.0, 0.0), &tri));
        assert!(!point_in_triangle(vec3(-1.0, 0.0, 0.0), &tri));
        assert!(!point_in_triangle(vec3(0.0, -1.0, 0.0), &tri));
        assert!(!point_in_triangle(vec3(1.0, 1.0, 0.0), &tri));
    }

    #[test]
    fn test_point_in_triangle_3d() {
        // Triangle not in a principal plane
        let tri = Triangle::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(0.0, 1.0, 1.0));

        // Point in the triangle's plane and inside
        let point = vec3(0.25, 0.25, 0.25);
        assert!(point_in_triangle(point, &tri));

        // Point outside the triangle's plane
        let point_out = vec3(0.25, 0.25, 5.0);
        assert!(!point_in_triangle(point_out, &tri));
    }

    #[test]
    #[should_panic(expected = "Triangle is degenerate")]
    fn test_barycentric_coordinates_degenerate() {
        let tri = Triangle::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(2.0, 0.0, 0.0));
        let _ = barycentric_coordinates(vec3(0.5, 0.5, 0.0), &tri);
    }

    #[test]
    fn test_triangle_area_is_positive() {
        // Area should always be non-negative
        let tri1 = Triangle::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));
        let area1 = triangle_area(tri1.p1, tri1.p2, tri1.p3);
        assert!(area1 >= 0.0);

        // Reversed vertex order should give same area
        let tri2 = Triangle::new(vec3(0.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0), vec3(1.0, 0.0, 0.0));
        let area2 = triangle_area(tri2.p1, tri2.p2, tri2.p3);
        assert!(Float::abs(area1 - area2) < 1e-10);
    }
}
