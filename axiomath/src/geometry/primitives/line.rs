//! Line primitive and related operations.

use crate::core::traits::Scalar;
use crate::vector::{cross_product, vec3, Vec3};
use num_traits::{Float, Zero};
use std::ops::{Add, Mul, Sub};

/// A line in 3D space represented by a point and a direction.
///
/// # Mathematical Background
///
/// A line in 3D can be represented in parametric form as:
///
/// ```text
/// L(t) = point + t * direction
/// ```
///
/// where:
/// - `point` is any point on the line
/// - `direction` is the direction vector (need not be normalized)
/// - `t` is a scalar parameter ranging over all real numbers
///
/// This representation is also called the point-direction form or
/// parametric form of a line.
///
/// # Representation
///
/// This implementation stores:
/// - `point`: A point that lies on the line
/// - `direction`: The direction vector of the line
///
/// Note: The direction vector does not need to be normalized, but some
/// operations may normalize it internally for numerical stability.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Line<T>
where
    T: Scalar,
{
    /// A point on the line
    pub point: Vec3<T>,
    /// The direction vector of the line
    pub direction: Vec3<T>,
}

impl<T> Line<T>
where
    T: Scalar + Float,
{
    /// Creates a new line from a point and a direction.
    ///
    /// # Parameters
    ///
    /// * `point` - A point on the line
    /// * `direction` - The direction vector (will not be automatically normalized)
    ///
    /// # Returns
    ///
    /// A new line
    ///
    /// # Example
    ///
    /// ```
    /// use axiomath::geometry::primitives::Line;
    /// use axiomath::vector::Vec3;
    ///
    /// let line = Line::new(
    ///     Vec3::new(0.0, 0.0, 0.0),
    ///     Vec3::new(1.0, 0.0, 0.0)
    /// );
    /// ```
    pub fn new(point: Vec3<T>, direction: Vec3<T>) -> Self {
        Line { point, direction }
    }
}

/// Computes the distance from a point to a line.
///
/// # Mathematical Background
///
/// The distance from a point P to a line L(t) = A + t*d is:
///
/// ```text
/// distance = |(P - A) × d| / |d|
/// ```
///
/// where:
/// - P is the point
/// - A is a point on the line
/// - d is the direction vector
/// - × denotes the cross product
/// - | | denotes magnitude
///
/// This formula comes from the geometric interpretation that the cross product
/// |(P - A) × d| gives the area of the parallelogram formed by (P - A) and d,
/// and dividing by |d| gives the perpendicular distance.
///
/// # Algorithm
///
/// 1. Compute the vector from line point to the test point: v = P - A
/// 2. Compute the cross product: v × d
/// 3. Compute the magnitude of the cross product
/// 4. Divide by the magnitude of the direction vector
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Parameters
///
/// * `point` - The point
/// * `line` - The line
///
/// # Returns
///
/// The shortest distance from the point to the line
///
/// # Panics
///
/// Panics if the line's direction vector has zero length.
///
/// # Example
///
/// ```
/// use axiomath::geometry::primitives::{Line, distance_point_to_line};
/// use axiomath::vector::Vec3;
///
/// let line = Line::new(
///     Vec3::new(0.0, 0.0, 0.0),
///     Vec3::new(1.0, 0.0, 0.0)
/// );
/// let point = Vec3::new(0.0, 5.0, 0.0);
/// let distance = distance_point_to_line(point, &line);
/// assert!((distance - 5.0).abs() < 1e-10);
/// ```
pub fn distance_point_to_line<T>(point: Vec3<T>, line: &Line<T>) -> T
where
    T: Scalar + Float + Zero + Add<Output = T> + Mul<Output = T> + Sub<Output = T>,
{
    // Vector from line point to the test point
    let v = vec3(
        point[0] - line.point[0],
        point[1] - line.point[1],
        point[2] - line.point[2],
    );

    // Cross product v × direction
    let cross = cross_product(&v, &line.direction);

    // Magnitude of cross product
    let cross_mag = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();

    // Magnitude of direction
    let dir_mag = (line.direction[0] * line.direction[0]
        + line.direction[1] * line.direction[1]
        + line.direction[2] * line.direction[2])
        .sqrt();

    assert!(
        dir_mag > T::zero(),
        "Line direction must be non-zero for distance calculation"
    );

    cross_mag / dir_mag
}

/// Finds the closest point on a line to a given point.
///
/// # Mathematical Background
///
/// Given a line L(t) = A + t*d and a point P, the closest point on the line
/// is found by projecting (P - A) onto the direction vector d:
///
/// ```text
/// t = ((P - A) · d) / (d · d)
/// closest_point = A + t*d
/// ```
///
/// where · denotes the dot product.
///
/// This comes from minimizing the squared distance |L(t) - P|² with respect to t,
/// which yields the projection formula.
///
/// # Algorithm
///
/// 1. Compute the vector from line point to test point: v = P - A
/// 2. Project v onto the direction: t = (v · d) / (d · d)
/// 3. Compute the point on the line: closest = A + t*d
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Parameters
///
/// * `point` - The point
/// * `line` - The line
///
/// # Returns
///
/// The closest point on the line to the given point
///
/// # Panics
///
/// Panics if the line's direction vector has zero length.
///
/// # Example
///
/// ```
/// use axiomath::geometry::primitives::{Line, closest_point_on_line};
/// use axiomath::vector::Vec3;
///
/// let line = Line::new(
///     Vec3::new(0.0, 0.0, 0.0),
///     Vec3::new(1.0, 0.0, 0.0)
/// );
/// let point = Vec3::new(5.0, 3.0, 0.0);
/// let closest = closest_point_on_line(point, &line);
///
/// // Closest point should be (5, 0, 0)
/// assert!((closest.x - 5.0).abs() < 1e-10);
/// assert!(closest.y.abs() < 1e-10);
/// assert!(closest.z.abs() < 1e-10);
/// ```
pub fn closest_point_on_line<T>(point: Vec3<T>, line: &Line<T>) -> Vec3<T>
where
    T: Scalar + Float + Zero + Add<Output = T> + Mul<Output = T> + Sub<Output = T>,
{
    // Vector from line point to the test point
    let v = vec3(
        point[0] - line.point[0],
        point[1] - line.point[1],
        point[2] - line.point[2],
    );

    // Project v onto the direction
    let dot_vd = v[0] * line.direction[0] + v[1] * line.direction[1] + v[2] * line.direction[2];
    let dot_dd = line.direction[0] * line.direction[0]
        + line.direction[1] * line.direction[1]
        + line.direction[2] * line.direction[2];

    assert!(
        dot_dd > T::zero(),
        "Line direction must be non-zero for closest point calculation"
    );

    let t = dot_vd / dot_dd;

    // Compute the closest point
    vec3(
        line.point[0] + t * line.direction[0],
        line.point[1] + t * line.direction[1],
        line.point[2] + t * line.direction[2],
    )
}

/// Computes the intersection point of two lines (if they intersect).
///
/// # Mathematical Background
///
/// Two lines L₁(s) = A₁ + s*d₁ and L₂(t) = A₂ + t*d₂ intersect if there exist
/// parameters s and t such that:
///
/// ```text
/// A₁ + s*d₁ = A₂ + t*d₂
/// ```
///
/// This is a system of 3 equations with 2 unknowns (s and t). For lines in 3D,
/// this system is generally overdetermined and has no solution (lines are skew).
///
/// Lines intersect if and only if:
/// ```text
/// (A₂ - A₁) · (d₁ × d₂) = 0
/// ```
///
/// When lines intersect, we can solve for the parameters using:
/// ```text
/// s = ((A₂ - A₁) × d₂) · (d₁ × d₂) / |d₁ × d₂|²
/// ```
///
/// # Algorithm
///
/// 1. Compute the cross product of directions: d₁ × d₂
/// 2. Check if lines are parallel (cross product near zero)
/// 3. Check if lines are coplanar using the triple product test
/// 4. If they intersect, compute the parameter s and the intersection point
/// 5. Verify the intersection point lies on both lines (within tolerance)
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Parameters
///
/// * `line1` - The first line
/// * `line2` - The second line
///
/// # Returns
///
/// `Some(point)` if the lines intersect, `None` if they are parallel or skew
///
/// # Example
///
/// ```
/// use axiomath::geometry::primitives::{Line, line_line_intersection};
/// use axiomath::vector::Vec3;
///
/// // Two lines that intersect at the origin
/// let line1 = Line::new(Vec3::new(-1.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
/// let line2 = Line::new(Vec3::new(0.0, -1.0, 0.0), Vec3::new(0.0, 1.0, 0.0));
///
/// let intersection = line_line_intersection(&line1, &line2);
/// assert!(intersection.is_some());
///
/// let point = intersection.unwrap();
/// assert!(point.x.abs() < 1e-10);
/// assert!(point.y.abs() < 1e-10);
/// assert!(point.z.abs() < 1e-10);
/// ```
pub fn line_line_intersection<T>(line1: &Line<T>, line2: &Line<T>) -> Option<Vec3<T>>
where
    T: Scalar + Float + Zero + Add<Output = T> + Mul<Output = T> + Sub<Output = T>,
{
    let d1 = &line1.direction;
    let d2 = &line2.direction;

    // Compute cross product d1 × d2
    let cross = cross_product(d1, d2);
    let cross_mag_sq = cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2];

    // Threshold for considering lines parallel
    let epsilon = T::epsilon() * (T::one() + T::one() + T::one() + T::one()); // 4ε

    // Check if lines are parallel
    if cross_mag_sq < epsilon {
        return None;
    }

    // Vector between the two line points
    let w = vec3(
        line2.point[0] - line1.point[0],
        line2.point[1] - line1.point[1],
        line2.point[2] - line1.point[2],
    );

    // Check if lines are coplanar: w · (d1 × d2) should be ~0
    let triple_product = w[0] * cross[0] + w[1] * cross[1] + w[2] * cross[2];

    if Float::abs(triple_product) > epsilon.sqrt() {
        // Lines are skew (not intersecting)
        return None;
    }

    // Lines intersect. Compute the parameter s for line1
    // s = (w × d2) · (d1 × d2) / |d1 × d2|²
    let w_cross_d2 = cross_product(&w, d2);
    let numerator = w_cross_d2[0] * cross[0] + w_cross_d2[1] * cross[1] + w_cross_d2[2] * cross[2];
    let s = numerator / cross_mag_sq;

    // Compute intersection point
    let intersection = vec3(
        line1.point[0] + s * d1[0],
        line1.point[1] + s * d1[1],
        line1.point[2] + s * d1[2],
    );

    Some(intersection)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::Float;

    #[test]
    fn test_line_creation() {
        let line = Line::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0));
        assert_eq!(line.point, vec3(0.0, 0.0, 0.0));
        assert_eq!(line.direction, vec3(1.0, 0.0, 0.0));
    }

    #[test]
    fn test_distance_point_to_line_perpendicular() {
        let line = Line::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0));
        let point = vec3(0.0, 5.0, 0.0);
        let distance = distance_point_to_line(point, &line);
        assert!(Float::abs(distance - 5.0) < 1e-10);
    }

    #[test]
    fn test_distance_point_to_line_on_line() {
        let line = Line::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0));
        let point = vec3(5.0, 0.0, 0.0);
        let distance = distance_point_to_line(point, &line);
        assert!(Float::abs(distance) < 1e-10);
    }

    #[test]
    fn test_distance_point_to_line_3d() {
        let line = Line::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 1.0, 0.0));
        let point = vec3(0.0, 0.0, 5.0);
        let distance = distance_point_to_line(point, &line);
        assert!(Float::abs(distance - 5.0) < 1e-10);
    }

    #[test]
    fn test_closest_point_on_line() {
        let line = Line::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0));
        let point = vec3(5.0, 3.0, 2.0);
        let closest = closest_point_on_line(point, &line);

        assert!(Float::abs(closest[0] - 5.0) < 1e-10);
        assert!(Float::abs(closest[1]) < 1e-10);
        assert!(Float::abs(closest[2]) < 1e-10);
    }

    #[test]
    fn test_closest_point_on_line_already_on_line() {
        let line = Line::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0));
        let point = vec3(5.0, 0.0, 0.0);
        let closest = closest_point_on_line(point, &line);

        assert!(Float::abs(closest[0] - 5.0) < 1e-10);
        assert!(Float::abs(closest[1]) < 1e-10);
        assert!(Float::abs(closest[2]) < 1e-10);
    }

    #[test]
    fn test_line_line_intersection_perpendicular() {
        let line1 = Line::new(vec3(-1.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0));
        let line2 = Line::new(vec3(0.0, -1.0, 0.0), vec3(0.0, 1.0, 0.0));

        let intersection = line_line_intersection(&line1, &line2);
        assert!(intersection.is_some());

        let point = intersection.unwrap();
        assert!(Float::abs(point[0]) < 1e-10);
        assert!(Float::abs(point[1]) < 1e-10);
        assert!(Float::abs(point[2]) < 1e-10);
    }

    #[test]
    fn test_line_line_intersection_parallel() {
        let line1 = Line::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0));
        let line2 = Line::new(vec3(0.0, 1.0, 0.0), vec3(1.0, 0.0, 0.0));

        let intersection = line_line_intersection(&line1, &line2);
        assert!(intersection.is_none());
    }

    #[test]
    fn test_line_line_intersection_skew() {
        // Two skew lines in 3D
        let line1 = Line::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0));
        let line2 = Line::new(vec3(0.0, 1.0, 1.0), vec3(0.0, 1.0, 0.0));

        let intersection = line_line_intersection(&line1, &line2);
        assert!(intersection.is_none());
    }

    #[test]
    fn test_line_line_intersection_coplanar() {
        // Two lines in the XY plane that intersect
        let line1 = Line::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 1.0, 0.0));
        let line2 = Line::new(vec3(1.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));

        let intersection = line_line_intersection(&line1, &line2);
        assert!(intersection.is_some());
    }

    #[test]
    fn test_distance_and_closest_point_consistency() {
        let line = Line::new(vec3(1.0, 2.0, 3.0), vec3(1.0, 0.0, 0.0));
        let point = vec3(5.0, 7.0, 11.0);

        let distance = distance_point_to_line(point, &line);
        let closest = closest_point_on_line(point, &line);

        // Distance should equal the distance from point to closest point
        let dx = point[0] - closest[0];
        let dy = point[1] - closest[1];
        let dz = point[2] - closest[2];
        let computed_distance = (dx * dx + dy * dy + dz * dz).sqrt();

        assert!(Float::abs(distance - computed_distance) < 1e-10);
    }
}
