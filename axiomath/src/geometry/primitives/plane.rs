//! Plane primitive and related operations.

use crate::core::traits::Scalar;
use crate::vector::{vec3, Vec3};
use num_traits::{Float, Zero};
use std::ops::{Add, Mul, Sub};

/// A plane in 3D space represented by a normal vector and a distance from origin.
///
/// # Mathematical Background
///
/// A plane in 3D can be represented by the equation:
///
/// ```text
/// n · P + d = 0
/// ```
///
/// where:
/// - `n` is the unit normal vector to the plane (perpendicular to the plane)
/// - `P` is any point [x, y, z] on the plane
/// - `d` is the signed distance from the origin to the plane
/// - · denotes the dot product
///
/// Expanding this gives the standard plane equation:
/// ```text
/// nx*x + ny*y + nz*z + d = 0
/// ```
///
/// The sign of d indicates which side of the origin the plane is on:
/// - d > 0: plane is on the opposite side of the origin from the normal direction
/// - d < 0: plane is in the direction of the normal from the origin
/// - d = 0: plane passes through the origin
///
/// # Representation
///
/// This implementation stores:
/// - `normal`: The unit normal vector (should be normalized)
/// - `distance`: The signed distance from the origin
///
/// # Alternative Representation
///
/// A plane can also be defined by a point on the plane and a normal vector.
/// The relationship is: d = -(n · point), where point is any point on the plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Plane<T>
where
    T: Scalar,
{
    /// The unit normal vector to the plane
    pub normal: Vec3<T>,
    /// The signed distance from the origin
    pub distance: T,
}

impl<T> Plane<T>
where
    T: Scalar + Float + Zero + Add<Output = T> + Mul<Output = T> + Sub<Output = T>,
{
    /// Creates a new plane from a normal vector and distance.
    ///
    /// # Parameters
    ///
    /// * `normal` - The normal vector (will be normalized)
    /// * `distance` - The signed distance from the origin
    ///
    /// # Returns
    ///
    /// A new plane with normalized normal vector
    ///
    /// # Panics
    ///
    /// Panics if the normal vector has zero length.
    ///
    /// # Example
    ///
    /// ```
    /// use axiomath::geometry::primitives::Plane;
    /// use axiomath::vector::Vec3;
    ///
    /// let plane = Plane::new(Vec3::new(0.0, 1.0, 0.0), 5.0);
    /// ```
    pub fn new(normal: Vec3<T>, distance: T) -> Self {
        // Normalize the normal vector
        let mag = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        assert!(
            mag > T::zero(),
            "Plane normal must be non-zero"
        );

        let normalized = vec3(normal[0] / mag, normal[1] / mag, normal[2] / mag);

        Plane {
            normal: normalized,
            distance,
        }
    }

    /// Creates a plane from a point on the plane and a normal vector.
    ///
    /// # Mathematical Background
    ///
    /// Given a point P on the plane and a normal vector n, we can compute
    /// the distance d using:
    ///
    /// ```text
    /// d = -(n · P)
    /// ```
    ///
    /// This ensures that the plane equation n · X + d = 0 is satisfied
    /// when X = P.
    ///
    /// # Parameters
    ///
    /// * `point` - A point on the plane
    /// * `normal` - The normal vector (will be normalized)
    ///
    /// # Returns
    ///
    /// A new plane
    ///
    /// # Panics
    ///
    /// Panics if the normal vector has zero length.
    ///
    /// # Example
    ///
    /// ```
    /// use axiomath::geometry::primitives::Plane;
    /// use axiomath::vector::Vec3;
    ///
    /// let plane = Plane::from_point_normal(
    ///     Vec3::new(0.0, 5.0, 0.0),
    ///     Vec3::new(0.0, 1.0, 0.0)
    /// );
    /// ```
    pub fn from_point_normal(point: Vec3<T>, normal: Vec3<T>) -> Self {
        // Normalize the normal vector
        let mag = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        assert!(
            mag > T::zero(),
            "Plane normal must be non-zero"
        );

        let normalized = vec3(normal[0] / mag, normal[1] / mag, normal[2] / mag);

        // Compute distance: d = -(n · point)
        let distance = -(normalized[0] * point[0] + normalized[1] * point[1] + normalized[2] * point[2]);

        Plane {
            normal: normalized,
            distance,
        }
    }
}

/// Computes the signed distance from a point to a plane.
///
/// # Mathematical Background
///
/// The signed distance from a point P to a plane with equation n · X + d = 0 is:
///
/// ```text
/// signed_distance = n · P + d
/// ```
///
/// where n is the unit normal vector.
///
/// The sign indicates which side of the plane the point is on:
/// - Positive: point is on the side the normal points to
/// - Negative: point is on the opposite side from the normal
/// - Zero: point is on the plane
///
/// The absolute value gives the perpendicular distance.
///
/// # Algorithm
///
/// 1. Compute the dot product n · P
/// 2. Add the distance d
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Parameters
///
/// * `point` - The point
/// * `plane` - The plane
///
/// # Returns
///
/// The signed distance from the point to the plane
///
/// # Example
///
/// ```
/// use axiomath::geometry::primitives::{Plane, distance_point_to_plane};
/// use axiomath::vector::Vec3;
///
/// let plane = Plane::new(Vec3::new(0.0, 1.0, 0.0), -5.0);
/// let point = Vec3::new(0.0, 10.0, 0.0);
/// let distance = distance_point_to_plane(point, &plane);
/// assert!((distance - 5.0).abs() < 1e-10);
/// ```
pub fn distance_point_to_plane<T>(point: Vec3<T>, plane: &Plane<T>) -> T
where
    T: Scalar + Float + Zero + Add<Output = T> + Mul<Output = T>,
{
    // Compute n · P + d
    plane.normal[0] * point[0] + plane.normal[1] * point[1] + plane.normal[2] * point[2] + plane.distance
}

/// A ray in 3D space, defined by an origin and a direction.
///
/// # Mathematical Background
///
/// A ray is a half-line that starts at an origin point and extends infinitely
/// in a direction. It is represented parametrically as:
///
/// ```text
/// R(t) = origin + t * direction,  t ≥ 0
/// ```
///
/// Unlike a line, t is restricted to non-negative values.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ray<T>
where
    T: Scalar,
{
    /// The origin point of the ray
    pub origin: Vec3<T>,
    /// The direction vector of the ray
    pub direction: Vec3<T>,
}

impl<T> Ray<T>
where
    T: Scalar,
{
    /// Creates a new ray from an origin and direction.
    ///
    /// # Parameters
    ///
    /// * `origin` - The origin point
    /// * `direction` - The direction vector
    ///
    /// # Example
    ///
    /// ```
    /// use axiomath::geometry::primitives::Ray;
    /// use axiomath::vector::Vec3;
    ///
    /// let ray = Ray::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
    /// ```
    pub fn new(origin: Vec3<T>, direction: Vec3<T>) -> Self {
        Ray { origin, direction }
    }
}

/// Computes the intersection point of a ray and a plane.
///
/// # Mathematical Background
///
/// A ray R(t) = O + t*d intersects a plane n · X + D = 0 when:
///
/// ```text
/// n · (O + t*d) + D = 0
/// n · O + t(n · d) + D = 0
/// t = -(n · O + D) / (n · d)
/// ```
///
/// The ray intersects the plane if:
/// 1. n · d ≠ 0 (ray is not parallel to plane)
/// 2. t ≥ 0 (intersection is in the positive direction of the ray)
///
/// If n · d = 0, the ray is parallel to the plane and either:
/// - Never intersects (if ray origin is not on plane)
/// - Is entirely contained in the plane (if ray origin is on plane)
///
/// # Algorithm
///
/// 1. Compute the denominator: n · d
/// 2. Check if ray is parallel to plane (denominator ≈ 0)
/// 3. Compute the parameter t
/// 4. Check if t ≥ 0 (intersection is in positive ray direction)
/// 5. Compute the intersection point: O + t*d
///
/// # Complexity
///
/// - Time: O(1)
/// - Space: O(1)
///
/// # Parameters
///
/// * `ray` - The ray
/// * `plane` - The plane
///
/// # Returns
///
/// `Some(point)` if the ray intersects the plane, `None` if parallel or behind ray origin
///
/// # Example
///
/// ```
/// use axiomath::geometry::primitives::{Ray, Plane, ray_plane_intersection};
/// use axiomath::vector::Vec3;
///
/// let ray = Ray::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0));
/// let plane = Plane::new(Vec3::new(0.0, 1.0, 0.0), -5.0);
///
/// let intersection = ray_plane_intersection(&ray, &plane);
/// assert!(intersection.is_some());
///
/// let point = intersection.unwrap();
/// assert!((point.y - 5.0).abs() < 1e-10);
/// ```
pub fn ray_plane_intersection<T>(ray: &Ray<T>, plane: &Plane<T>) -> Option<Vec3<T>>
where
    T: Scalar + Float + Zero + Add<Output = T> + Mul<Output = T> + Sub<Output = T>,
{
    // Compute n · d (denominator)
    let denom = plane.normal[0] * ray.direction[0]
        + plane.normal[1] * ray.direction[1]
        + plane.normal[2] * ray.direction[2];

    // Threshold for considering ray parallel to plane
    let epsilon = T::epsilon() * (T::one() + T::one() + T::one() + T::one()); // 4ε

    // Check if ray is parallel to plane
    if Float::abs(denom) < epsilon {
        return None;
    }

    // Compute n · O + D (numerator)
    let numer = plane.normal[0] * ray.origin[0]
        + plane.normal[1] * ray.origin[1]
        + plane.normal[2] * ray.origin[2]
        + plane.distance;

    // Compute t
    let t = -numer / denom;

    // Check if intersection is in positive direction (t ≥ 0)
    if t < T::zero() {
        return None;
    }

    // Compute intersection point
    let intersection = vec3(
        ray.origin[0] + t * ray.direction[0],
        ray.origin[1] + t * ray.direction[1],
        ray.origin[2] + t * ray.direction[2],
    );

    Some(intersection)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::Float;

    #[test]
    fn test_plane_creation() {
        let plane = Plane::new(vec3(0.0, 1.0, 0.0), 5.0);
        assert!(Float::abs(plane.normal[1] - 1.0) < 1e-10);
        assert_eq!(plane.distance, 5.0);
    }

    #[test]
    fn test_plane_from_point_normal() {
        let point = vec3(0.0, 5.0, 0.0);
        let normal = vec3(0.0, 1.0, 0.0);
        let plane = Plane::from_point_normal(point, normal);

        // The plane equation should be satisfied by the point
        let dist = distance_point_to_plane(point, &plane);
        assert!(Float::abs(dist) < 1e-10);
    }

    #[test]
    fn test_plane_normalizes_normal() {
        let plane = Plane::new(vec3(0.0, 2.0, 0.0), 5.0);
        let mag = (plane.normal[0] * plane.normal[0]
            + plane.normal[1] * plane.normal[1]
            + plane.normal[2] * plane.normal[2])
            .sqrt();
        assert!(Float::abs(mag - 1.0) < 1e-10);
    }

    #[test]
    fn test_distance_point_to_plane_on_plane() {
        let plane = Plane::new(vec3(0.0, 1.0, 0.0), -5.0);
        let point = vec3(0.0, 5.0, 0.0);
        let distance = distance_point_to_plane(point, &plane);
        assert!(Float::abs(distance) < 1e-10);
    }

    #[test]
    fn test_distance_point_to_plane_above() {
        let plane = Plane::new(vec3(0.0, 1.0, 0.0), -5.0);
        let point = vec3(0.0, 10.0, 0.0);
        let distance = distance_point_to_plane(point, &plane);
        assert!(Float::abs(distance - 5.0) < 1e-10);
    }

    #[test]
    fn test_distance_point_to_plane_below() {
        let plane = Plane::new(vec3(0.0, 1.0, 0.0), -5.0);
        let point = vec3(0.0, 0.0, 0.0);
        let distance = distance_point_to_plane(point, &plane);
        assert!(Float::abs(distance + 5.0) < 1e-10);
    }

    #[test]
    fn test_distance_point_to_plane_origin() {
        let plane = Plane::new(vec3(0.0, 1.0, 0.0), 0.0);
        let point = vec3(0.0, 0.0, 0.0);
        let distance = distance_point_to_plane(point, &plane);
        assert!(Float::abs(distance) < 1e-10);
    }

    #[test]
    fn test_ray_creation() {
        let ray = Ray::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0));
        assert_eq!(ray.origin, vec3(0.0, 0.0, 0.0));
        assert_eq!(ray.direction, vec3(1.0, 0.0, 0.0));
    }

    #[test]
    fn test_ray_plane_intersection_perpendicular() {
        let ray = Ray::new(vec3(0.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));
        let plane = Plane::new(vec3(0.0, 1.0, 0.0), -5.0);

        let intersection = ray_plane_intersection(&ray, &plane);
        assert!(intersection.is_some());

        let point = intersection.unwrap();
        assert!(Float::abs(point[0]) < 1e-10);
        assert!(Float::abs(point[1] - 5.0) < 1e-10);
        assert!(Float::abs(point[2]) < 1e-10);
    }

    #[test]
    fn test_ray_plane_intersection_angled() {
        let ray = Ray::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 1.0, 0.0));
        let plane = Plane::new(vec3(0.0, 1.0, 0.0), -5.0);

        let intersection = ray_plane_intersection(&ray, &plane);
        assert!(intersection.is_some());

        let point = intersection.unwrap();
        assert!(Float::abs(point[1] - 5.0) < 1e-10);
        assert!(Float::abs(point[0] - 5.0) < 1e-10);
    }

    #[test]
    fn test_ray_plane_intersection_parallel() {
        let ray = Ray::new(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0));
        let plane = Plane::new(vec3(0.0, 1.0, 0.0), -5.0);

        let intersection = ray_plane_intersection(&ray, &plane);
        assert!(intersection.is_none());
    }

    #[test]
    fn test_ray_plane_intersection_behind_ray() {
        // Ray pointing away from plane
        let ray = Ray::new(vec3(0.0, 0.0, 0.0), vec3(0.0, -1.0, 0.0));
        let plane = Plane::new(vec3(0.0, 1.0, 0.0), -5.0);

        let intersection = ray_plane_intersection(&ray, &plane);
        assert!(intersection.is_none());
    }

    #[test]
    fn test_ray_plane_intersection_origin_on_plane() {
        // Ray starts on plane and points into it (not parallel)
        let ray = Ray::new(vec3(0.0, 5.0, 0.0), vec3(1.0, 1.0, 0.0));
        let plane = Plane::new(vec3(0.0, 1.0, 0.0), -5.0);

        let intersection = ray_plane_intersection(&ray, &plane);
        // Ray origin is on the plane, so t should be 0 or very small
        assert!(intersection.is_some());

        let point = intersection.unwrap();
        // Point should be very close to the ray origin
        assert!(Float::abs(point[0] - ray.origin[0]) < 1e-8);
        assert!(Float::abs(point[1] - ray.origin[1]) < 1e-8);
    }

    #[test]
    fn test_plane_consistency() {
        // Test that from_point_normal is consistent with distance calculation
        let point_on_plane = vec3(3.0, 4.0, 5.0);
        let normal = vec3(1.0, 1.0, 1.0);
        let plane = Plane::from_point_normal(point_on_plane, normal);

        // Any point we declared to be on the plane should have distance 0
        let dist = distance_point_to_plane(point_on_plane, &plane);
        assert!(Float::abs(dist) < 1e-10);
    }
}
