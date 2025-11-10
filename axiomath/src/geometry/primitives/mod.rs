//! Geometric primitives for 3D space.
//!
//! This module provides fundamental geometric primitives and operations
//! for working with 3D geometry, including lines, planes, and rays.
//!
//! # Lines
//!
//! Lines are infinite straight paths in 3D space, represented by a point
//! and a direction vector. Available operations:
//! - [`distance_point_to_line`]: Compute distance from a point to a line
//! - [`closest_point_on_line`]: Find the closest point on a line to a given point
//! - [`line_line_intersection`]: Find intersection of two lines
//!
//! # Planes
//!
//! Planes are infinite flat surfaces in 3D space, represented by a normal
//! vector and a distance from the origin. Available operations:
//! - [`distance_point_to_plane`]: Compute signed distance from a point to a plane
//!
//! # Rays
//!
//! Rays are half-lines that start at an origin and extend infinitely in
//! a direction. Available operations:
//! - [`ray_plane_intersection`]: Find where a ray intersects a plane
//!
//! # Example: Working with Lines
//!
//! ```
//! use axiomath::geometry::primitives::{Line, distance_point_to_line, closest_point_on_line};
//! use axiomath::vector::Vec3;
//!
//! // Create a line along the X-axis
//! let line = Line::new(
//!     Vec3::new(0.0, 0.0, 0.0),
//!     Vec3::new(1.0, 0.0, 0.0)
//! );
//!
//! // Find distance from a point to the line
//! let point = Vec3::new(5.0, 3.0, 4.0);
//! let distance = distance_point_to_line(point, &line);
//!
//! // Find the closest point on the line
//! let closest = closest_point_on_line(point, &line);
//! ```
//!
//! # Example: Working with Planes
//!
//! ```
//! use axiomath::geometry::primitives::{Plane, Ray, ray_plane_intersection};
//! use axiomath::vector::Vec3;
//!
//! // Create a plane at y = 5
//! let plane = Plane::new(Vec3::new(0.0, 1.0, 0.0), -5.0);
//!
//! // Create a ray from origin pointing up
//! let ray = Ray::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0));
//!
//! // Find intersection
//! let intersection = ray_plane_intersection(&ray, &plane);
//! ```
//!
//! # Example: Line-Line Intersection
//!
//! ```
//! use axiomath::geometry::primitives::{Line, line_line_intersection};
//! use axiomath::vector::Vec3;
//!
//! // Two perpendicular lines that intersect at origin
//! let line1 = Line::new(Vec3::new(-5.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
//! let line2 = Line::new(Vec3::new(0.0, -5.0, 0.0), Vec3::new(0.0, 1.0, 0.0));
//!
//! let intersection = line_line_intersection(&line1, &line2);
//! assert!(intersection.is_some());
//! ```

mod line;
mod plane;
mod triangle;

pub use line::*;
pub use plane::*;
pub use triangle::*;
