//! Quaternions for representing 3D rotations.
//!
//! This module provides a complete quaternion implementation for 3D rotation
//! representation and manipulation. Quaternions avoid gimbal lock and provide
//! smooth interpolation capabilities that Euler angles cannot match.
//!
//! # Overview
//!
//! A quaternion is a four-component number q = w + xi + yj + zk that can represent
//! rotations in 3D space. Unit quaternions (|q| = 1) form a double cover of the
//! 3D rotation group SO(3).
//!
//! # Main Types and Functions
//!
//! - [`Quaternion`]: The quaternion type with w, x, y, z components
//! - [`multiply`]: Hamilton product (quaternion multiplication)
//! - [`conjugate`]: Quaternion conjugate
//! - [`inverse`]: Quaternion inverse
//! - [`magnitude`]: Quaternion magnitude/norm
//! - [`normalize`]: Normalize to unit quaternion
//! - [`from_axis_angle`]: Create quaternion from axis-angle representation
//! - [`to_rotation_matrix_3x3`]: Convert to 3x3 rotation matrix
//! - [`to_rotation_matrix_4x4`]: Convert to 4x4 rotation matrix
//! - [`slerp`]: Spherical linear interpolation between quaternions
//!
//! # Example: Basic Quaternion Operations
//!
//! ```
//! use axiomath::geometry::quaternion::{Quaternion, multiply, magnitude, normalize};
//!
//! let q1 = Quaternion::new(1.0, 0.0, 0.0, 0.0);
//! let q2 = Quaternion::new(0.0, 1.0, 0.0, 0.0);
//!
//! let product = multiply(q1, q2);
//! let mag = magnitude(product);
//! let normalized = normalize(product);
//! ```
//!
//! # Example: Creating Rotations
//!
//! ```
//! use axiomath::geometry::quaternion::from_axis_angle;
//! use axiomath::vector::Vec3;
//! use std::f64::consts::PI;
//!
//! // 90 degree rotation around Z axis
//! let axis = Vec3::new(0.0, 0.0, 1.0);
//! let rotation = from_axis_angle(axis, PI / 2.0);
//! ```
//!
//! # Example: Interpolating Rotations
//!
//! ```
//! use axiomath::geometry::quaternion::{from_axis_angle, slerp};
//! use axiomath::vector::Vec3;
//! use std::f64::consts::PI;
//!
//! let axis = Vec3::new(0.0, 1.0, 0.0);
//! let start = from_axis_angle(axis, 0.0);
//! let end = from_axis_angle(axis, PI);
//!
//! // Smooth interpolation at 50%
//! let halfway = slerp(start, end, 0.5);
//! ```
//!
//! # Example: Converting to Matrix
//!
//! ```
//! use axiomath::geometry::quaternion::{from_axis_angle, to_rotation_matrix_4x4};
//! use axiomath::geometry::transform::apply_transform_3d;
//! use axiomath::vector::Vec3;
//! use std::f64::consts::PI;
//!
//! let axis = Vec3::new(0.0, 0.0, 1.0);
//! let q = from_axis_angle(axis, PI / 2.0);
//! let matrix = to_rotation_matrix_4x4(q);
//!
//! // Apply rotation to a point
//! let point = Vec3::new(1.0, 0.0, 0.0);
//! let rotated = apply_transform_3d(point, &matrix);
//! ```

mod operations;
mod slerp;

pub use operations::*;
pub use slerp::*;
