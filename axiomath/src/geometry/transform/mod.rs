//! Transformation matrices for 2D and 3D geometry.
//!
//! This module provides functions to create transformation matrices for
//! common geometric operations including translation, rotation, scaling,
//! and shearing. It uses homogeneous coordinates to represent transformations
//! as matrix multiplications.
//!
//! # 2D Transformations
//!
//! 2D transformations use 3x3 matrices with homogeneous coordinates:
//! - Translation: [`translation_matrix_2d`]
//! - Rotation: [`rotation_matrix_2d`]
//! - Scaling: [`scaling_matrix_2d`]
//! - Shearing: [`shearing_matrix_2d`]
//!
//! # 3D Transformations
//!
//! 3D transformations use 4x4 matrices with homogeneous coordinates:
//! - Translation: [`translation_matrix_3d`]
//! - Rotation around axes: [`rotation_matrix_x`], [`rotation_matrix_y`], [`rotation_matrix_z`]
//! - Rotation around arbitrary axis: [`rotation_matrix_axis_angle`]
//! - Scaling: [`scaling_matrix_3d`]
//!
//! # Applying Transformations
//!
//! Use the helper functions to apply transformations to points:
//! - 2D: [`apply_transform_2d`]
//! - 3D: [`apply_transform_3d`]
//!
//! # Example
//!
//! ```
//! use axiomath::geometry::transform::{translation_matrix_2d, rotation_matrix_2d, apply_transform_2d};
//! use axiomath::vector::Vec2;
//! use std::f64::consts::PI;
//!
//! // Create a point
//! let point = Vec2::new(1.0, 0.0);
//!
//! // Rotate 90 degrees counterclockwise
//! let rotation = rotation_matrix_2d(PI / 2.0);
//! let rotated = apply_transform_2d(point, &rotation);
//!
//! // Then translate by (5, 3)
//! let translation = translation_matrix_2d(5.0, 3.0);
//! let final_point = apply_transform_2d(rotated, &translation);
//!
//! // Or compose the transformations first
//! let composed = translation.mul_mat(&rotation);
//! let result = apply_transform_2d(point, &composed);
//! ```

mod transform_2d;
mod transform_3d;
mod homogeneous;

pub use transform_2d::*;
pub use transform_3d::*;
pub use homogeneous::*;
