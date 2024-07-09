#![allow(non_snake_case)]

pub mod math;
pub mod mesh;
pub mod physics;

/// (H/m) vacuum magnetic permeability.
/// Value from 2022 CODATA recommended values, [NIST SPI 961](https://physics.nist.gov/cuu/pdf/wall_2022.pdf).
pub const MU_0: f64 = 0.999_999_999_87 * core::f64::consts::PI * 4e-7; // [H/m]

/// (H/m) Recurring constant multiple of `mu_0`
pub const MU0_OVER_4PI: f64 = MU_0 / (4.0 * core::f64::consts::PI);
