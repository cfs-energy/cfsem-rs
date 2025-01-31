#![allow(non_snake_case)]

pub mod math;
pub mod mesh;
pub mod physics;

#[cfg(test)]
pub(crate) mod testing;

/// (H/m) vacuum magnetic permeability.
/// Value from 2022 CODATA recommended values, [NIST SPI 961](https://physics.nist.gov/cuu/pdf/wall_2022.pdf).
pub const MU_0: f64 = 0.999_999_999_87 * core::f64::consts::PI * 4e-7; // [H/m]

/// (H/m) Recurring constant multiple of `mu_0`
pub const MU0_OVER_4PI: f64 = MU_0 / (4.0 * core::f64::consts::PI);

#[macro_use]
pub(crate) mod macros {
    macro_rules! check_length {
        ($n:expr, $($y:expr),+) => {
            $(  // Repeat for all y
                if $y.len() != $n {
                    return Err("Length mismatch");
                }
            )+
        };
    }

    macro_rules! check_length_3tup {
        ($n:expr, $x:expr) => {
            if $x.0.len() != $n || x.1.len() != $n || x.2.len() != $n {
                return Err("Length mismatch");
            }
        };
    }

    // Publish macros within crate
    pub(crate) use check_length;
    pub(crate) use check_length_3tup;
}
