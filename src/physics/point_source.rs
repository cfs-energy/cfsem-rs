//! Calculations for 0D field sources such as dipoles.

use crate::{
    math::{dot3, rss3},
    MU0_OVER_4PI,
};

/// Magnetic flux density of a dipole in cartesian coordiantes.
/// 
/// Arguments
/// 
/// * loc: location of the dipole
/// * moment: magnetic moment of the dipole
/// * obs: observation point to examine
#[inline]
pub fn flux_density_dipole_scalar(
    loc: (f64, f64, f64),
    moment: (f64, f64, f64),
    obs: (f64, f64, f64),
) -> (f64, f64, f64) {
    // Radius vector decomposed into direction and magnitude
    let r = (obs.0 - loc.0, obs.1 - loc.1, obs.2 - loc.2); // [m]
    let rmag = rss3(r.0, r.1, r.2); // [m]
    let rhat = (r.0 / rmag, r.1 / rmag, r.2 / rmag); // [dimensionless]
    let rinv3 = rmag.powf(-3.0);

    // r(dot(m, r))/|r|^5 reordered to avoid computing the 5th power for improved float resolution
    let m_dot_r = dot3(moment.0, moment.1, moment.2, rhat.0, rhat.1, rhat.2);
    let rmr = (rhat.0 * m_dot_r, rhat.1 * m_dot_r, rhat.2 * m_dot_r);

    // Assemble components
    let c = 3.0 * rinv3;
    let term1 = (rmr.0 * c, rmr.1 * c, rmr.2 * c);
    let term2 = (-moment.0 * rinv3, -moment.1 * rinv3, -moment.2 * rinv3);
    let tsum = (term1.0 + term2.0, term1.1 + term2.1, term1.2 + term2.2);

    let (bx, by, bz) = (
        MU0_OVER_4PI * tsum.0,
        MU0_OVER_4PI * tsum.1,
        MU0_OVER_4PI * tsum.2,
    );

    (bx, by, bz) // [T]
}
