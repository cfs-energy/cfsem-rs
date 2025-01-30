//! Calculations for 0D field sources such as dipoles.

use rayon::{
    iter::{IndexedParallelIterator, ParallelIterator},
    slice::{ParallelSlice, ParallelSliceMut},
};

use crate::{
    math::{dot3, rss3},
    MU0_OVER_4PI,
};

/// Magnetic flux density of a dipole in cartesian coordiantes.
///
/// Arguments
///
/// * loc: (m) location of the point source
/// * moment: (A-m^2) magnetic moment vector of the point source
/// * obs: (m) observation point to examine
///
/// Returns
///
/// * (bx, by, bz) [T] magnetic field components at observation point
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

/// Magnetic flux density of a dipole in cartesian coordiantes.
///
/// Arguments
///
/// * loc: (m) location of the point source
/// * moment: (A-m^2) magnetic moment vector of the point source
/// * obs: (m) observation point to examine
/// * out: (T) storage for B-field
///
/// Returns
///
/// * (bx, by, bz) [T] magnetic field components at observation point
#[inline]
pub fn flux_density_dipole(
    loc: (&[f64], &[f64], &[f64]),
    moment: (&[f64], &[f64], &[f64]),
    obs: (&[f64], &[f64], &[f64]),
    out: (&mut [f64], &mut [f64], &mut [f64]),
) -> Result<(), &'static str> {
    // Check lengths
    let m = loc.0.len();
    let n = obs.0.len();

    if loc.0.len() != m
        || loc.1.len() != m
        || loc.2.len() != m
        || moment.0.len() != m
        || moment.1.len() != m
        || moment.2.len() != m
        || obs.0.len() != n
        || obs.1.len() != n
        || obs.2.len() != n
        || out.0.len() != n
        || out.1.len() != n
        || out.2.len() != n
    {
        return Err("Input length mismatch");
    }

    // Do calcs
    for i in 0..n {
        for j in 0..m {
            let obsi = (obs.0[i], obs.1[i], obs.2[i]);
            let locj = (loc.0[j], loc.1[j], loc.2[j]);
            let momentj = (moment.0[j], moment.1[j], moment.2[j]);
            let (bx, by, bz) = flux_density_dipole_scalar(locj, momentj, obsi);
            out.0[i] += bx;
            out.1[i] += by;
            out.2[i] += bz;
        }
    }

    Ok(())
}

/// Magnetic flux density of a dipole in cartesian coordiantes.
/// Parallelized over chunks of observation points and vectorized over source points.
///
/// Arguments
///
/// * loc: (m) location of the point source
/// * moment: (A-m^2) magnetic moment vector of the point source
/// * obs: (m) observation point to examine
/// * out: (T) storage for B-field
///
/// Returns
///
/// * (bx, by, bz) [T] magnetic field components at observation point
#[inline]
pub fn flux_density_dipole_par(
    loc: (&[f64], &[f64], &[f64]),
    moment: (&[f64], &[f64], &[f64]),
    obs: (&[f64], &[f64], &[f64]),
    out: (&mut [f64], &mut [f64], &mut [f64]),
) -> Result<(), &'static str> {
    // Chunk inputs
    let ncores = std::thread::available_parallelism()
        .unwrap_or(NonZeroUsize::MIN)
        .get();

    let chunk_size = (obs.0.len() / ncores).max(1);

    let obsxc = obs.0.par_chunks(chunk_size);
    let obsyc = obs.1.par_chunks(chunk_size);
    let obszc = obs.2.par_chunks(chunk_size);

    let outxc = out.0.par_chunks_mut(chunk_size);
    let outyc = out.1.par_chunks_mut(chunk_size);
    let outzc = out.2.par_chunks_mut(chunk_size);

    outxc
        .zip(outyc.zip(outzc.zip(obsxc.zip(obsyc.zip(obszc)))))
        .try_for_each(|(outx, (outy, (outz, (obsx, (obsy, obsz)))))| {
            flux_density_dipole(loc, moment, (obsx, obsy, obsz), (outx, outy, outz))
        });

    Ok(())
}
