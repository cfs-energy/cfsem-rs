//! Biot-Savart calculations for B-field from current filaments.
use std::num::NonZeroUsize;

use rayon::{
    iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator},
    slice::{ParallelSlice, ParallelSliceMut},
};

use crate::math::{cross3, dot3};

/// Biot-Savart calculation for B-field contribution from many current filament
/// segments to many observation points. This variant of the function is
/// parallelized over chunks of observation points.
///
/// # Arguments
///
/// * `xyzp`:     (m) Observation point coords, each length `n`
/// * `xyzfil`:   (m) Filament origin coords (start of segment), each length `m`
/// * `dlxyzfil`: (m) Filament segment length deltas, each length `m`
/// * `ifil`:     (A) Filament current, length `m`
/// * `out`:      (T) bx, by, bz at observation points, each length `n`
pub fn flux_density_biot_savart_par(
    xyzp: (&[f64], &[f64], &[f64]),
    xyzfil: (&[f64], &[f64], &[f64]),
    dlxyzfil: (&[f64], &[f64], &[f64]),
    ifil: &[f64],
    out: (&mut [f64], &mut [f64], &mut [f64]),
) -> Result<(), &'static str> {
    // Chunk inputs
    let ncores = std::thread::available_parallelism()
        .unwrap_or(NonZeroUsize::MIN)
        .get();

    let n = (xyzp.0.len() / ncores).max(1);

    fn c(v: &[f64], n: usize) -> Vec<&[f64]> {
        v.par_chunks(n).collect::<Vec<&[f64]>>()
    }

    //    Some very small allocations that could be eliminated by using .chunk()
    //    but would produce an even more cumbersome zip chain
    let xpc = c(xyzp.0, n);
    let ypc = c(xyzp.1, n);
    let zpc = c(xyzp.2, n);

    let bxc = out.0.par_chunks_mut(n);
    let byc = out.1.par_chunks_mut(n);
    let bzc = out.2.par_chunks_mut(n);

    let nchunks = xpc.len();

    if ypc.len() != nchunks
        || zpc.len() != nchunks
        || bxc.len() != nchunks
        || byc.len() != nchunks
        || bzc.len() != nchunks
    {
        return Err("Dimension mismatch");
    }

    // Run calcs
    (0..nchunks)
        .into_par_iter()
        .zip(bxc.zip(byc.zip(bzc)))
        .try_for_each(|(i, (bx, (by, bz)))| {
            flux_density_biot_savart(
                (xpc[i], ypc[i], zpc[i]),
                xyzfil,
                dlxyzfil,
                ifil,
                (bx, by, bz),
            )
        })?;

    Ok(())
}

/// Biot-Savart calculation for B-field contribution from many current filament
/// segments to many observation points.
///
/// # Arguments
///
/// * `xyzp`:     (m) Observation point coords, each length `n`
/// * `xyzfil`:   (m) Filament origin coords (start of segment), each length `m`
/// * `dlxyzfil`: (m) Filament segment length deltas, each length `m`
/// * `ifil`:     (A) Filament current, length `m`
/// * `out`:      (T) bx, by, bz at observation points, each length `n`
pub fn flux_density_biot_savart(
    xyzp: (&[f64], &[f64], &[f64]),
    xyzfil: (&[f64], &[f64], &[f64]),
    dlxyzfil: (&[f64], &[f64], &[f64]),
    ifil: &[f64],
    out: (&mut [f64], &mut [f64], &mut [f64]),
) -> Result<(), &'static str> {
    // Unpack
    let (xp, yp, zp) = xyzp;
    let (xfil, yfil, zfil) = xyzfil;
    let (dlxfil, dlyfil, dlzfil) = dlxyzfil;

    let (bx, by, bz) = out;

    // Check lengths; if there is any possibility of mismatch,
    // the compiler will bypass vectorization
    let n = xfil.len();
    let m = xp.len();

    if xp.len() != m
        || yp.len() != m
        || zp.len() != m
        || xfil.len() != n
        || yfil.len() != n
        || zfil.len() != n
        || dlxfil.len() != n
        || dlyfil.len() != n
        || dlzfil.len() != n
        || ifil.len() != n
    {
        return Err("Input length mismatch");
    }

    let mu0_over_4pi: f64 = 1e-7; // [H/m] Collapse some algebra to reduce float error

    // For each filament, evaluate the contribution to each observation point
    for i in 0..n {
        // Get filament midpoint
        let dlxi = dlxfil[i]; // [m]
        let dlyi = dlyfil[i]; // [m]
        let dlzi = dlzfil[i]; // [m]
        let xmid = dlxi.mul_add(0.5, xfil[i]); // [m]
        let ymid = dlyi.mul_add(0.5, yfil[i]); // [m]
        let zmid = dlzi.mul_add(0.5, zfil[i]); // [m]

        // Get filament current and bake in the constant factor
        let ifil_scaled = mu0_over_4pi * ifil[i];

        for j in 0..m {
            // Get distance from middle of the filament segment to the observation point
            let rx = xp[j] - xmid; // [m]
            let ry = yp[j] - ymid; // [m]
            let rz = zp[j] - zmid; // [m]

            // Do 1/r^3 operation with an ordering that improves float error by eliminating
            // the actual cube operation and using fused multiply-add to reduce roundoff events,
            // then rolling the result into the factor that is constant between all contributions.
            let sumsq = dot3(rx, ry, rz, rx, ry, rz);
            let rnorm3_inv = sumsq.powf(-1.5); // [m^-3]

            // This factor is constant across all x, y, and z components
            let c = rnorm3_inv * ifil_scaled;

            // Evaluate the cross products for each axis component
            // separately using mul_add which would not be assumed usable
            // in a more general implementation.
            let (cx, cy, cz) = cross3(dlxi, dlyi, dlzi, rx, ry, rz);

            // Sum up the contributions at each observation point on each axis
            // using fused multiply-add again to reduce roundoff error and slightly improve speed
            bx[j] = c.mul_add(cx, bx[j]);
            by[j] = c.mul_add(cy, by[j]);
            bz[j] = c.mul_add(cz, bz[j]);
        }
    }
    Ok(())
}
