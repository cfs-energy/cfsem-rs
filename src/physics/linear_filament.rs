//! Magnetics calculations for piecewise-linear current filaments.
use std::num::NonZeroUsize;

use rayon::{
    iter::{IndexedParallelIterator, ParallelIterator},
    slice::{ParallelSlice, ParallelSliceMut},
};

use crate::math::{cross3, dot3, rss3};
use crate::MU0_OVER_4PI;

/// Estimate the mutual inductance between two piecewise-linear current filaments.
///
/// # Arguments
///
/// * `xyzfil0`:         (m) filament origin coordinates for first path, length `n`
/// * `dlxyzfil0`:       (m) filament segment lengths for first path, length `n`
/// * `xyzfil1`:         (m) filament origin coordinates for second path, length `n`
/// * `dlxyzfil1`:       (m) filament segment lengths for second path, length `n`
/// * `self_inductance`: Flag for whether this calc is being used for self-inductance,
///                      in which case segment self-field terms are replaced with a hand-calc
///
/// # Commentary
///
/// Uses Neumann's Formula for the mutual inductance of arbitrary loops, which is
/// originally from \[2\] and can be found in a more friendly format on wikipedia.
///
/// When `self_inductance` flag is set, zeroes-out the contributions from self-pairings
/// to resolve the thin-filament self-inductance singularity and replaces the
/// segment self-inductance term with an analytic value from \[3\].
///
/// # Assumptions
///
/// * Thin, well-behaved filaments
/// * Uniform current distribution within segments
///     * Low frequency operation; no skin effect
///       (which would reduce the segment self-field term)
/// * Vacuum permeability everywhere
/// * Each filament has a constant current in all segments
///   (otherwise we need an inductance matrix)
///
/// # References
///
///   \[1\] “Inductance,” Wikipedia. Dec. 12, 2022. Accessed: Jan. 23, 2023. \[Online\].
///         Available: <https://en.wikipedia.org/w/index.php?title=Inductance>
///
///   \[2\] F. E. Neumann, “Allgemeine Gesetze der inducirten elektrischen Ströme,”
///         Jan. 1846, doi: [10.1002/andp.18461430103](https://doi.org/10.1002/andp.18461430103).
///
///   \[3\] R. Dengler, “Self inductance of a wire loop as a curve integral,”
///         AEM, vol. 5, no. 1, p. 1, Jan. 2016, doi: [10.7716/aem.v5i1.331](https://doi.org/10.7716/aem.v5i1.331).
pub fn inductance_piecewise_linear_filaments(
    xyzfil0: (&[f64], &[f64], &[f64]),
    dlxyzfil0: (&[f64], &[f64], &[f64]),
    xyzfil1: (&[f64], &[f64], &[f64]),
    dlxyzfil1: (&[f64], &[f64], &[f64]),
    self_inductance: bool,
) -> Result<f64, &'static str> {
    // Unpack
    let (xfil0, yfil0, zfil0) = xyzfil0;
    let (dlxfil0, dlyfil0, dlzfil0) = dlxyzfil0;
    let (xfil1, yfil1, zfil1) = xyzfil1;
    let (dlxfil1, dlyfil1, dlzfil1) = dlxyzfil1;

    // Check lengths; Error if they do not match
    let n = xfil0.len();
    if xfil0.len() != n
        || yfil0.len() != n
        || zfil0.len() != n
        || dlxfil0.len() != n
        || dlyfil0.len() != n
        || dlzfil0.len() != n
    {
        return Err("Input length mismatch");
    }

    let m = xfil1.len();
    if xfil1.len() != m
        || yfil1.len() != m
        || zfil1.len() != m
        || dlxfil1.len() != m
        || dlyfil1.len() != m
        || dlzfil1.len() != m
    {
        return Err("Input length mismatch");
    }

    if self_inductance {
        if m != n {
            return Err("For self-inductance runs, the two paths must be the same length and should be identical");
        }
    }

    let mut inductance: f64 = 0.0; // [H], although it is in [m] until the final calc
    let mut total_length: f64 = 0.0; // [m]
    for i in 0..n {
        // Filament i midpoint
        let dlxi = dlxfil0[i]; // [m]
        let dlyi = dlyfil0[i]; // [m]
        let dlzi = dlzfil0[i]; // [m]
        let xmidi = dlxi.mul_add(0.5, xfil0[i]); // [m]
        let ymidi = dlyi.mul_add(0.5, yfil0[i]); // [m]
        let zmidi = dlzi.mul_add(0.5, zfil0[i]); // [m]

        // Accumulate total length if we need it
        if self_inductance {
            total_length += rss3(dlxi, dlyi, dlzi);
        }

        for j in 0..m {
            // Skip self-interaction terms which are handled separately
            if self_inductance && i == j {
                continue;
            }

            // Filament j midpoint
            let dlxj = dlxfil1[j]; // [m]
            let dlyj = dlyfil1[j]; // [m]
            let dlzj = dlzfil1[j]; // [m]
            let xmidj = dlxj.mul_add(0.5, xfil1[j]); // [m]
            let ymidj = dlyj.mul_add(0.5, yfil1[j]); // [m]
            let zmidj = dlzj.mul_add(0.5, zfil1[j]); // [m]

            // Distance between midpoints
            let rx = xmidi - xmidj;
            let ry = ymidi - ymidj;
            let rz = zmidi - zmidj;
            let dist = rss3(rx, ry, rz);

            // Dot product of segment vectors
            let dxdot = dot3(dlxi, dlyi, dlzi, dlxj, dlyj, dlzj);

            inductance += dxdot / dist;
        }
    }

    // Add self-inductance of individual filament segments
    // if this is a self-inductance calc
    if self_inductance {
        inductance += 0.5 * total_length;
    }

    // Finally, do the shared constant factor
    inductance *= MU0_OVER_4PI;

    Ok(inductance)
}

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
pub fn flux_density_linear_filament_par(
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

    let xpc = xyzp.0.par_chunks(n);
    let ypc = xyzp.1.par_chunks(n);
    let zpc = xyzp.2.par_chunks(n);

    let bxc = out.0.par_chunks_mut(n);
    let byc = out.1.par_chunks_mut(n);
    let bzc = out.2.par_chunks_mut(n);

    // Run calcs
    bxc.zip(byc.zip(bzc.zip(xpc.zip(ypc.zip(zpc)))))
        .try_for_each(|(bx, (by, (bz, (xp, (yp, zp)))))| {
            flux_density_linear_filament((xp, yp, zp), xyzfil, dlxyzfil, ifil, (bx, by, bz))
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
pub fn flux_density_linear_filament(
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

    for j in 0..m {
        bx[j] = 0.0;
        by[j] = 0.0;
        bz[j] = 0.0;
    }

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
        let current = ifil[i];

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
            let c = current * rnorm3_inv;

            // Evaluate the cross products for each axis component
            // separately using mul_add which would not be assumed usable
            // in a more general implementation.
            let (cx, cy, cz) = cross3(dlxi, dlyi, dlzi, rx, ry, rz);

            // Sum up the contributions at each observation point on each axis
            // using fused multiply-add again to reduce roundoff error and slightly improve speed
            // The B-field contributions added here *do not* have units of tesla,
            // but the sum will have units of tesla after it is multiplied by mu_0 / (4 * pi) below.
            bx[j] = c.mul_add(cx, bx[j]);
            by[j] = c.mul_add(cy, by[j]);
            bz[j] = c.mul_add(cz, bz[j]);
        }
    }

    for j in 0..m {
        bx[j] *= MU0_OVER_4PI;
        by[j] *= MU0_OVER_4PI;
        bz[j] *= MU0_OVER_4PI;
    }

    Ok(())
}

/// Vector potential calculation for A-field contribution from many current filament
/// segments to many observation points. This variant of the function is
/// parallelized over chunks of observation points.
///
/// # Arguments
///
/// * `xyzp`:     (m) Observation point coords, each length `n`
/// * `xyzfil`:   (m) Filament origin coords (start of segment), each length `m`
/// * `dlxyzfil`: (m) Filament segment length deltas, each length `m`
/// * `ifil`:     (A) Filament current, length `m`
/// * `out`:      (V-s/m) ax, ay, az at observation points, each length `n`
pub fn vector_potential_linear_filament_par(
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

    let xpc = xyzp.0.par_chunks(n);
    let ypc = xyzp.1.par_chunks(n);
    let zpc = xyzp.2.par_chunks(n);

    let bxc = out.0.par_chunks_mut(n);
    let byc = out.1.par_chunks_mut(n);
    let bzc = out.2.par_chunks_mut(n);

    // Run calcs
    bxc.zip(byc.zip(bzc.zip(xpc.zip(ypc.zip(zpc)))))
        .try_for_each(|(bx, (by, (bz, (xp, (yp, zp)))))| {
            vector_potential_linear_filament((xp, yp, zp), xyzfil, dlxyzfil, ifil, (bx, by, bz))
        })?;

    Ok(())
}

/// Vector potential calculation for A-field contribution from many current filament
/// segments to many observation points.
///
/// # Arguments
///
/// * `xyzp`:     (m) Observation point coords, each length `n`
/// * `xyzfil`:   (m) Filament origin coords (start of segment), each length `m`
/// * `dlxyzfil`: (m) Filament segment length deltas, each length `m`
/// * `ifil`:     (A) Filament current, length `m`
/// * `out`:      (V-s/m) ax, ay, az at observation points, each length `n`
pub fn vector_potential_linear_filament(
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

    let (ax, ay, az) = out;

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

    for j in 0..m {
        ax[j] = 0.0;
        ay[j] = 0.0;
        az[j] = 0.0;
    }

    // For each filament, evaluate the contribution to each observation point
    for i in 0..n {
        // Get filament midpoint
        let dlxi = dlxfil[i]; // [m]
        let dlyi = dlyfil[i]; // [m]
        let dlzi = dlzfil[i]; // [m]
        let xmid = dlxi.mul_add(0.5, xfil[i]); // [m]
        let ymid = dlyi.mul_add(0.5, yfil[i]); // [m]
        let zmid = dlzi.mul_add(0.5, zfil[i]); // [m]
        let current = ifil[i];

        for j in 0..m {
            // Get distance from middle of the filament segment to the observation point
            let rx = xp[j] - xmid; // [m]
            let ry = yp[j] - ymid; // [m]
            let rz = zp[j] - zmid; // [m]
            let rnorm = rss3(rx, ry, rz); // [m]

            let c = current / rnorm;

            // The contributions added here have units of amperes,
            // but the sum will have units of volt second per meter after being multiplied
            // by mu_0 / (4 * pi) below.
            ax[j] = c.mul_add(dlxi, ax[j]);
            ay[j] = c.mul_add(dlyi, ay[j]);
            az[j] = c.mul_add(dlzi, az[j])
        }
    }

    for j in 0..m {
        ax[j] *= MU0_OVER_4PI;
        ay[j] *= MU0_OVER_4PI;
        az[j] *= MU0_OVER_4PI;
    }

    Ok(())
}

#[cfg(test)]
mod test {
    use std::f64::consts::PI;

    use super::*;

    /// Div/0-resistant approximate comparison
    fn approx(truth: f64, val: f64, rtol: f64, atol: f64) -> bool {
        let abs_err = (val - truth).abs();
        let lim = rtol * truth.abs() + atol;
        abs_err < lim
    }

    /// Check that B = curl(A)
    #[test]
    fn test_vector_potential() {
        // One super basic filament as the source
        let xyz = [0.0];
        let dlxyz = [1.0];

        // Build a second scattering of filament locations as the target
        const NFIL: usize = 10;
        let xfil2: Vec<f64> = (0..NFIL).map(|i| (i as f64).sin() + PI).collect();
        let yfil2: Vec<f64> = (0..NFIL).map(|i| (i as f64).cos() - PI).collect();
        let zfil2: Vec<f64> = (0..NFIL)
            .map(|i| (i as f64) - (NFIL as f64) / 2.0 + PI)
            .collect();
        let xyzfil2 = (
            &xfil2[..=NFIL - 2],
            &yfil2[..=NFIL - 2],
            &zfil2[..=NFIL - 2],
        );

        let dlxfil2: Vec<f64> = (0..=NFIL - 2).map(|i| xfil2[i + 1] - xfil2[i]).collect();
        let dlyfil2: Vec<f64> = (0..=NFIL - 2).map(|i| yfil2[i + 1] - yfil2[i]).collect();
        let dlzfil2: Vec<f64> = (0..=NFIL - 2).map(|i| zfil2[i + 1] - zfil2[i]).collect();
        let dlxyzfil2 = (&dlxfil2[..], &dlyfil2[..], &dlzfil2[..]);

        let xmid2: Vec<f64> = xfil2
            .iter()
            .zip(dlxfil2.iter())
            .map(|(x, dx)| x + dx / 2.0)
            .collect();
        let ymid2: Vec<f64> = yfil2
            .iter()
            .zip(dlyfil2.iter())
            .map(|(x, dx)| x + dx / 2.0)
            .collect();
        let zmid2: Vec<f64> = zfil2
            .iter()
            .zip(dlzfil2.iter())
            .map(|(x, dx)| x + dx / 2.0)
            .collect();

        // Check against Neumann's formula for mutual inductance
        let outx = &mut [0.0; NFIL - 1];
        let outy = &mut [0.0; NFIL - 1];
        let outz = &mut [0.0; NFIL - 1];
        vector_potential_linear_filament(
            (&xmid2, &ymid2, &zmid2),
            (&xyz, &xyz, &xyz),
            (&dlxyz, &dlxyz, &dlxyz),
            &[1.0],
            (outx, outy, outz),
        )
        .unwrap();
        // Here the mutual inductance of the two filaments is calculated from the
        // vector potential at filament 2 due to 1 ampere of current flowing in filament 1.
        // By Stokes' theorem, the line integral of A over filament 2 is equal to the
        // magnetic flux through a surface bounded by filament 2. The flux through
        // filament 2 due to 1 ampere of current in filament 1 is the mutual inductance.
        // (We are stretching the applicability of Stokes' therorem because the filaments
        // are not closed loops)
        let a_dot_dl: Vec<f64> = (0..NFIL - 1)
            .map(|i| outx[i] * dlxfil2[i] + outy[i] * dlyfil2[i] + outz[i] * dlzfil2[i])
            .collect();
        let m_from_a = a_dot_dl.iter().sum();
        let m = inductance_piecewise_linear_filaments(
            (&xyz, &xyz, &xyz),
            (&dlxyz, &dlxyz, &dlxyz),
            xyzfil2,
            dlxyzfil2,
            false,
        )
        .unwrap();
        assert!(approx(m, m_from_a, 1e-10, 1e-15));

        let vp = |x: f64, y: f64, z: f64| {
            let mut outx = [0.0];
            let mut outy = [0.0];
            let mut outz = [0.0];

            vector_potential_linear_filament(
                (&[x], &[y], &[z]),
                (&xyz, &xyz, &xyz),
                (&dlxyz, &dlxyz, &dlxyz),
                &[1.0],
                (&mut outx, &mut outy, &mut outz),
            )
            .unwrap();

            (outx[0], outy[0], outz[0])
        };

        let vals = [
            0.25, 0.5, 2.5, 10.0, 100.0, 1000.0, -1000.0, -100.0, -10.0, -2.5, -0.5, -0.25,
        ];
        // finite diff delta needs to be small enough to be accurate
        // but large enough that we can tell the difference between adjacent points
        // that are very far from the origin
        let eps = 1e-7;
        for x in vals.iter() {
            for y in vals.iter() {
                for z in vals.iter() {
                    let x = &(x + 1e-2); // Slightly adjust to avoid nans
                    let y = &(y + 1e-2);
                    let z = &(z - 1e-2);

                    // Brute-force jac because we're only using it once
                    let mut da = [[0.0; 3]; 3];
                    // da/dx
                    let (ax0, ay0, az0) = vp(*x - eps, *y, *z);
                    let (ax1, ay1, az1) = vp(*x + eps, *y, *z);
                    da[0][0] = (ax1 - ax0) / (2.0 * eps);
                    da[0][1] = (ay1 - ay0) / (2.0 * eps);
                    da[0][2] = (az1 - az0) / (2.0 * eps);

                    // da/dy
                    let (ax0, ay0, az0) = vp(*x, *y - eps, *z);
                    let (ax1, ay1, az1) = vp(*x, *y + eps, *z);
                    da[1][0] = (ax1 - ax0) / (2.0 * eps);
                    da[1][1] = (ay1 - ay0) / (2.0 * eps);
                    da[1][2] = (az1 - az0) / (2.0 * eps);

                    // da/dz
                    let (ax0, ay0, az0) = vp(*x, *y, *z - eps);
                    let (ax1, ay1, az1) = vp(*x, *y, *z + eps);
                    da[2][0] = (ax1 - ax0) / (2.0 * eps);
                    da[2][1] = (ay1 - ay0) / (2.0 * eps);
                    da[2][2] = (az1 - az0) / (2.0 * eps);

                    // B = curl(A)
                    let daz_dy = da[1][2];
                    let day_dz = da[2][1];

                    let daz_dx = da[0][2];
                    let dax_dz = da[2][0];

                    let day_dx = da[0][1];
                    let dax_dy = da[1][0];

                    let ca = [daz_dy - day_dz, dax_dz - daz_dx, day_dx - dax_dy];

                    // B via biot-savart
                    let mut bx = [0.0];
                    let mut by = [0.0];
                    let mut bz = [0.0];
                    flux_density_linear_filament(
                        (&[*x], &[*y], &[*z]),
                        (&xyz, &xyz, &xyz),
                        (&dlxyz, &dlxyz, &dlxyz),
                        &[1.0],
                        (&mut bx, &mut by, &mut bz),
                    )
                    .unwrap();

                    assert!(approx(bx[0], ca[0], 1e-6, 1e-15));
                    assert!(approx(by[0], ca[1], 1e-6, 1e-15));
                    assert!(approx(bz[0], ca[2], 1e-6, 1e-15));
                }
            }
        }
    }

    /// Check that parallel variants of functions produce the same result as serial.
    /// This also incidentally tests defensive zeroing of input slices.
    #[test]
    fn test_serial_vs_parallel() {
        const NFIL: usize = 10;
        const NOBS: usize = 100;

        // Build a scattering of filament locations
        let xfil: Vec<f64> = (0..NFIL).map(|i| (i as f64).sin()).collect();
        let yfil: Vec<f64> = (0..NFIL).map(|i| (i as f64).cos()).collect();
        let zfil: Vec<f64> = (0..NFIL)
            .map(|i| (i as f64) - (NFIL as f64) / 2.0)
            .collect();
        let xyzfil = (&xfil[..=NFIL - 2], &yfil[..=NFIL - 2], &zfil[..=NFIL - 2]);

        let dlxfil: Vec<f64> = (0..=NFIL - 2).map(|i| xfil[i + 1] - xfil[i]).collect();
        let dlyfil: Vec<f64> = (0..=NFIL - 2).map(|i| yfil[i + 1] - yfil[i]).collect();
        let dlzfil: Vec<f64> = (0..=NFIL - 2).map(|i| zfil[i + 1] - zfil[i]).collect();
        let dlxyzfil = (&dlxfil[..], &dlyfil[..], &dlzfil[..]);

        let ifil: &[f64] = &(0..NFIL - 1).map(|i| (i as f64)).collect::<Vec<f64>>()[..];

        // Build a scattering of observation locations
        let xp: Vec<f64> = (0..NOBS).map(|i| 2.0 * (i as f64).sin() + 2.1).collect();
        let yp: Vec<f64> = (0..NOBS).map(|i| 4.0 * (2.0 * i as f64).cos()).collect();
        let zp: Vec<f64> = (0..NOBS).map(|i| (0.1 * i as f64).exp()).collect();
        let xyzp = (&xp[..], &yp[..], &zp[..]);

        // Some output storage
        // Initialize with different values for each buffer to test zeroing
        let out0 = &mut [0.0; NOBS];
        let out1 = &mut [1.0; NOBS];
        let out2 = &mut [2.0; NOBS];
        let out3 = &mut [3.0; NOBS];
        let out4 = &mut [4.0; NOBS];
        let out5 = &mut [5.0; NOBS];

        // Flux density
        flux_density_linear_filament(xyzp, xyzfil, dlxyzfil, ifil, (out0, out1, out2)).unwrap();
        flux_density_linear_filament_par(xyzp, xyzfil, dlxyzfil, ifil, (out3, out4, out5)).unwrap();
        for i in 0..NOBS {
            assert_eq!(out0[i], out3[i]);
            assert_eq!(out1[i], out4[i]);
            assert_eq!(out2[i], out5[i]);
        }

        // Reinit to test zeroing
        let out0 = &mut [0.0; NOBS];
        let out1 = &mut [1.0; NOBS];
        let out2 = &mut [2.0; NOBS];
        let out3 = &mut [3.0; NOBS];
        let out4 = &mut [4.0; NOBS];
        let out5 = &mut [5.0; NOBS];

        // Vector potential
        vector_potential_linear_filament(xyzp, xyzfil, dlxyzfil, ifil, (out0, out1, out2)).unwrap();
        vector_potential_linear_filament_par(xyzp, xyzfil, dlxyzfil, ifil, (out3, out4, out5))
            .unwrap();
        for i in 0..NOBS {
            assert_eq!(out0[i], out3[i]);
            assert_eq!(out1[i], out4[i]);
            assert_eq!(out2[i], out5[i]);
        }
    }
}
