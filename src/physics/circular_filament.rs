//! Magnetics calculations for circular current filaments.
use std::num::NonZeroUsize;

use rayon::{
    iter::{IndexedParallelIterator, ParallelIterator},
    slice::{ParallelSlice, ParallelSliceMut},
};

use crate::math::{ellipe, ellipk};

use crate::{MU0_OVER_4PI, MU_0};

/// Flux contributions from some circular filaments to some observation points, which happens to be
/// the Green's function for the Grad-Shafranov elliptic operator, $\Delta^{\*}$.
/// This variant of the function is parallelized over chunks of observation points.
///
/// # Arguments
///
/// * `ifil`: (A) current in each filament, length `m`
/// * `rfil`:    (m) r-coord of each filament, length `m`
/// * `zfil`:    (m) z-coord of each filament, length `m`
/// * `rprime`:  (m) r-coord of each observation point, length `n`
/// * `zprime`:  (m) z-coord of each observation point, length `n`
/// * `out_psi`: (Wb) or (H-A) or (T-m^2) or (V-s), poloidal flux at observation locations, length `n`
///
/// # Commentary
///
/// Represents contribution from a current at (R, Z) to an observation point at (Rprime, Zprime)
///
/// Note Jardin's 4.61-4.66 presents it with a different definition of
/// the elliptic integrals from what is used here and in scipy.
///
/// # References
///
///   \[1\] D. Kaltsas, A. Kuiroukidis, and G. Throumoulopoulos, “A tokamak pertinent analytic equilibrium with plasma flow of arbitrary direction,”
///         Physics of Plasmas, vol. 26, p. 124501, Dec. 2019,
///         doi: [10.1063/1.5120341](https://doi.org/10.1063/1.5120341).
///
///   \[2\] S. Jardin, *Computational Methods in Plasma Physics*, 1st ed. USA: CRC Press, Inc., 2010.
///
///   \[3\] J. Huang and J. Menard, “Development of an Auto-Convergent Free-Boundary Axisymmetric Equilibrium Solver,”
///         Journal of Undergraduate Research, vol. 6, Jan. 2006, Accessed: May 05, 2021. \[Online\].
///         Available: <https://www.osti.gov/biblio/1051805-development-auto-convergent-free-boundary-axisymmetric-equilibrium-solver>
///
///   \[4\] J. C. Simpson, J. E. Lane, C. D. Immer, R. C. Youngquist, and T. Steinrock,
///         “Simple Analytic Expressions for the Magnetic Field of a Circular Current Loop,”
///         Jan. 01, 2001. Accessed: Sep. 06, 2022. [Online]. Available: <https://ntrs.nasa.gov/citations/20010038494>
pub fn flux_circular_filament_par(
    ifil: &[f64],
    rfil: &[f64],
    zfil: &[f64],
    rprime: &[f64],
    zprime: &[f64],
    out_psi: &mut [f64],
) -> Result<(), &'static str> {
    // Chunk inputs
    let ncores = std::thread::available_parallelism()
        .unwrap_or(NonZeroUsize::MIN)
        .get();

    let n = (rprime.len() / ncores).max(1);

    let rprimec = rprime.par_chunks(n);
    let zprimec = zprime.par_chunks(n);

    let outc = out_psi.par_chunks_mut(n);

    // Run calcs
    outc.zip(rprimec.zip(zprimec))
        .try_for_each(|(outc, (rc, zc))| flux_circular_filament(ifil, rfil, zfil, rc, zc, outc))?;

    Ok(())
}

/// Flux contributions from some circular filaments to some observation points, which happens to be
/// the Green's function for the Grad-Shafranov elliptic operator, $\Delta^{\*}$.
///
/// # Arguments
///
/// * `ifil`: (A) current in each filament, length `m`
/// * `rfil`:    (m) r-coord of each filament, length `m`
/// * `zfil`:    (m) z-coord of each filament, length `m`
/// * `rprime`:  (m) r-coord of each observation point, length `n`
/// * `zprime`:  (m) z-coord of each observation point, length `n`
/// * `out_psi`: (Wb) or (H-A) or (T-m^2) or (V-s), poloidal flux at observation locations, length `n`
///
/// # Commentary
///
/// Represents contribution from a current at (R, Z) to an observation point at (Rprime, Zprime)
///
/// Note Jardin's 4.61-4.66 presents it with a different definition of
/// the elliptic integrals from what is used here and in scipy.
///
/// # References
///
///   \[1\] D. Kaltsas, A. Kuiroukidis, and G. Throumoulopoulos, “A tokamak pertinent analytic equilibrium with plasma flow of arbitrary direction,”
///         Physics of Plasmas, vol. 26, p. 124501, Dec. 2019,
///         doi: [10.1063/1.5120341](https://doi.org/10.1063/1.5120341).
///
///   \[2\] S. Jardin, *Computational Methods in Plasma Physics*, 1st ed. USA: CRC Press, Inc., 2010.
///
///   \[3\] J. Huang and J. Menard, “Development of an Auto-Convergent Free-Boundary Axisymmetric Equilibrium Solver,”
///         Journal of Undergraduate Research, vol. 6, Jan. 2006, Accessed: May 05, 2021. \[Online\].
///         Available: <https://www.osti.gov/biblio/1051805-development-auto-convergent-free-boundary-axisymmetric-equilibrium-solver>
///
///   \[4\] J. C. Simpson, J. E. Lane, C. D. Immer, R. C. Youngquist, and T. Steinrock,
///         “Simple Analytic Expressions for the Magnetic Field of a Circular Current Loop,”
///         Jan. 01, 2001. Accessed: Sep. 06, 2022. [Online]. Available: <https://ntrs.nasa.gov/citations/20010038494>
pub fn flux_circular_filament(
    ifil: &[f64],
    rfil: &[f64],
    zfil: &[f64],
    rprime: &[f64],
    zprime: &[f64],
    out_psi: &mut [f64],
) -> Result<(), &'static str> {
    // Check lengths; Error if they do not match
    let m: usize = ifil.len();
    let n: usize = rprime.len();
    if rfil.len() != m || zfil.len() != m || zprime.len() != n || out_psi.len() != n {
        return Err("Length mismatch");
    }

    for i in 0..n {
        out_psi[i] = 0.0;
    }

    for i in 0..n {
        for j in 0..m {
            let rrprime = rfil[j] * rprime[i];
            let r_plus_rprime = rfil[j] + rprime[i];
            let z_minus_zprime = zfil[j] - zprime[i];
            let k2 = 4.0 * rrprime / (r_plus_rprime.powi(2) + z_minus_zprime.powi(2));

            // The contributions added here have units of ampere meter,
            // but the sum will have units of weber after it is multiplied by mu_0 below.
            out_psi[i] +=
                ifil[j] * (rrprime / k2).sqrt() * ((2.0 - k2) * ellipk(k2) - 2.0 * ellipe(k2));
        }
    }

    for i in 0..n {
        out_psi[i] *= MU_0;
    }

    Ok(())
}

/// Off-axis Br,Bz components for a circular current filament in vacuum.
/// This variant of the function is parallelized over chunks of observation points.
///
/// # Arguments
///
/// * `ifil`:    (A) current in each filament, length `m`
/// * `rfil`:    (m) r-coord of each filament, length `m`
/// * `zfil`:    (m) z-coord of each filament, length `m`
/// * `rprime`:  (m) r-coord of each observation point, length `n`
/// * `zprime`:  (m) z-coord of each observation point, length `n`
/// * `out_r`:   (T), r-component of magnetic flux density at observation locations, length `n`
/// * `out_z`:   (T), z-component of magnetic flux density at observation locations, length `n`
///
/// # Commentary
///
/// Near-exact formula (except numerically-evaluated elliptic integrals).
/// See eqns. 12,13 pg. 34 in \[1\], eqn 9.8.7 in \[2\], and all of \[3\].
///
/// Note the formula for Br as given by \[1\] is incorrect and does not satisfy the
/// constraints of the calculation without correcting by a factor of (z / r).
///
/// # References
///
///   \[1\] D. B. Montgomery and J. Terrell,
///         “Some Useful Information For The Design Of Aircore Solenoids,
///         Part I. Relationships Between Magnetic Field, Power, Ampere-Turns
///         And Current Density. Part II. Homogeneous Magnetic Fields,”
///         Massachusetts Inst. Of Tech. Francis Bitter National Magnet Lab, Cambridge, MA,
///         Nov. 1961. Accessed: May 18, 2021. \[Online\].
///         Available: <https://apps.dtic.mil/sti/citations/tr/AD0269073>
///
///   \[2\] 8.02 Course Notes. Available: <https://web.mit.edu/8.02t/www/802TEAL3D/visualizations/coursenotes/modules/guide09.pdf>
///
///   \[3\] Eric Dennyson, "Magnet Formulas". Available: <https://tiggerntatie.github.io/emagnet-py/offaxis/off_axis_loop.html>
///
///   \[4\] J. C. Simpson, J. E. Lane, C. D. Immer, R. C. Youngquist, and T. Steinrock,
///         “Simple Analytic Expressions for the Magnetic Field of a Circular Current Loop,”
///         Jan. 01, 2001. Accessed: Sep. 06, 2022. [Online]. Available: <https://ntrs.nasa.gov/citations/20010038494>
pub fn flux_density_circular_filament_par(
    ifil: &[f64],
    rfil: &[f64],
    zfil: &[f64],
    rprime: &[f64],
    zprime: &[f64],
    out_r: &mut [f64],
    out_z: &mut [f64],
) -> Result<(), &'static str> {
    // Chunk inputs
    let ncores = std::thread::available_parallelism()
        .unwrap_or(NonZeroUsize::MIN)
        .get();

    let n = (rprime.len() / ncores).max(1);

    let rprimec = rprime.par_chunks(n);
    let zprimec = zprime.par_chunks(n);

    let outrc = out_r.par_chunks_mut(n);
    let outzc = out_z.par_chunks_mut(n);

    // Run calcs
    outrc
        .zip(outzc.zip(rprimec.zip(zprimec)))
        .try_for_each(|(orc, (ozc, (rc, zc)))| {
            flux_density_circular_filament(ifil, rfil, zfil, rc, zc, orc, ozc)
        })?;

    Ok(())
}

/// Off-axis Br,Bz components for a circular current filament in vacuum.
///
/// # Arguments
///
/// * `ifil`:    (A) current in each filament, length `m`
/// * `rfil`:    (m) r-coord of each filament, length `m`
/// * `zfil`:    (m) z-coord of each filament, length `m`
/// * `rprime`:  (m) r-coord of each observation point, length `n`
/// * `zprime`:  (m) z-coord of each observation point, length `n`
/// * `out_r`:   (T), r-component of magnetic flux density at observation locations, length `n`
/// * `out_z`:   (T), z-component of magnetic flux density at observation locations, length `n`
///
/// # Commentary
///
/// Near-exact formula (except numerically-evaluated elliptic integrals).
/// See eqns. 12,13 pg. 34 in \[1\], eqn 9.8.7 in \[2\], and all of \[3\].
///
/// Note the formula for Br as given by \[1\] is incorrect and does not satisfy the
/// constraints of the calculation without correcting by a factor of (z / r).
///
/// # References
///
///   \[1\] D. B. Montgomery and J. Terrell,
///         “Some Useful Information For The Design Of Aircore Solenoids,
///         Part I. Relationships Between Magnetic Field, Power, Ampere-Turns
///         And Current Density. Part II. Homogeneous Magnetic Fields,”
///         Massachusetts Inst. Of Tech. Francis Bitter National Magnet Lab, Cambridge, MA,
///         Nov. 1961. Accessed: May 18, 2021. \[Online\].
///         Available: <https://apps.dtic.mil/sti/citations/tr/AD0269073>
///
///   \[2\] 8.02 Course Notes. Available: <https://web.mit.edu/8.02t/www/802TEAL3D/visualizations/coursenotes/modules/guide09.pdf>
///
///   \[3\] Eric Dennyson, "Magnet Formulas". Available: <https://tiggerntatie.github.io/emagnet-py/offaxis/off_axis_loop.html>
///
///   \[4\] J. C. Simpson, J. E. Lane, C. D. Immer, R. C. Youngquist, and T. Steinrock,
///         “Simple Analytic Expressions for the Magnetic Field of a Circular Current Loop,”
///         Jan. 01, 2001. Accessed: Sep. 06, 2022. [Online]. Available: <https://ntrs.nasa.gov/citations/20010038494>
pub fn flux_density_circular_filament(
    ifil: &[f64],
    rfil: &[f64],
    zfil: &[f64],
    rprime: &[f64],
    zprime: &[f64],
    out_r: &mut [f64],
    out_z: &mut [f64],
) -> Result<(), &'static str> {
    let n = ifil.len();
    let m = rprime.len();

    // Check lengths; Error if they do not match
    if rfil.len() != n
        || zfil.len() != n
        || zprime.len() != m
        || out_r.len() != m
        || out_z.len() != m
    {
        return Err("Length mismatch");
    }

    for j in 0..m {
        out_r[j] = 0.0;
        out_z[j] = 0.0;
    }

    // There aren't necessarily more observation points or filaments, depending on the use case.
    // The more common extreme is to see a very large number of filaments evaluated at a smaller
    // number of observation points. However, this particular calc suffers badly when iterating
    // over observation points first, so to capture a 50% speedup for cases with >=10 observation
    // points at the expense of a 30% slowdown for evaluating single observation points, we
    // iterate over filaments first here.
    for i in 0..n {
        for j in 0..m {
            let z = zprime[j] - zfil[i]; // [m]

            let z2 = z * z; // [m^2]
            let r2 = rprime[j] * rprime[j]; // [m^2]

            let rpr = rfil[i] + rprime[j];

            let q = rpr.mul_add(rpr, z2); // [m^2]
            let k2 = 4.0 * rfil[i] * rprime[j] / q; // [nondim]

            let a0 = 2.0 * ifil[i] / q.sqrt(); // [A/m]

            let f = ellipk(k2); // [nondim]
            let s = ellipe(k2) / (1.0 - k2); // [nondim]

            // Bake some reusable values
            let s_over_q = s / q; // [m^-2]
            let rfil2 = rfil[i] * rfil[i]; // [m^2]

            // Magnetic field intensity, less the factor of 4pi that we have adjusted out of mu_0
            let hr = (z / rprime[j]) * a0 * s_over_q.mul_add(rfil2 + r2 + z2, -f);
            let hz = a0 * s_over_q.mul_add(rfil2 - r2 - z2, f);

            // Magnetic flux density assuming vacuum permeability
            // The contributions added here have units of ampere per meter,
            // but the result will have units of tesla after it is multiplied by mu_0 / (4 * pi) below.
            out_r[j] += hr;
            out_z[j] += hz;
        }
    }

    for j in 0..m {
        out_r[j] *= MU0_OVER_4PI;
        out_z[j] *= MU0_OVER_4PI;
    }

    Ok(())
}

/// Off-axis A_phi component for a circular current filament in vacuum.
/// This variant of the function is parallelized over chunks of observation points.
///
/// # Arguments
///
/// * `ifil`:    (A) current in each filament, length `m`
/// * `rfil`:    (m) r-coord of each filament, length `m`
/// * `zfil`:    (m) z-coord of each filament, length `m`
/// * `rprime`:  (m) r-coord of each observation point, length `n`
/// * `zprime`:  (m) z-coord of each observation point, length `n`
/// * `out_phi`: (V-s/m), phi-component of magnetic vector potential at observation locations, length `n`
///
/// # Commentary
///
/// Near-exact formula (except numerically-evaluated elliptic integrals).
/// The vector potential of a loop has zero r- and z- components due to symmetry,
/// and does not vary in the phi-direction.
///
/// # References
///
///   \[1\] J. C. Simpson, J. E. Lane, C. D. Immer, R. C. Youngquist, and T. Steinrock,
///         “Simple Analytic Expressions for the Magnetic Field of a Circular Current Loop,”
///         Jan. 01, 2001. Accessed: Sep. 06, 2022. [Online]. Available: <https://ntrs.nasa.gov/citations/20010038494>
pub fn vector_potential_circular_filament_par(
    ifil: &[f64],
    rfil: &[f64],
    zfil: &[f64],
    rprime: &[f64],
    zprime: &[f64],
    out_phi: &mut [f64],
) -> Result<(), &'static str> {
    // Chunk inputs
    let ncores = std::thread::available_parallelism()
        .unwrap_or(NonZeroUsize::MIN)
        .get();

    let n = (rprime.len() / ncores).max(1);

    let rprimec = rprime.par_chunks(n);
    let zprimec = zprime.par_chunks(n);

    let outc = out_phi.par_chunks_mut(n);

    // Run calcs
    outc.zip(rprimec.zip(zprimec))
        .try_for_each(|(outc, (rc, zc))| {
            vector_potential_circular_filament(ifil, rfil, zfil, rc, zc, outc)
        })?;

    Ok(())
}

/// Off-axis A_phi component for a circular current filament in vacuum.
///
/// # Arguments
///
/// * `ifil`:    (A) current in each filament, length `m`
/// * `rfil`:    (m) r-coord of each filament, length `m`
/// * `zfil`:    (m) z-coord of each filament, length `m`
/// * `rprime`:  (m) r-coord of each observation point, length `n`
/// * `zprime`:  (m) z-coord of each observation point, length `n`
/// * `out_phi`: (V-s/m), phi-component of magnetic vector potential at observation locations, length `n`
///
/// # Commentary
///
/// Near-exact formula (except numerically-evaluated elliptic integrals).
/// The vector potential of a loop has zero r- and z- components due to symmetry,
/// and does not vary in the phi-direction.
///
/// # References
///
///   \[1\] J. C. Simpson, J. E. Lane, C. D. Immer, R. C. Youngquist, and T. Steinrock,
///         “Simple Analytic Expressions for the Magnetic Field of a Circular Current Loop,”
///         Jan. 01, 2001. Accessed: Sep. 06, 2022. [Online]. Available: <https://ntrs.nasa.gov/citations/20010038494>
pub fn vector_potential_circular_filament(
    ifil: &[f64],
    rfil: &[f64],
    zfil: &[f64],
    rprime: &[f64],
    zprime: &[f64],
    out_phi: &mut [f64],
) -> Result<(), &'static str> {
    let n = ifil.len();
    let m = rprime.len();

    // Check lengths; Error if they do not match
    if rfil.len() != n || zfil.len() != n || zprime.len() != m || out_phi.len() != m {
        return Err("Length mismatch");
    }

    for j in 0..m {
        out_phi[j] = 0.0;
    }

    for i in 0..n {
        for j in 0..m {
            // Eq. 1 and 2 of Simpson2001 give a formula for the vector potential of a loop in spherical coordinates.
            // Here, we use that formula adjusted to cylindrical coordinates.
            // r_spherical*sin(theta) = r_cylindrical
            // r_spherical^2 = r_cylindrical^2 + z^2
            let z = zprime[j] - zfil[i]; // [m]

            // Assemble argument to elliptic integrals
            let rpr = rfil[i] + rprime[j];
            let denom = z.mul_add(z, rpr * rpr);
            let numer = 4.0 * rfil[i] * rprime[j];
            let k2 = numer / denom;

            // Elliptic integral terms
            let c0 = ((2.0 - k2) * ellipk(k2) - 2.0 * ellipe(k2)) / k2;

            // Factor multiplied into elliptic integral terms
            let c1 = ifil[i] * 4.0 * rfil[i] / denom.sqrt();

            out_phi[j] = c0.mul_add(c1, out_phi[j]);
        }
    }

    for j in 0..m {
        out_phi[j] *= MU0_OVER_4PI;
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
    /// and that psi = integral(dot(A, dL)) =  2pi * r * a
    #[test]
    fn test_vector_potential() {
        let rfil = 1.0 / core::f64::consts::PI; // [m] some number
        let zfil = 1.0 / core::f64::consts::E; // [m] some number

        let vp = |r: f64, z: f64| {
            let mut out = [0.0];

            vector_potential_circular_filament(&[1.0], &[rfil], &[zfil], &[r], &[z], &mut out)
                .unwrap();

            out[0]
        };

        let zvals = [0.25, 0.5, 2.5, 10.0, 0.0, -10.0, -2.5, -0.5, -0.25];
        let rvals = [0.25, 0.5, 2.5, 10.0];
        // finite diff delta needs to be small enough to be accurate
        // but large enough that we can tell the difference between adjacent points
        // that are very far from the origin
        let eps = 1e-7;
        for r in rvals.iter() {
            for z in zvals.iter() {
                // Finite-difference curl of the vector potential in cylindrical coordinates.
                // The radial and z components of the vector potential are zero.
                let mut ca = [0.0; 3];
                // curl(A)[0] = - d(A_phi) / dz
                let a0 = vp(*r, *z - eps);
                let a1 = vp(*r, *z + eps);
                ca[0] = -(a1 - a0) / (2.0 * eps);
                // curl(A)[2] = (1 / rho ) d(rho A_phi) / d(rho)
                let ra0 = (*r - eps) * vp(*r - eps, *z);
                let ra1 = (*r + eps) * vp(*r + eps, *z);
                ca[2] = (ra1 - ra0) / (2.0 * eps) / *r;

                // B via biot-savart
                let mut br = [0.0];
                let mut bz = [0.0];
                flux_density_circular_filament(
                    &[1.0],
                    &[rfil],
                    &[zfil],
                    &[*r],
                    &[*z],
                    &mut br,
                    &mut bz,
                )
                .unwrap();

                assert!(approx(br[0], ca[0], 1e-7, 1e-13));
                assert!(approx(bz[0], ca[2], 1e-7, 1e-13));

                // Flux via analytic formula
                // psi = integral(dot(A, dL)) =  2pi * r * a
                let psi_from_a = 2.0 * PI * *r * vp(*r, *z);
                let mut psi = [0.0];
                flux_circular_filament(&[1.0], &[rfil], &[zfil], &[*r], &[*z], &mut psi).unwrap();
                println!("{psi:?}, {psi_from_a}");
                assert!(approx(psi_from_a, psi[0], 1e-10, 0.0)); // Should be very close to float roundoff
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
        let rfil: Vec<f64> = (0..NFIL).map(|i| (i as f64).sin() + 1.2).collect();
        let zfil: Vec<f64> = (0..NFIL)
            .map(|i| (i as f64) - (NFIL as f64) / 2.0)
            .collect();
        let ifil: Vec<f64> = (0..NFIL).map(|i| (i as f64)).collect();

        // Build a scattering of observation locations
        let rprime: Vec<f64> = (0..NOBS).map(|i| 2.0 * (i as f64).sin() + 2.1).collect();
        let zprime: Vec<f64> = (0..NOBS).map(|i| 4.0 * (2.0 * i as f64).cos()).collect();

        // Some output storage
        // Initialize with different values for each buffer to test zeroing
        let out0 = &mut [0.0; NOBS];
        let out1 = &mut [1.0; NOBS];
        let out2 = &mut [2.0; NOBS];
        let out3 = &mut [3.0; NOBS];

        // Flux
        flux_circular_filament(&ifil, &rfil, &zfil, &rprime, &zprime, out0).unwrap();
        flux_circular_filament_par(&ifil, &rfil, &zfil, &rprime, &zprime, out1).unwrap();
        for i in 0..NOBS {
            assert_eq!(out0[i], out1[i]);
        }

        // Flux density
        flux_density_circular_filament(&ifil, &rfil, &zfil, &rprime, &zprime, out0, out1).unwrap();
        flux_density_circular_filament_par(&ifil, &rfil, &zfil, &rprime, &zprime, out2, out3)
            .unwrap();
        for i in 0..NOBS {
            assert_eq!(out0[i], out2[i]);
            assert_eq!(out1[i], out3[i]);
        }

        // Vector potential
        let out0 = &mut [0.0; NOBS]; // Reinit with different values to test zeroing
        let out1 = &mut [1.0; NOBS];
        vector_potential_circular_filament(&ifil, &rfil, &zfil, &rprime, &zprime, out0).unwrap();
        vector_potential_circular_filament_par(&ifil, &rfil, &zfil, &rprime, &zprime, out1)
            .unwrap();
        for i in 0..NOBS {
            assert_eq!(out0[i], out1[i]);
        }
    }
}
