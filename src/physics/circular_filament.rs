//! Magnetics calculations for circular current filaments.
use std::num::NonZeroUsize;

use rayon::{
    iter::{IndexedParallelIterator, ParallelIterator},
    slice::{ParallelSlice, ParallelSliceMut},
};

use crate::math::{ellipe, ellipk, rss3};

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
            // The inner function is inlined, so values that are reused between iterations
            // are pulled to the outer scope by the compiler and do not affect performance
            out_psi[i] +=
                flux_circular_filament_scalar(ifil[j], rfil[j], zfil[j], rprime[i], zprime[i]);
        }
    }

    Ok(())
}

/// Flux contributions from some circular filaments to some observation points, which happens to be
/// the Green's function for the Grad-Shafranov elliptic operator, $\Delta^{\*}$.
///
/// # Arguments
///
/// * `ifil`: (A) current in filament
/// * `rfil`:    (m) r-coord of filament
/// * `zfil`:    (m) z-coord of filament
/// * `rprime`:  (m) r-coord of observation point
/// * `zprime`:  (m) z-coord of observation point
///
/// # Returns
///
/// * `psi`: (Wb) or (H-A) or (T-m^2) or (V-s), poloidal flux at observation location
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
#[inline]
pub fn flux_circular_filament_scalar(
    ifil: f64,
    rfil: f64,
    zfil: f64,
    rprime: f64,
    zprime: f64,
) -> f64 {
    let rrprime = rfil * rprime;
    let r_plus_rprime = rfil + rprime;
    let z_minus_zprime = zfil - zprime;
    let k2 = 4.0 * rrprime / (r_plus_rprime.powi(2) + z_minus_zprime.powi(2));
    let psi = MU_0 * ifil * (rrprime / k2).sqrt() * ((2.0 - k2) * ellipk(k2) - 2.0 * ellipe(k2)); // [V-s]
    psi
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
            // The inner function is inlined, so values that are reused between iterations
            // are pulled to the outer scope by the compiler and do not affect performance
            let (br, bz) = flux_density_circular_filament_scalar(
                ifil[i], rfil[i], zfil[i], rprime[j], zprime[j],
            );
            out_r[j] += br;
            out_z[j] += bz;
        }
    }

    Ok(())
}

/// Off-axis Br,Bz components for a circular current filament in vacuum.
///
/// # Arguments
///
/// * `ifil`:    (A) current in filament
/// * `rfil`:    (m) r-coord of filament
/// * `zfil`:    (m) z-coord of filament
/// * `rprime`:  (m) r-coord of observation point
/// * `zprime`:  (m) z-coord of observation point
///
/// # Returns
///
/// * `br`:   (T), r-component of magnetic flux density at observation location
/// * `bz`:   (T), z-component of magnetic flux density at observation location
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
#[inline]
pub fn flux_density_circular_filament_scalar(
    ifil: f64,
    rfil: f64,
    zfil: f64,
    rprime: f64,
    zprime: f64,
) -> (f64, f64) {
    let z = zprime - zfil; // [m]

    let z2 = z * z; // [m^2]
    let r2 = rprime * rprime; // [m^2]

    let rpr = rfil + rprime;

    let q = rpr.mul_add(rpr, z2); // [m^2]
    let k2 = 4.0 * rfil * rprime / q; // [nondim]

    let a0 = 2.0 * ifil / q.sqrt(); // [A/m]

    let f = ellipk(k2); // [nondim]
    let s = ellipe(k2) / (1.0 - k2); // [nondim]

    // Bake some reusable values
    let s_over_q = s / q; // [m^-2]
    let rfil2 = rfil * rfil; // [m^2]

    // Magnetic field intensity, less the factor of 4pi that we have adjusted out of mu_0
    let hr = (z / rprime) * a0 * s_over_q.mul_add(rfil2 + r2 + z2, -f);
    let hz = a0 * s_over_q.mul_add(rfil2 - r2 - z2, f);

    // Magnetic flux density assuming vacuum permeability
    let br = MU0_OVER_4PI * hr;
    let bz = MU0_OVER_4PI * hz;

    (br, bz)
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
            // The inner function is inlined, so values that are reused between iterations
            // are pulled to the outer scope by the compiler and do not affect performance
            out_phi[j] += vector_potential_circular_filament_scalar(
                ifil[i], rfil[i], zfil[i], rprime[j], zprime[j],
            );
        }
    }

    Ok(())
}

/// Off-axis A_phi component for a circular current filament in vacuum.
///
/// # Arguments
///
/// * `ifil`:    (A) current in filament
/// * `rfil`:    (m) r-coord of filament
/// * `zfil`:    (m) z-coord of filament
/// * `rprime`:  (m) r-coord of observation point
/// * `zprime`:  (m) z-coord of observation point
///
/// # Returns
/// * `a_phi`: (V-s/m), phi-component of magnetic vector potential at observation location
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
#[inline]
pub fn vector_potential_circular_filament_scalar(
    ifil: f64,
    rfil: f64,
    zfil: f64,
    rprime: f64,
    zprime: f64,
) -> f64 {
    // Eq. 1 and 2 of Simpson2001 give a formula for the vector potential of a loop in spherical coordinates.
    // Here, we use that formula adjusted to cylindrical coordinates.
    // r_spherical*sin(theta) = r_cylindrical
    // r_spherical^2 = r_cylindrical^2 + z^2
    let z = zprime - zfil; // [m]

    // Assemble argument to elliptic integrals
    let rpr2 = (rfil + rprime).powf(2.0);
    let denom = z.mul_add(z, rpr2);
    let numer = 4.0 * rfil * rprime;
    let k2 = numer / denom;

    // Elliptic integral terms
    let c0 = ((2.0 - k2) * ellipk(k2) - 2.0 * ellipe(k2)) / k2;

    // Factor multiplied into elliptic integral terms
    let c1 = MU0_OVER_4PI * ifil * 4.0 * rfil / denom.sqrt();

    let a_phi = c0 * c1; // [V-s/m] phi-component of vector potential
    a_phi // Other components are zero
}

/// Mutual inductance between a circular filament and a linear filament.
/// This method is much faster (~100x typically) than discretizing the circular loop
/// into linear segments and using Neumann's formula.
///
/// # Arguments
///
/// * `rznfil`:    (m, m, nondim) r,z-coord and number of turns of circular filament
/// * `xyzfil0`:   (m) (x, y, z) coordinates of start of linear segment
/// * `xyzfil1`:   (m) (x, y, z) coordinates of end of linear segment
///
/// # Returns
///
/// * `m`: (H), mutual inductance
#[inline]
pub fn mutual_inductance_circular_to_linear_scalar(
    rznfil: (f64, f64, f64),
    xyzfil0: (f64, f64, f64),
    xyzfil1: (f64, f64, f64),
) -> f64 {
    // First, get the filament vector
    let dlxfil = xyzfil1.0 - xyzfil0.0; // [m]
    let dlyfil = xyzfil1.1 - xyzfil0.1;
    let dlzfil = xyzfil1.2 - xyzfil0.2;
    // Next, we need to map the linear filament into cylindrical coordinates
    //    r = (x^2 + y^2)^0.5 in cylindrical
    let path_r = rss3(xyzfil0.0, xyzfil0.1, 0.0); // [m]
    let path_dr = rss3(dlxfil, dlyfil, 0.0); // [m]

    //    phi = tan^-1(y/x)
    let path_phi0 = f64::atan2(xyzfil0.1, xyzfil0.0);
    let path_phi1 = f64::atan2(xyzfil1.1, xyzfil1.0);
    let path_dphi = path_phi1 - path_phi0;

    //    midpoint is best for capturing curvature in piecewise-linear paths properly
    let path_r_mid = path_r + path_dr / 2.0; // [m]
    let path_z_mid = xyzfil0.2 + dlzfil / 2.0; // [m]
    let path_dlphi = path_r_mid * path_dphi; // [m] length in phi-direction; 2*pi cancels out

    // Get cylindrical vector potential at linear segment midpoint
    // for a unit current, which is equivalent to mutual inductance per unit length
    // [H/m]
    let a_phi_per_A = rznfil.2
        * vector_potential_circular_filament_scalar(
            1.0, rznfil.0, rznfil.1, path_r_mid, path_z_mid,
        );

    // Recover mutual inductance as dot(A, dL)/I
    let m = a_phi_per_A * path_dlphi; // [H]
    m
}

/// Mutual inductance between a collection of circular filaments and a piecewise-linear filament.
/// This method is much faster (~100x typically) than discretizing the circular loop
/// into linear segments and using Neumann's formula.
///
/// # Arguments
///
/// * `rznfil`:  (m, m, nondim) r,z-coord and number of turns of each circular filament, length `m`
/// * `xyzfil`:  (m) filament origin coordinates for linear path, length `n`, including endpoint
///
/// # Returns
///
/// * `m`: (V-s/m), phi-component of magnetic vector potential at observation locations
pub fn mutual_inductance_circular_to_linear(
    rznfil: (&[f64], &[f64], &[f64]),
    xyzfil: (&[f64], &[f64], &[f64]),
) -> Result<f64, &'static str> {
    // Check lengths; Error if they do not match
    let n = xyzfil.0.len();
    if xyzfil.0.len() != n || xyzfil.1.len() != n || xyzfil.2.len() != n || n < 2
    // Need at least 2 points to form a piecewise linear path
    {
        return Err("Input length mismatch");
    }

    // Check lengths; Error if they do not match
    let m = rznfil.0.len();
    if rznfil.0.len() != m || rznfil.1.len() != m || rznfil.2.len() != m {
        return Err("Length mismatch");
    }

    let mut mutual_inductance = 0.0;

    for i in 0..n - 1 {
        for j in 0..m {
            // The inner function is inlined, so values that are reused between iterations
            // are pulled to the outer scope by the compiler and do not affect performance
            let xyzfil0 = (xyzfil.0[i], xyzfil.1[i], xyzfil.2[i]);
            let xyzfil1 = (xyzfil.0[i + 1], xyzfil.1[i + 1], xyzfil.2[i + 1]);
            mutual_inductance += mutual_inductance_circular_to_linear_scalar(
                (rznfil.0[j], rznfil.1[j], rznfil.2[j]),
                xyzfil0,
                xyzfil1,
            );
        }
    }

    Ok(mutual_inductance)
}

pub fn mutual_inductance_circular_to_linear_par(
    rznfil: (&[f64], &[f64], &[f64]),
    xyzfil: (&[f64], &[f64], &[f64]),
) -> Result<f64, &'static str> {
    // Unpack
    let (rfil, zfil, nfil) = rznfil;

    // Chunk inputs
    let ncores = std::thread::available_parallelism()
        .unwrap_or(NonZeroUsize::MIN)
        .get();

    let n = (rfil.len() / ncores).max(1);

    let rfilc = rfil.par_chunks(n);
    let zfilc = zfil.par_chunks(n);
    let nfilc = nfil.par_chunks(n);

    // Run calcs
    // We have to sum over contributions that are each individually fallible,
    // which results in a bit of clutter with the fold-reduce pattern
    let mutual_inductance = nfilc
        .zip(rfilc.zip(zfilc))
        .try_fold(
            || 0.0,
            |acc, (nc, (rc, zc))| {
                let m_contrib = mutual_inductance_circular_to_linear((rc, zc, nc), xyzfil)?;
                Ok::<f64, &'static str>(acc + m_contrib)
            },
        )
        .try_reduce(|| 0.0, |acc, v| Ok(acc + v))?;

    Ok(mutual_inductance)
}

#[cfg(test)]
mod test {
    use core::f64::consts::{E, PI};

    use super::*;

    /// Div/0-resistant approximate comparison
    fn approx(truth: f64, val: f64, rtol: f64, atol: f64) -> bool {
        let abs_err = (val - truth).abs();
        let lim = rtol * truth.abs() + atol;
        abs_err < lim
    }

    /// Make sure the circular-to-linear mutual inductance calc matches
    /// the result achieved by discretizing the circular filament
    /// into linear segments
    #[test]
    fn test_mutual_inductance_to_linear() {
        let linspace = |start, end, n| {
            (0..n)
                .map(|i| start + (i as f64 / (n - 1) as f64) * (end - start))
                .collect::<Vec<f64>>()
        };

        let diff = |v: &[f64]| {
            v[1..]
                .iter()
                .zip(v[0..v.len() - 1].iter())
                .map(|(&b, &a)| b - a)
                .collect::<Vec<f64>>()
        };

        let discretize_circular_filament = |r: f64, z, ndiscr| {
            let x: Vec<f64> = linspace(0.0, 2.0 * PI, ndiscr)
                .iter()
                .map(|v| r * v.cos())
                .collect();
            let y: Vec<f64> = linspace(0.0, 2.0 * PI, ndiscr)
                .iter()
                .map(|v| r * v.sin())
                .collect();
            let z: Vec<f64> = (0..ndiscr).map(|_| z).collect();
            (x, y, z)
        };

        // Make some circular filaments
        let r = 1.0 / PI; // [m] some number
        let z = 1.0 / E; // [m] some number

        let rfil = [r, r + E / 4.0];
        let zfil = [z, -z];
        let nfil = [PI, E];

        // Make a slightly tilted helical piecewise-linear filament
        let n = 10_000;
        let xc = [0.1, -0.1]; // Start and end of centerline path
        let yc = [-0.05, 0.2];
        let zc = [-2.0 * z, 2.0 * z];

        let xc: Vec<f64> = linspace(xc[0], xc[1], n);
        let yc: Vec<f64> = linspace(yc[0], yc[1], n);
        let zc: Vec<f64> = linspace(zc[0], zc[1], n);

        let mut x = xc.clone();
        let mut y = xc.clone();
        let mut z = xc.clone();
        crate::mesh::filament_helix_path(
            (&xc, &yc, &zc),
            (2.0 * E / 3.0, 0.0, 0.0),
            0.5,
            0.0,
            (&mut x, &mut y, &mut z),
        )
        .unwrap();
        let dlxfil1 = diff(&x);
        let dlyfil1 = diff(&y);
        let dlzfil1 = diff(&z);
        let dlxyzfil1 = (&dlxfil1[..], &dlyfil1[..], &dlzfil1[..]);

        // Get mutual inductance by purpose-made calc
        // [H]
        let mutual_inductance =
            mutual_inductance_circular_to_linear_par((&rfil, &zfil, &nfil), (&x, &y, &z)).unwrap();

        // Get mutual inductance by brute-force calc
        let mut mutual_inductance_2 = 0.0;
        for i in 0..rfil.len() {
            let ndiscr = 1000;
            let (xfil, yfil, zfil) = discretize_circular_filament(rfil[i], zfil[i], ndiscr);
            let dlxfil0 = diff(&xfil);
            let dlyfil0 = diff(&yfil);
            let dlzfil0 = diff(&zfil);
            let dlxyzfil0 = (&dlxfil0[..], &dlyfil0[..], &dlzfil0[..]);
            mutual_inductance_2 +=
                crate::physics::linear_filament::inductance_piecewise_linear_filaments(
                    (
                        &xfil[0..ndiscr - 1],
                        &yfil[0..ndiscr - 1],
                        &zfil[0..ndiscr - 1],
                    ),
                    dlxyzfil0,
                    (&x[0..n - 1], &y[0..n - 1], &z[0..n - 1]),
                    dlxyzfil1,
                    false,
                ).unwrap();
        }

        println!("{mutual_inductance} {mutual_inductance_2}");
        assert!(approx(mutual_inductance_2, mutual_inductance, 1e-6, 1e-12));
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
