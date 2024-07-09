//! Magnetics calculations for circular current filaments.
use crate::math::{ellipe, ellipk};

use crate::{MU0_OVER_4PI, MU_0};

/// Flux contributions from some circular filaments to some observation points, which happens to be
/// the Green's function for the Grad-Shafranov elliptic operator, $\Delta^{\*}$.
///
/// # Arguments
///
/// * `current`: (A) current in each filament, length `m`
/// * `rfil`:    (m) r-coord of each filament, length `m`
/// * `zfil`:    (m) z-coord of each filament, length `m`
/// * `rprime`:  (m) r-coord of each observation point, length `n`
/// * `zprime`:  (m) z-coord of each observation point, length `n`
/// * `out_psi`:     (Wb) or (T-m^2) or (V-s), poloidal flux at observation locations, length `n`
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
pub fn flux_circular_filament(
    current: &[f64],
    rfil: &[f64],
    zfil: &[f64],
    rprime: &[f64],
    zprime: &[f64],
    out_psi: &mut [f64],
) -> Result<(), &'static str> {
    // Check lengths; Error if they do not match
    let m: usize = current.len();
    let n: usize = rprime.len();
    if rfil.len() != m || zfil.len() != m || zprime.len() != n || out_psi.len() != n {
        return Err("Length mismatch");
    }

    for i in 0..n {
        for j in 0..m {
            let rrprime = rfil[j] * rprime[i];
            let r_plus_rprime = rfil[j] + rprime[i];
            let z_minus_zprime = zfil[j] - zprime[i];
            let k2 = 4.0 * rrprime / (r_plus_rprime.powi(2) + z_minus_zprime.powi(2));

            out_psi[i] += current[j]
                * MU_0
                * (rrprime / k2).sqrt()
                * ((2.0 - k2) * ellipk(k2) - 2.0 * ellipe(k2));
        }
    }

    Ok(())
}

/// Off-axis Br,Bz components for a circular current filament in vacuum.
///
/// # Arguments
///
/// * `current`: (A) current in each filament, length `m`
/// * `rfil`:    (m) r-coord of each filament, length `m`
/// * `zfil`:    (m) z-coord of each filament, length `m`
/// * `rprime`:  (m) r-coord of each observation point, length `n`
/// * `zprime`:  (m) z-coord of each observation point, length `n`
/// * `out_br`:  (T), r-component of magnetic flux density at observation locations, length `n`
/// * `out_bz`:  (T), z-component of magnetic flux density at observation locations, length `n`
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
pub fn flux_density_circular_filament(
    current: &[f64],
    rfil: &[f64],
    zfil: &[f64],
    rprime: &[f64],
    zprime: &[f64],
    out_br: &mut [f64],
    out_bz: &mut [f64],
) -> Result<(), &'static str> {
    let n = current.len();
    let m = rprime.len();

    // Check lengths; Error if they do not match
    if rfil.len() != n
        || zfil.len() != n
        || zprime.len() != m
        || out_br.len() != m
        || out_bz.len() != m
    {
        return Err("Length mismatch");
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

            let a0 = 2.0 * current[i] / q.sqrt(); // [A/m]

            let f = ellipk(k2); // [nondim]
            let s = ellipe(k2) / (1.0 - k2); // [nondim]

            // Bake some reusable values
            let s_over_q = s / q; // [m^-2]
            let rfil2 = rfil[i] * rfil[i]; // [m^2]

            // Magnetic field intensity, less the factor of 4pi that we have adjusted out of mu_0
            let hr = (z / rprime[j]) * a0 * s_over_q.mul_add(rfil2 + r2 + z2, -f);
            let hz = a0 * s_over_q.mul_add(rfil2 - r2 - z2, f);

            // Magnetic flux density assuming vacuum permeability
            out_br[j] = MU0_OVER_4PI.mul_add(hr, out_br[j]);
            out_bz[j] = MU0_OVER_4PI.mul_add(hz, out_bz[j]);
        }
    }

    Ok(())
}
