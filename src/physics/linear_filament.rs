//! Magnetics calculations for linear current filaments.
use crate::math::{dot3, rss3};
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
