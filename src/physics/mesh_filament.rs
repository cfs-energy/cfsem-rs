//! Methods for performing linear filament calcs on inputs
//! defined in mesh edge list format.
use rayon::{
    iter::{IndexedParallelIterator, ParallelIterator},
    slice::ParallelSliceMut,
};

use num_traits::Float;

use crate::{math::rss3, mesh::MeshSegmentList, MU0_OVER_4PI};
use crate::{
    math::{decompose_filament, dot3},
    physics::linear_filament::vector_potential_linear_filament_scalar,
};

/// Convert a point to f64 values
fn convert_point<T>(p: (T, T, T)) -> (f64, f64, f64)
where
    T: Into<f64>,
{
    (p.0.into(), p.1.into(), p.2.into())
}

/// Mutual inductance from each edge in mesh 1 to each edge in mesh 2.
/// If mesh 2 is not populated, the self-inductance of mesh 1 is taken,
/// using the thin-filament scalar self-inductance for self-terms.
pub fn mesh_inductance<T>(m1: &MeshSegmentList<T>, m2: Option<&MeshSegmentList<T>>) -> Vec<T>
where
    T: Float + Into<f64> + From<f64> + Send + Sync,
{
    // If there is no second mesh, we're doing self inductance
    let self_inductance = m2.is_none();
    let m2 = m2.unwrap_or(m1);

    // Allocate for the output
    let n1 = m1.edges().len();
    let n2 = m2.edges().len();
    let n_out = n1 * n2;
    let mut out = vec![T::zero(); n_out];

    // Chunk output,
    // taking each chunk as all the contributions of an edge in mesh 1
    // to each edge in mesh 2
    let outc = out.par_chunks_mut(n1);

    // Loop over pairs of edges, taking M = dot((A/I), dL) .
    outc.enumerate().for_each(|(i, o)| {
        let (first, second) = m1.edges()[i];
        let xyzifil1 = (
            convert_point(m1.nodes()[first]),
            convert_point(m1.nodes()[second]),
            1.0,
        );
        for j in 0..m2.edges().len() {
            // Handle segment self-inductance case
            if self_inductance && i == j {
                let (start, end, _) = xyzifil1;
                let length = rss3(end.0 - start.0, end.1 - start.1, end.2 - start.2);
                o[j] = (0.5 * MU0_OVER_4PI * length).into();
                continue;
            }

            // Handle mutual-inductance case
            let edge2 = m2.edges()[j];
            let (start2, end2) = (
                convert_point(m2.nodes()[edge2.0]),
                convert_point(m2.nodes()[edge2.1]),
            );

            //    First, get vector potential from edge 1 (with its midpoint as the source) to the midpoint of edge 2
            //    with unit current in edge 1 in order to extract A/I.
            let (midpoint2, dl2) = decompose_filament(start2, end2);
            let (ax_per_amp, ay_per_amp, az_per_amp) =
                vector_potential_linear_filament_scalar(xyzifil1, midpoint2);

            //    Take M = dot((A/I), dL)
            let m = dot3(ax_per_amp, ay_per_amp, az_per_amp, dl2.0, dl2.1, dl2.2);
            o[j] = m.into();
        }
    });

    out
}
