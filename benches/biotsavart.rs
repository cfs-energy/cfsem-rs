#![allow(clippy::all)] // Clippy will attempt to remove black_box() internals

use cfsem::physics::biotsavart::{flux_density_biot_savart, flux_density_biot_savart_par};
use criterion::*;
use std::time::Duration;

fn bench_biotsavart(c: &mut Criterion) {
    let mut group = c.benchmark_group("Biot-Savart");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(5));

    // Examine logspace with fixed total throughput
    for nfac in [1, 10, 100, 1000].iter() {
        for nfils in (0_usize..=5).map(|i| 10_usize.pow(i as u32)) {
            // Filament inputs
            let nfils = nfils * nfac;
            let xfil = vec![1.0 / 7.0_f64; nfils];
            let yfil = vec![1.0 / 9.0_f64; nfils];
            let zfil = vec![1.0 / 11.0_f64; nfils];
            let dlxfil = vec![1.0 / 1.3_f64; nfils];
            let dlyfil = vec![1.0 / 2.3_f64; nfils];
            let dlzfil = vec![1.0 / 3.3_f64; nfils];
            let ifil = vec![0.5_f64; nfils];

            // Observation points
            let nobs = 1000;
            let nobs = nobs / nfac;
            let xobs = vec![-1.0 / 7.0_f64; nobs];
            let yobs = vec![-1.0 / 9.0_f64; nobs];
            let zobs = vec![-1.0 / 11.0_f64; nobs];

            let ntot = nobs * nfils;
            group.throughput(Throughput::Elements(ntot as u64));
            group.bench_with_input(
                BenchmarkId::new(
                    format!("Biot-Savart w/ Alloc\n{} Obs. Point(s)", nobs),
                    ntot,
                ),
                &ntot,
                |b, &_| {
                    b.iter(|| {
                        let n = xobs.len();
                        let (mut bx, mut by, mut bz) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
                        black_box(
                            flux_density_biot_savart(
                                (&xobs[..], &yobs[..], &zobs[..]),
                                (&xfil[..], &yfil[..], &zfil[..]),
                                (&dlxfil[..], &dlyfil[..], &dlzfil[..]),
                                &ifil[..],
                                (&mut bx, &mut by, &mut bz),
                            )
                            .unwrap(),
                        )
                    });
                },
            );
            group.bench_with_input(
                BenchmarkId::new(
                    format!("Biot-Savart Parallel w/ Alloc\n{} Obs. Point(s)", nobs),
                    ntot,
                ),
                &ntot,
                |b, &_| {
                    b.iter(|| {
                        let n = xobs.len();
                        let (mut bx, mut by, mut bz) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
                        black_box(
                            flux_density_biot_savart_par(
                                (&xobs[..], &yobs[..], &zobs[..]),
                                (&xfil[..], &yfil[..], &zfil[..]),
                                (&dlxfil[..], &dlyfil[..], &dlzfil[..]),
                                &ifil[..],
                                (&mut bx, &mut by, &mut bz),
                            )
                            .unwrap(),
                        )
                    });
                },
            );
        }
    }

    group.finish();
}

criterion_group!(group_bench_biotsavart, bench_biotsavart);
criterion_main!(group_bench_biotsavart);
