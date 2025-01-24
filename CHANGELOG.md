# Changelog

## 2.0.0 2025-01-21

### Added

* Add libm dep for reproducible trig functions
    * Rust std/core defers to libc or other platform-dependent math libraries for trig functions, which can cause platform-dependent results
    * libm is a pure rust implementation of most functions from libc, and is platform-independent to the extent that the processor's implementation of floating point math conforms to IEEE-754
* Add scalar calculations of vector potential, magnetic field, and poloidal flux of a circular filament extracted from vector loop
    * Scalar calculations now used as the inner function in the vector loops
    * This produces no performance regression, and some improvements in a few cases (up to a 4x speedup for some edge cases in the parallel variants)
* Add `mutual_inductance_circular_to_linear` family of functions for calculating mutual inductance between circular filaments and piecewise-linear paths
* Add `flux_density_circular_filament_cartesian` family of functions for calculating B-field from circular filaments to points in cartesian coordinates
* Add `cartesian_to_cylindrical` and `cylindrical_to_cartesian` conversion functions
* Add `decompose_filament` function for converting the start and end points of a filament to the midpoint and length vector
* Add `body_force_density_linear_filament` family of functions for calculating JxB force density due to a linear filament

### Changed

* Use `Slice::fill(0.0)` instead of manually zeroing output arrays
* !Consolidate function signatures of `flux_density_circular_filament[_par]`, `flux_circular_filament[_par]`, `vector_potential_circular_filament[_par]`, `flux_density_linear_filament[_par]`, `vector_potential_linear_filament[_par]`,
and `inductance_piecewise_linear_filaments`
    * For linear filament methods, this change formalizes the filament input as a point series describing piecewise-continuous segments with the same value of current, which is a functionality-breaking change
* !Remove deprecated `biot_savart` module, which has been superceded by the `linear_filament` module

# Changelog

## 1.1.0 2024-08-20

### Added

* Add vector potential calcs for linear and circular filaments w/ parallel variants
* Add parallel variants of circular filament flux and flux density calcs
* Add tests of serial and parallel variants to make sure they produce the same result
* Add tests of equivalence between flux/inductance, flux density, and vector potential calcs

### Changed

* Move Biot-Savart calcs to linear filament module and rename appropriately
  * Leave use-as references to prevent breaking change to API
* Eliminate small allocations from parallel variant of Biot-Savart to reduce overhead when running with a large number of cores
  * 40%-100% speedup for small numbers of observation points
* Defensively zero-out output slices
* Convert some `#[inline(always)]` directives to plain `#[inline]`

## 1.0.0 2024-07-09

### Added

* Initial release
