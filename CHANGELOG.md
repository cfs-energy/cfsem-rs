# Changelog

## 1.2.0 2025-01-21

### Added

* Add scalar calculations of vector potential, magnetic field, and poloidal flux of a circular filament extracted from vector loop
    * Scalar calculations now used as the inner function in the vector loops
    * This produces no performance regression, and some improvements in a few cases (up to a 4x speedup for some edge cases in the parallel variants)
* Add `mutual_inductance_circular_to_linear` family of functions for calculating mutual inductance between circular filaments and piecewise-linear paths
* Add `cartesian_to_cylindrical` and `cylindrical_to_cartesian` conversion functions


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
