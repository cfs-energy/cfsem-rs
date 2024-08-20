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
