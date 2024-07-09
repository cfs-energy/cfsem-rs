//! Electromagnetics calculations.
pub mod biotsavart;
pub mod circular_filament;
pub mod gradshafranov;
pub mod linear_filament;

pub use circular_filament::{flux_circular_filament, flux_density_circular_filament};
