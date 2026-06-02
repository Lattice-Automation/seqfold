//! seqfold — predict the minimum free energy structure of nucleic acids.
//!
//! The numerical engine lives in [`core`] (pure Rust, no Python). When built
//! with the `python` feature, [`python`] exposes it as the `seqfold._core`
//! extension module via PyO3.

pub mod core;

#[cfg(feature = "python")]
mod python;
