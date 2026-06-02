//! Energy maps for DNA and RNA, assembled from the generated `data` tables.

use std::collections::HashMap;
use std::sync::OnceLock;

use super::data;

/// Enthalpy/entropy parameters for one nucleic-acid type.
///
/// Mirrors the Python `seqfold.types.Energies` container.
pub struct Energies {
    pub bulge_loops: HashMap<i64, (f64, f64)>,
    pub complement: HashMap<u8, u8>,
    pub de: HashMap<String, (f64, f64)>,
    pub hairpin_loops: HashMap<i64, (f64, f64)>,
    pub multibranch: (f64, f64, f64, f64),
    pub internal_loops: HashMap<i64, (f64, f64)>,
    pub internal_mm: HashMap<String, (f64, f64)>,
    pub nn: HashMap<String, (f64, f64)>,
    pub terminal_mm: HashMap<String, (f64, f64)>,
    pub tri_tetra_loops: Option<HashMap<String, (f64, f64)>>,
}

fn bp(pairs: &[(&str, (f64, f64))]) -> HashMap<String, (f64, f64)> {
    pairs.iter().map(|(k, v)| (k.to_string(), *v)).collect()
}

fn loops(pairs: &[(i64, (f64, f64))]) -> HashMap<i64, (f64, f64)> {
    pairs.iter().copied().collect()
}

fn comp(pairs: &[(u8, u8)]) -> HashMap<u8, u8> {
    pairs.iter().copied().collect()
}

pub fn dna() -> &'static Energies {
    static DNA: OnceLock<Energies> = OnceLock::new();
    DNA.get_or_init(|| Energies {
        bulge_loops: loops(data::DNA_BULGE_LOOPS),
        complement: comp(data::DNA_COMPLEMENT),
        de: bp(data::DNA_DE),
        hairpin_loops: loops(data::DNA_HAIRPIN_LOOPS),
        multibranch: data::DNA_MULTIBRANCH,
        internal_loops: loops(data::DNA_INTERNAL_LOOPS),
        internal_mm: bp(data::DNA_INTERNAL_MM),
        nn: bp(data::DNA_NN),
        terminal_mm: bp(data::DNA_TERMINAL_MM),
        tri_tetra_loops: Some(bp(data::DNA_TRI_TETRA_LOOPS)),
    })
}

pub fn rna() -> &'static Energies {
    static RNA: OnceLock<Energies> = OnceLock::new();
    RNA.get_or_init(|| Energies {
        bulge_loops: loops(data::RNA_BULGE_LOOPS),
        complement: comp(data::RNA_COMPLEMENT),
        de: bp(data::RNA_DE),
        hairpin_loops: loops(data::RNA_HAIRPIN_LOOPS),
        multibranch: data::RNA_MULTIBRANCH,
        internal_loops: loops(data::RNA_INTERNAL_LOOPS),
        internal_mm: bp(data::RNA_INTERNAL_MM),
        nn: bp(data::RNA_NN),
        terminal_mm: bp(data::RNA_TERMINAL_MM),
        tri_tetra_loops: None,
    })
}
