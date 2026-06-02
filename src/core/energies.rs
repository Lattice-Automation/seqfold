//! Energy maps for DNA and RNA, assembled from the generated `data` tables.

use std::collections::HashMap;
use std::sync::OnceLock;

use super::data;

/// Maps a nucleotide byte to a 3-bit symbol code; 7 marks "invalid".
const fn build_sym() -> [u8; 256] {
    let mut t = [7u8; 256];
    t[b'A' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b'U' as usize] = 4;
    t[b'N' as usize] = 5;
    t[b'.' as usize] = 6;
    t
}

static SYM: [u8; 256] = build_sym();

/// Number of slots in a dense pair table: 4 symbols x 3 bits = 12 bits.
const TABLE_LEN: usize = 1 << 12;

/// Pack the four nucleotides of a stack/pair key ("b0 b1 / b2 b3") into a
/// 12-bit index. `.` (a dangling end) is a valid symbol.
#[inline(always)]
pub fn encode4(b0: u8, b1: u8, b2: u8, b3: u8) -> usize {
    ((SYM[b0 as usize] as usize) << 9)
        | ((SYM[b1 as usize] as usize) << 6)
        | ((SYM[b2 as usize] as usize) << 3)
        | (SYM[b3 as usize] as usize)
}

/// A dense lookup table indexed by [`encode4`], used to avoid string
/// allocation + hashing on the hot folding path.
pub type PairTable = Vec<Option<(f64, f64)>>;

/// Enthalpy/entropy parameters for one nucleic-acid type.
///
/// Mirrors the Python `seqfold.types.Energies` container. The `*_t` fields and
/// `comp_lut` are dense-array equivalents of the string-keyed maps, used by the
/// folding hot path; the string maps are retained for the (cold) Tm path and
/// for the variable-length tri/tetra-loop keys.
pub struct Energies {
    pub bulge_loops: HashMap<i64, (f64, f64)>,
    pub complement: HashMap<u8, u8>,
    pub comp_lut: [u8; 256],
    pub de: HashMap<String, (f64, f64)>,
    pub de_t: PairTable,
    pub hairpin_loops: HashMap<i64, (f64, f64)>,
    pub multibranch: (f64, f64, f64, f64),
    pub internal_loops: HashMap<i64, (f64, f64)>,
    pub internal_mm: HashMap<String, (f64, f64)>,
    pub internal_mm_t: PairTable,
    pub nn: HashMap<String, (f64, f64)>,
    pub nn_t: PairTable,
    pub terminal_mm: HashMap<String, (f64, f64)>,
    pub terminal_mm_t: PairTable,
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

fn comp_lut(pairs: &[(u8, u8)]) -> [u8; 256] {
    let mut t = [0u8; 256];
    for &(k, v) in pairs {
        t[k as usize] = v;
    }
    t
}

/// Build a dense pair table from string-keyed data, keeping only entries whose
/// key is a 5-char pair ("XY/ZW"); non-pair keys (e.g. "init") are skipped.
fn dense(pairs: &[(&str, (f64, f64))]) -> PairTable {
    let mut t: PairTable = vec![None; TABLE_LEN];
    for (k, v) in pairs {
        let kb = k.as_bytes();
        if kb.len() == 5 && kb[2] == b'/' {
            t[encode4(kb[0], kb[1], kb[3], kb[4])] = Some(*v);
        }
    }
    t
}

pub fn dna() -> &'static Energies {
    static DNA: OnceLock<Energies> = OnceLock::new();
    DNA.get_or_init(|| Energies {
        bulge_loops: loops(data::DNA_BULGE_LOOPS),
        complement: comp(data::DNA_COMPLEMENT),
        comp_lut: comp_lut(data::DNA_COMPLEMENT),
        de: bp(data::DNA_DE),
        de_t: dense(data::DNA_DE),
        hairpin_loops: loops(data::DNA_HAIRPIN_LOOPS),
        multibranch: data::DNA_MULTIBRANCH,
        internal_loops: loops(data::DNA_INTERNAL_LOOPS),
        internal_mm: bp(data::DNA_INTERNAL_MM),
        internal_mm_t: dense(data::DNA_INTERNAL_MM),
        nn: bp(data::DNA_NN),
        nn_t: dense(data::DNA_NN),
        terminal_mm: bp(data::DNA_TERMINAL_MM),
        terminal_mm_t: dense(data::DNA_TERMINAL_MM),
        tri_tetra_loops: Some(bp(data::DNA_TRI_TETRA_LOOPS)),
    })
}

pub fn rna() -> &'static Energies {
    static RNA: OnceLock<Energies> = OnceLock::new();
    RNA.get_or_init(|| Energies {
        bulge_loops: loops(data::RNA_BULGE_LOOPS),
        complement: comp(data::RNA_COMPLEMENT),
        comp_lut: comp_lut(data::RNA_COMPLEMENT),
        de: bp(data::RNA_DE),
        de_t: dense(data::RNA_DE),
        hairpin_loops: loops(data::RNA_HAIRPIN_LOOPS),
        multibranch: data::RNA_MULTIBRANCH,
        internal_loops: loops(data::RNA_INTERNAL_LOOPS),
        internal_mm: bp(data::RNA_INTERNAL_MM),
        internal_mm_t: dense(data::RNA_INTERNAL_MM),
        nn: bp(data::RNA_NN),
        nn_t: dense(data::RNA_NN),
        terminal_mm: bp(data::RNA_TERMINAL_MM),
        terminal_mm_t: dense(data::RNA_TERMINAL_MM),
        tri_tetra_loops: None,
    })
}
