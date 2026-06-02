//! Nucleic-acid secondary-structure prediction (Zuker & Stiegler, 1981).
//!
//! A faithful port of the original `seqfold/fold.py`. Indices are kept as
//! `i64` because several helpers are intentionally called with `-1` (a dangling
//! end) or with indices one past the sequence end; sequence access converts to
//! `usize` only where the index is known to be valid.

use std::collections::HashSet;

use rayon::prelude::*;
use smallvec::{smallvec, SmallVec};

use super::energies::{self, Energies};
use super::pyfloat::pyround;

const INF: f64 = f64::INFINITY;
const NEG_INF: f64 = f64::NEG_INFINITY;

/// The basepairs of a structure, stored in 32-bit (sequences are short) and
/// inline for the common 0-1 element case. Multi-branch structures (>1 pair)
/// spill to the heap.
type Ij = SmallVec<[(i32, i32); 1]>;

/// A working list of branches inside `mb`; inline for up to 4 branches (the
/// typical multi-branch fan-out), arithmetic kept in `i64`.
type Branches = SmallVec<[(i64, i64); 4]>;

/// Sequences at least this long fill each anti-diagonal in parallel; shorter
/// ones run sequentially (thread overhead outweighs the work).
const PARALLEL_THRESHOLD: usize = 64;

/// A compact descriptor of a structure, stored in the DP caches in place of a
/// formatted label. The human-readable `desc` string is rendered only for the
/// handful of structures that end up in the traceback (see [`render_desc`]),
/// which avoids O(n^3) `format!` allocations (and the cross-thread allocator
/// contention they cause) during the fill.
#[derive(Clone, Debug, PartialEq)]
enum Desc {
    /// No label (sentinels, the isolated-pair penalty).
    None,
    /// A hairpin closed at (i, j).
    Hairpin,
    /// A stack / bulge / interior loop; the exact label is reconstructed from
    /// (i, j) and the inner pair (ij[0]) at render time.
    Pair,
    /// A multi-branch structure; carries the counts the label needs.
    Bifurcation { unpaired: i32, count: u32 },
}

/// A single structure with a free energy, description, and inward children.
#[derive(Clone, Debug)]
pub struct Struct {
    pub e: f64,
    /// Rendered label. Empty for cache cells; filled in for traceback output.
    pub desc: String,
    pub ij: Ij,
    /// Compact descriptor used during the fill (cache cells only).
    tag: Desc,
}

/// Python's `Struct.__eq__` compares only `e` and `ij` (not `desc`).
impl PartialEq for Struct {
    fn eq(&self, other: &Self) -> bool {
        self.e == other.e && self.ij == other.ij
    }
}

impl Struct {
    /// Construct an output structure with an already-rendered label.
    pub fn new(e: f64, desc: &str, ij: Ij) -> Struct {
        Struct {
            e,
            desc: desc.to_string(),
            ij,
            tag: Desc::None,
        }
    }

    /// Construct a cache cell carrying a compact tag (no label allocation).
    fn tagged(e: f64, tag: Desc, ij: Ij) -> Struct {
        Struct {
            e,
            desc: String::new(),
            ij,
            tag,
        }
    }

    /// The `STRUCT_DEFAULT` sentinel: `Struct(-inf)` with an empty `ij`.
    pub fn default_() -> Struct {
        Struct::tagged(NEG_INF, Desc::None, SmallVec::new())
    }

    /// The `STRUCT_NULL` sentinel: `Struct(inf)` with an empty `ij`.
    pub fn null() -> Struct {
        Struct::tagged(INF, Desc::None, SmallVec::new())
    }

    /// Python `__bool__`: false when the energy is +/- infinity.
    pub fn truthy(&self) -> bool {
        self.e != INF && self.e != NEG_INF
    }
}

/// Reconstruct the human-readable label for a structure at cell (i, j).
fn render_desc(tag: &Desc, s: &[u8], i: usize, j: usize, inner: Option<(i32, i32)>, n: usize) -> String {
    match tag {
        Desc::None => String::new(),
        Desc::Hairpin => format!(
            "HAIRPIN:{}",
            pair_str(s, i as i64, i as i64 + 1, j as i64, j as i64 - 1)
        ),
        Desc::Bifurcation { unpaired, count } => {
            format!("BIFURCATION:{}n/{}h", unpaired, count)
        }
        Desc::Pair => {
            let (i1i, j1i) = inner.expect("pair desc needs an inner pair");
            let (i1, j1) = (i1i as usize, j1i as usize);
            if i1 == i + 1 && j1 == j - 1 {
                // stacking pair (possibly with a dangling end)
                let pair = pair_str(s, i as i64, i1i as i64, j as i64, j1i as i64);
                if (i > 0 && j == n - 1) || (i == 0 && j < n - 1) {
                    format!("STACK_DE:{}", pair)
                } else {
                    format!("STACK:{}", pair)
                }
            } else if i1 > i + 1 && j1 == j - 1 {
                format!("BULGE:{}", i1 - i)
            } else if i1 == i + 1 && j1 < j - 1 {
                format!("BULGE:{}", j - j1)
            } else if i1 - i == 2 && j - j1 == 2 {
                // technically an interior loop of 1; really a 1bp mismatch
                let loop_left = std::str::from_utf8(&s[i..=i1]).unwrap();
                let mut loop_right_rev: Vec<u8> = s[j1..=j].to_vec();
                loop_right_rev.reverse();
                let loop_right_rev = String::from_utf8(loop_right_rev).unwrap();
                format!("STACK:{}/{}", loop_left, loop_right_rev)
            } else {
                format!("INTERIOR_LOOP:{}/{}", i1 - i, j - j1)
            }
        }
    }
}

/// Error raised for invalid sequences (maps to Python `RuntimeError`).
#[derive(Debug)]
pub struct FoldError(pub String);

type Cache = Vec<Vec<Struct>>;

/// The character at `i`, or `.` for a negative (dangling) index.
fn ch(s: &[u8], i: i64) -> char {
    if i < 0 {
        '.'
    } else {
        s[i as usize] as char
    }
}

#[inline(always)]
fn comp(emap: &Energies, b: u8) -> u8 {
    emap.comp_lut[b as usize]
}

/// Pack the four (possibly dangling) indices into a dense-table code.
#[inline(always)]
fn code4(s: &[u8], i: i64, i1: i64, j: i64, j1: i64) -> usize {
    let g = |x: i64| if x < 0 { b'.' } else { s[x as usize] };
    energies::encode4(g(i), g(i1), g(j), g(j1))
}

/// Fold the sequence and return the list of minimum-free-energy structures.
pub fn fold(seq: &str, temp: f64) -> Result<Vec<Struct>, FoldError> {
    let (v_cache, w_cache) = cache(seq, temp)?;
    // traceback renders labels from the (upper-cased) sequence, matching cache()
    let upper = seq.to_uppercase();
    let s = upper.as_bytes();
    Ok(traceback(s, 0, s.len() - 1, &v_cache, &w_cache))
}

/// Fold and return just the rounded delta-G of the structure.
pub fn dg(seq: &str, temp: f64) -> Result<f64, FoldError> {
    let structs = fold(seq, temp)?;
    let dg_sum: f64 = structs.iter().map(|s| s.e).sum();
    Ok(pyround(dg_sum, 2))
}

/// Fold and return the (i, j) -> min-free-energy matrix (the W cache energies).
pub fn dg_cache(seq: &str, temp: f64) -> Result<Vec<Vec<f64>>, FoldError> {
    let (_v, w_cache) = cache(seq, temp)?;
    Ok(w_cache
        .iter()
        .map(|row| row.iter().map(|s| s.e).collect())
        .collect())
}

/// Dot-bracket notation for a folded structure.
pub fn dot_bracket(seq: &str, structs: &[Struct]) -> String {
    let mut result = vec![b'.'; seq.len()];
    for s in structs {
        if s.ij.len() == 1 {
            let (i, j) = s.ij[0];
            result[i as usize] = b'(';
            result[j as usize] = b')';
        }
    }
    String::from_utf8(result).expect("ascii")
}

/// Build the V and W caches for a sequence.
pub fn cache(seq: &str, temp: f64) -> Result<(Cache, Cache), FoldError> {
    let seq = seq.to_uppercase();
    let temp = temp + 273.15; // kelvin
    let s = seq.as_bytes();

    // figure out whether it's DNA or RNA, choose energy map
    let bps: HashSet<u8> = s.iter().copied().collect();
    if bps.contains(&b'U') && bps.contains(&b'T') {
        return Err(FoldError(
            "Both T and U in sequence. Provide one or the other for DNA OR RNA.".to_string(),
        ));
    }

    let mut dna = true;
    if bps.iter().all(|b| b"AUCG".contains(b)) {
        dna = false;
    } else if bps.iter().any(|b| !b"ATGC".contains(b)) {
        let diff: Vec<char> = bps
            .iter()
            .filter(|b| !b"ATUGC".contains(b))
            .map(|&b| b as char)
            .collect();
        return Err(FoldError(format!(
            "Unknown bp: {:?}. Only DNA/RNA foldable",
            diff
        )));
    }
    let emap: &Energies = if dna { energies::dna() } else { energies::rna() };

    Ok(fill(s, temp, emap))
}

/// Fill the V and W caches bottom-up by anti-diagonal (span `d = j - i`).
///
/// Every cell on a given anti-diagonal depends only on cells of strictly
/// smaller span, so each diagonal can be computed in parallel. A cell's value
/// is a pure function of already-finalized smaller cells, so the result is
/// independent of evaluation order — identical to the recursive version.
fn fill(s: &[u8], temp: f64, emap: &Energies) -> (Cache, Cache) {
    let n = s.len();
    let mut v_cache: Cache = vec![vec![Struct::default_(); n]; n];
    let mut w_cache: Cache = vec![vec![Struct::default_(); n]; n];

    // W over spans < 4 is NULL (the recursive base case).
    for d in 0..4 {
        if d >= n {
            break;
        }
        for i in 0..(n - d) {
            w_cache[i][i + d] = Struct::null();
        }
    }

    let parallel = n >= PARALLEL_THRESHOLD;
    for d in 4..n {
        let row: Vec<(Struct, Struct)> = if parallel {
            (0..(n - d))
                .into_par_iter()
                .map(|i| compute_cell(s, i, i + d, temp, &v_cache, &w_cache, emap))
                .collect()
        } else {
            (0..(n - d))
                .map(|i| compute_cell(s, i, i + d, temp, &v_cache, &w_cache, emap))
                .collect()
        };
        for (i, (vs, ws)) in row.into_iter().enumerate() {
            v_cache[i][i + d] = vs;
            w_cache[i][i + d] = ws;
        }
    }

    (v_cache, w_cache)
}

/// Compute the V and W structs for one cell, reading only finalized cells.
#[inline]
fn compute_cell(
    s: &[u8],
    i: usize,
    j: usize,
    temp: f64,
    v_cache: &Cache,
    w_cache: &Cache,
    emap: &Energies,
) -> (Struct, Struct) {
    let vs = compute_v(s, i, j, temp, v_cache, w_cache, emap);
    let ws = compute_w(s, i, j, temp, &vs, w_cache, emap);
    (vs, ws)
}

/// Lowest free energy structure in the Sij subsequence (Fig. 2B), given that
/// all smaller-span cells are already in the caches.
fn compute_w(
    s: &[u8],
    i: usize,
    j: usize,
    temp: f64,
    v_ij: &Struct,
    w_cache: &Cache,
    emap: &Energies,
) -> Struct {
    let w1 = w_cache[i + 1][j].clone();
    let w2 = w_cache[i][j - 1].clone();
    let w3 = v_ij.clone();

    let mut w4 = Struct::null();
    for k in (i + 1)..(j - 1) {
        let w4_test = mb(s, i, k, j, temp, w_cache, emap, false);
        if w4_test.truthy() && w4_test.e < w4.e {
            w4 = w4_test;
        }
    }

    min_struct(&[w1, w2, w3, w4])
}

/// Minimum free energy of the structure between i and j (Fig. 2B), given that
/// all smaller-span cells are already in the caches.
fn compute_v(
    s: &[u8],
    i: usize,
    j: usize,
    temp: f64,
    v_cache: &Cache,
    w_cache: &Cache,
    emap: &Energies,
) -> Struct {
    let n = s.len();

    // the ends must basepair for V(i,j)
    if comp(emap, s[i]) != s[j] {
        return Struct::null();
    }

    // if the basepair is isolated, and the seq large, penalize at 1,600 kcal/mol
    let mut isolated_outer = true;
    if i != 0 && j < n - 1 {
        isolated_outer = comp(emap, s[i - 1]) != s[j + 1];
    }
    let isolated_inner = comp(emap, s[i + 1]) != s[j - 1];

    if isolated_outer && isolated_inner {
        return Struct::tagged(1600.0, Desc::None, SmallVec::new());
    }

    // E1 = FH(i, j); hairpin. The label is reconstructed at traceback time.
    let e1 = Struct::tagged(hairpin(s, i, j, temp, emap), Desc::Hairpin, SmallVec::new());
    if j - i == 4 {
        return e1;
    }

    // E2 = min{FL(i, j, i', j') + V(i', j')}, i<i'<j'<j
    //
    // Only the energy and the winning inner pair are tracked; the exact label
    // (STACK / BULGE / INTERIOR_LOOP / ...) is reconstructed from those indices
    // at traceback time, so no strings are built during the search.
    let mut e2 = Struct::tagged(INF, Desc::None, SmallVec::new());
    let mut best: Option<(f64, usize, usize)> = None;
    if j >= 4 {
        let pair_inner_left =
            emap.nn_t[code4(s, i as i64, i as i64 + 1, j as i64, j as i64 - 1)].is_some();

        for i1 in (i + 1)..(j - 4) {
            for j1 in (i1 + 4)..j {
                if comp(emap, s[i1]) != s[j1] {
                    continue;
                }

                let stack = i1 == i + 1 && j1 == j - 1;
                let bulge_left = i1 > i + 1;
                let bulge_right = j1 < j - 1;

                let mut e2_test: f64;
                if stack {
                    e2_test = stack_e(s, i as i64, i1 as i64, j as i64, j1 as i64, temp, emap);
                } else if bulge_left && bulge_right && {
                    let pair_inner = pair_inner_left
                        || emap.nn_t[code4(s, i1 as i64 - 1, i1 as i64, j1 as i64 + 1, j1 as i64)]
                            .is_some();
                    !pair_inner
                } {
                    e2_test =
                        internal_loop(s, i as i64, i1 as i64, j as i64, j1 as i64, temp, emap);
                } else if bulge_left && !bulge_right {
                    e2_test = bulge(s, i as i64, i1 as i64, j as i64, j1 as i64, temp, emap);
                } else if !bulge_left && bulge_right {
                    e2_test = bulge(s, i as i64, i1 as i64, j as i64, j1 as i64, temp, emap);
                } else {
                    continue;
                }

                e2_test += v_cache[i1][j1].e;
                let cur = best.map_or(INF, |b| b.0);
                if e2_test != NEG_INF && e2_test < cur {
                    best = Some((e2_test, i1, j1));
                }
            }
        }

        if let Some((e2_test, i1, j1)) = best {
            e2 = Struct::tagged(e2_test, Desc::Pair, smallvec![(i1 as i32, j1 as i32)]);
        }
    }

    // E3 = min{W(i+1,i') + W(i'+1,j-1)}, i+1<i'<j-2
    let mut e3 = Struct::null();
    if !isolated_outer || i == 0 || j == n - 1 {
        for k in (i + 1)..(j - 1) {
            let e3_test = mb(s, i, k, j, temp, w_cache, emap, true);
            if e3_test.truthy() && e3_test.e < e3.e {
                e3 = e3_test;
            }
        }
    }

    min_struct(&[e1, e2, e3])
}

/// A stack representation string ("AT/TG"), the key for the NN maps.
fn pair_str(s: &[u8], i: i64, i1: i64, j: i64, j1: i64) -> String {
    format!("{}{}/{}{}", ch(s, i), ch(s, i1), ch(s, j), ch(s, j1))
}

/// The struct with the lowest free energy that isn't -inf (undefined).
fn min_struct(structs: &[Struct]) -> Struct {
    let mut s = Struct::null();
    for st in structs {
        if st.e != NEG_INF && st.e < s.e {
            s = st.clone();
        }
    }
    s
}

/// Free energy from delta-h, delta-s, temp.
fn d_g(d_h: f64, d_s: f64, temp: f64) -> f64 {
    d_h - temp * (d_s / 1000.0)
}

/// Jacobson-Stockmayer length extrapolation.
fn j_s(query_len: i64, known_len: i64, d_g_x: f64, temp: f64) -> f64 {
    let gas_constant = 1.9872e-3;
    d_g_x + 2.44 * gas_constant * temp * (query_len as f64 / known_len as f64).ln()
}

/// Free energy of a stack / dangling end / terminal mismatch.
fn stack_e(s: &[u8], i: i64, i1: i64, j: i64, j1: i64, temp: f64, emap: &Energies) -> f64 {
    let n = s.len() as i64;
    if i >= n || i1 >= n || j >= n || j1 >= n {
        return 0.0;
    }

    let p = code4(s, i, i1, j, j1);
    if i == -1 || i1 == -1 || j == -1 || j1 == -1 {
        // it's a dangling end
        let (d_h, d_s) = emap.de_t[p].expect("de");
        return d_g(d_h, d_s, temp);
    }

    if i > 0 && j < n - 1 {
        // it's internal
        let (d_h, d_s) = emap.nn_t[p].or_else(|| emap.internal_mm_t[p]).expect("nn/imm");
        return d_g(d_h, d_s, temp);
    }

    if i == 0 && j == n - 1 {
        // it's terminal
        let (d_h, d_s) = emap.nn_t[p].or_else(|| emap.terminal_mm_t[p]).expect("nn/tmm");
        return d_g(d_h, d_s, temp);
    }

    if i > 0 && j == n - 1 {
        // it's dangling on the left
        let (d_h, d_s) = emap.nn_t[p].or_else(|| emap.terminal_mm_t[p]).expect("nn/tmm");
        let mut dgv = d_g(d_h, d_s, temp);

        let pd = energies::encode4(s[(i - 1) as usize], s[i as usize], b'.', s[j as usize]);
        if let Some((d_h2, d_s2)) = emap.de_t[pd] {
            dgv += d_g(d_h2, d_s2, temp);
        }
        return dgv;
    }

    if i == 0 && j < n - 1 {
        // it's dangling on the right
        let (d_h, d_s) = emap.nn_t[p].or_else(|| emap.terminal_mm_t[p]).expect("nn/tmm");
        let mut dgv = d_g(d_h, d_s, temp);

        let pd = energies::encode4(b'.', s[i as usize], s[(j + 1) as usize], s[j as usize]);
        if let Some((d_h2, d_s2)) = emap.de_t[pd] {
            dgv += d_g(d_h2, d_s2, temp);
        }
        return dgv;
    }

    0.0
}

/// Free energy of a hairpin closed at (i, j).
fn hairpin(s: &[u8], i: usize, j: usize, temp: f64, emap: &Energies) -> f64 {
    if j - i < 4 {
        return INF;
    }

    let hp = &s[i..=j];
    let hairpin_len = (hp.len() - 2) as i64;
    let p = code4(s, i as i64, i as i64 + 1, j as i64, j as i64 - 1);

    if comp(emap, hp[0]) != hp[hp.len() - 1] {
        panic!("hairpin: no closing pair");
    }

    let mut d_gv = 0.0;
    if let Some(tt) = &emap.tri_tetra_loops {
        let key = std::str::from_utf8(hp).unwrap();
        if let Some(&(d_h, d_s)) = tt.get(key) {
            d_gv = d_g(d_h, d_s, temp);
        }
    }

    // size penalty
    if let Some(&(d_h, d_s)) = emap.hairpin_loops.get(&hairpin_len) {
        d_gv += d_g(d_h, d_s, temp);
    } else {
        let (d_h, d_s) = emap.hairpin_loops[&30];
        let d_g_inc = d_g(d_h, d_s, temp);
        d_gv += j_s(hairpin_len, 30, d_g_inc, temp);
    }

    // terminal mismatch
    if hairpin_len > 3 {
        if let Some((d_h, d_s)) = emap.terminal_mm_t[p] {
            d_gv += d_g(d_h, d_s, temp);
        }
    }

    // length-3 AT closing penalty
    if hairpin_len == 3 && (hp[0] == b'A' || hp[hp.len() - 1] == b'A') {
        d_gv += 0.5;
    }

    d_gv
}

/// Free energy of a bulge.
fn bulge(s: &[u8], i: i64, i1: i64, j: i64, j1: i64, temp: f64, emap: &Energies) -> f64 {
    let loop_len = (i1 - i - 1).max(j - j1 - 1);
    if loop_len <= 0 {
        panic!("bulge: non-positive loop");
    }

    let mut d_gv = if let Some(&(d_h, d_s)) = emap.bulge_loops.get(&loop_len) {
        d_g(d_h, d_s, temp)
    } else {
        let (d_h, d_s) = emap.bulge_loops[&30];
        let g = d_g(d_h, d_s, temp);
        j_s(loop_len, 30, g, temp)
    };

    if loop_len == 1 {
        debug_assert!(emap.nn_t[code4(s, i, i1, j, j1)].is_some());
        d_gv += stack_e(s, i, i1, j, j1, temp, emap);
    }

    if [i, i1, j, j1].iter().any(|&k| s[k as usize] == b'A') {
        d_gv += 0.5;
    }

    d_gv
}

/// Free energy of an internal loop.
fn internal_loop(s: &[u8], i: i64, i1: i64, j: i64, j1: i64, temp: f64, emap: &Energies) -> f64 {
    let loop_left = i1 - i - 1;
    let loop_right = j - j1 - 1;
    let loop_len = loop_left + loop_right;

    if loop_left < 1 || loop_right < 1 {
        panic!("internal_loop: bad loop");
    }

    if loop_left == 1 && loop_right == 1 {
        let mm_left = stack_e(s, i, i1, j, j1, temp, emap);
        let mm_right = stack_e(s, i1 - 1, i1, j1 + 1, j1, temp, emap);
        return mm_left + mm_right;
    }

    let mut d_gv = if let Some(&(d_h, d_s)) = emap.internal_loops.get(&loop_len) {
        d_g(d_h, d_s, temp)
    } else {
        let (d_h, d_s) = emap.internal_loops[&30];
        let g = d_g(d_h, d_s, temp);
        j_s(loop_len, 30, g, temp)
    };

    let loop_asymmetry = (loop_left - loop_right).abs();
    d_gv += 0.3 * loop_asymmetry as f64;

    let (d_h, d_s) = emap.terminal_mm_t[code4(s, i, i + 1, j, j - 1)].expect("tmm");
    d_gv += d_g(d_h, d_s, temp);

    let (d_h, d_s) = emap.terminal_mm_t[code4(s, i1 - 1, i1, j1 + 1, j1)].expect("tmm");
    d_gv += d_g(d_h, d_s, temp);

    d_gv
}

/// Recursively gather basepairing branches into `branches` (reads only the
/// already-finalized W cache).
fn add_branch(st: &Struct, branches: &mut Branches, w_cache: &Cache) {
    if !st.truthy() || st.ij.is_empty() {
        return;
    }
    if st.ij.len() == 1 {
        let (a, b) = st.ij[0];
        branches.push((a as i64, b as i64));
        return;
    }
    let ij = st.ij.clone();
    for (i1, j1) in ij {
        let sub = w_cache[i1 as usize][j1 as usize].clone();
        add_branch(&sub, branches, w_cache);
    }
}

/// Multi-branch energy penalty (linear; Jaeger, Turner & Zuker, 1989). Reads
/// only finalized (smaller-span) cells of the W cache.
fn mb(
    s: &[u8],
    i: usize,
    k: usize,
    j: usize,
    temp: f64,
    w_cache: &Cache,
    emap: &Energies,
    helix: bool,
) -> Struct {
    let (left, right) = if helix {
        (
            w_cache[i + 1][k].clone(),
            w_cache[k + 1][j - 1].clone(),
        )
    } else {
        (w_cache[i][k].clone(), w_cache[k + 1][j].clone())
    };

    if !left.truthy() || !right.truthy() {
        return Struct::null();
    }

    let mut branches: Branches = SmallVec::new();
    add_branch(&left, &mut branches, w_cache);
    add_branch(&right, &mut branches, w_cache);

    // this isn't multi-branched
    if branches.len() < 2 {
        return Struct::null();
    }

    let (ii, jj) = (i as i64, j as i64);
    if helix {
        branches.push((ii, jj));
    }

    let branches_count = branches.len();
    let blen = branches.len();
    let mut unpaired: i64 = 0;
    let mut e_sum = 0.0;
    for index in 0..blen {
        let (i2, j2) = branches[index];
        let (_, j1) = branches[(index + blen - 1) % blen];
        let (i3, j3) = branches[(index + 1) % blen];

        let mut unpaired_left = 0i64;
        let mut unpaired_right = 0i64;
        let mut de = 0.0;
        if index == blen - 1 && !helix {
            // pass
        } else if (i3, j3) == (ii, jj) {
            unpaired_left = i2 - j1 - 1;
            unpaired_right = j3 - j2 - 1;

            if unpaired_left != 0 && unpaired_right != 0 {
                de = stack_e(s, i2 - 1, i2, j2 + 1, j2, temp, emap);
            } else if unpaired_right != 0 {
                de = stack_e(s, -1, i2, j2 + 1, j2, temp, emap);
                if unpaired_right == 1 {
                    de = stack_e(s, i3, -1, j3, j3 - 1, temp, emap).min(de);
                }
            }
        } else if (i2, j2) == (ii, jj) {
            unpaired_left = j2 - j1 - 1;
            unpaired_right = i3 - i2 - 1;

            if unpaired_left != 0 && unpaired_right != 0 {
                de = stack_e(s, i2 - 1, i2, j2 + 1, j2, temp, emap);
            } else if unpaired_right != 0 {
                de = stack_e(s, i2, i2 + 1, j2, -1, temp, emap);
                if unpaired_right == 1 {
                    de = stack_e(s, i3 - 1, i3, -1, j3, temp, emap).min(de);
                }
            }
        } else {
            unpaired_left = i2 - j1 - 1;
            unpaired_right = i3 - j2 - 1;

            if unpaired_left != 0 && unpaired_right != 0 {
                de = stack_e(s, i2 - 1, i2, j2 + 1, j2, temp, emap);
            } else if unpaired_right != 0 {
                de = stack_e(s, -1, i2, j2 + 1, j2, temp, emap);
                if unpaired_right == 1 {
                    de = stack_e(s, i2 - 1, i2, j2 + 1, j2, temp, emap).min(de);
                }
            }
        }
        let _ = unpaired_left;

        e_sum += de;
        unpaired += unpaired_right;
        assert!(unpaired_right >= 0);

        if (i2, j2) != (ii, jj) {
            e_sum += w_cache[i2 as usize][j2 as usize].e;
        }
    }

    assert!(unpaired >= 0);

    let (a, b, c, d) = emap.multibranch;
    let mut e_multibranch = a + b * branches.len() as f64 + c * unpaired as f64;
    if unpaired == 0 {
        e_multibranch = a + d;
    }

    let e = e_multibranch + e_sum;

    if helix {
        branches.pop();
    }

    let ij: Ij = branches.iter().map(|&(a, b)| (a as i32, b as i32)).collect();
    Struct::tagged(
        e,
        Desc::Bifurcation {
            unpaired: unpaired as i32,
            count: branches_count as u32,
        },
        ij,
    )
}

/// Traceback through the V and W caches to build the structure list.
fn traceback(s: &[u8], mut i: usize, mut j: usize, v_cache: &Cache, w_cache: &Cache) -> Vec<Struct> {
    let n = s.len();
    let s_w = w_cache[i][j].clone();
    if s_w.tag != Desc::Hairpin {
        while w_cache[i + 1][j] == s_w {
            i += 1;
        }
        while w_cache[i][j - 1] == s_w {
            j -= 1;
        }
    }

    let mut structs: Vec<Struct> = Vec::new();
    loop {
        let mut s_v = v_cache[i][j].clone();

        // multibranch structures are only in the w_cache
        if s_w.ij.len() > 1 {
            s_v = s_w.clone();
        }

        // render the label now, while both (i, j) and the inner pair are known
        let desc = render_desc(&s_v.tag, s, i, j, s_v.ij.first().copied(), n);
        structs.push(Struct::new(s_v.e, &desc, smallvec![(i as i32, j as i32)]));

        // it's a multibranch
        if s_v.ij.len() > 1 {
            let mut e_sum = 0.0;
            structs = trackback_energy(&structs);
            let mut branches: Vec<Struct> = Vec::new();
            for &(i1, j1) in &s_v.ij {
                let tb = traceback(s, i1 as usize, j1 as usize, v_cache, w_cache);
                if !tb.is_empty() && !tb[0].ij.is_empty() {
                    let (i2, j2) = tb[0].ij[0];
                    e_sum += w_cache[i2 as usize][j2 as usize].e;
                    branches.extend(tb);
                }
            }

            let last = structs.last().unwrap().clone();
            let m = structs.len();
            structs[m - 1] = Struct::new(pyround(last.e - e_sum, 1), &last.desc, last.ij.clone());
            structs.extend(branches);
            return structs;
        }

        // it's a stack, bulge, etc -- another single structure beyond this
        if s_v.ij.len() == 1 {
            let (ni, nj) = s_v.ij[0];
            i = ni as usize;
            j = nj as usize;
            continue;
        }

        // it's a hairpin, end of structure
        return trackback_energy(&structs);
    }
}

/// Add energy to each structure based on how its W(i,j) differs from the next.
fn trackback_energy(structs: &[Struct]) -> Vec<Struct> {
    let len = structs.len();
    let mut out: Vec<Struct> = Vec::with_capacity(len);
    for (index, st) in structs.iter().enumerate() {
        let e_next = if index == len - 1 {
            0.0
        } else {
            structs[index + 1].e
        };
        let e_corrected = pyround(st.e - e_next, 1);
        out.push(Struct::new(e_corrected, &st.desc, st.ij.clone()));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::energies;

    #[test]
    fn test_pair() {
        let seq = b"ATGGAATAGTG";
        assert_eq!(pair_str(seq, 0, 1, 9, 10), "AT/TG");
    }

    #[test]
    fn test_stack() {
        let seq = b"GCUCAGCUGGGAGAGC";
        let temp = 310.15;
        let e = stack_e(seq, 1, 2, 14, 13, temp, energies::rna());
        assert!((e - -2.1).abs() <= 0.1, "got {}", e);
    }

    #[test]
    fn test_bulge() {
        let seq = b"ACCCCCATCCTTCCTTGAGTCAAGGGGCTCAA";
        let e = bulge(seq, 5, 7, 18, 17, 310.15, energies::dna());
        assert!((3.22 - e).abs() <= 0.4, "got {}", e);
    }

    #[test]
    fn test_hairpin() {
        let temp = 310.15;
        let e = hairpin(b"ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA", 11, 16, temp, energies::dna());
        assert!((e - 4.3).abs() <= 1.0, "got {}", e);

        let e = hairpin(
            b"ACCCGCAAGCCCTCCTTCCTTGGATCAAGGGGCTCAA",
            3,
            8,
            temp,
            energies::dna(),
        );
        assert!((0.67 - e).abs() <= 0.1, "got {}", e);

        let e = hairpin(b"CUUUGCACG", 0, 8, temp, energies::rna());
        assert!((4.5 - e).abs() <= 0.2, "got {}", e);
    }

    #[test]
    fn test_internal_loop() {
        let seq = b"ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA";
        let (i, j) = (6i64, 21i64);
        let e = internal_loop(seq, i, i + 4, j, j - 4, 310.15, energies::dna());
        assert!((e - 3.5).abs() <= 0.1, "got {}", e);
    }

    #[test]
    fn test_w() {
        let cases: &[(&[u8], usize, &Energies)] = &[
            (b"GCUCAGCUGGGAGAGC", 15, energies::rna()),
            (b"CCUGCUUUGCACGCAGG", 16, energies::rna()),
            (b"GCGGUUCGAUCCCGC", 14, energies::rna()),
        ];
        let expected = [-3.8, -6.4, -4.2];
        for (idx, (seq, j, emap)) in cases.iter().enumerate() {
            let (_v, wc) = fill(seq, 310.15, emap);
            let e = wc[0][*j].e;
            assert!(
                (e - expected[idx]).abs() <= 0.2,
                "case {}: got {}",
                idx,
                e
            );
        }
    }
}
