//! Melting-temperature (Tm) calculation. A port of `seqfold/tm.py`.

use super::energies::{self, Energies};
use super::pyfloat::pyround;

/// Error raised for invalid Tm inputs (maps to Python `ValueError`).
#[derive(Debug)]
pub struct TmError(pub String);

fn dna() -> &'static Energies {
    energies::dna()
}

/// Annealing temperature between `seq1` and `seq2` (its complement if empty).
pub fn tm(seq1: &str, seq2: &str, pcr: bool) -> Result<f64, TmError> {
    let (seq1, seq2) = parse_input(seq1, seq2)?;
    let s1 = seq1.as_bytes();
    let s2 = seq2.as_bytes();
    let emap = dna();

    // start with initiation enthalpy and entropy
    let (mut dh, mut ds) = emap.nn["init"];

    // add in initial A/T and initial G/Cs
    let init = [s1[0], s1[s1.len() - 1]];
    let init_at = init.iter().filter(|&&b| b == b'A' || b == b'T').count() as f64;
    let init_gc = init.iter().filter(|&&b| b == b'G' || b == b'C').count() as f64;
    let (init_at_h, init_at_s) = emap.nn["init_A/T"];
    let (init_gc_h, init_gc_s) = emap.nn["init_G/C"];
    dh += init_at * init_at_h + init_gc * init_gc_h;
    ds += init_at * init_at_s + init_gc * init_gc_s;

    // work through each nearest neighbor pair
    for i in 0..(s1.len() - 1) {
        let pair = format!(
            "{}{}/{}{}",
            s1[i] as char,
            s1[i + 1] as char,
            s2[i] as char,
            s2[i + 1] as char
        );

        let mut pair_dh = 0.0;
        let mut pair_ds = 0.0;
        if let Some(&(h, s)) = emap.nn.get(&pair) {
            pair_dh = h;
            pair_ds = s;
        } else if let Some(&(h, s)) = emap.internal_mm.get(&pair) {
            pair_dh = h;
            pair_ds = s;
        }

        // overwrite if it's a terminal pair
        if i == 0 || i == s1.len() - 2 {
            if let Some(&(h, s)) = emap.terminal_mm.get(&pair) {
                pair_dh = h;
                pair_ds = s;
            }
        }

        dh += pair_dh;
        ds += pair_ds;
    }

    let gc = gc_ratio(s1);
    Ok(calc_tm(dh, ds, pcr, gc, s1.len()))
}

/// A Tm matrix where (i, j) is the Tm of the subspan from i to j inclusive.
pub fn tm_cache(seq1: &str, seq2: &str, pcr: bool) -> Result<Vec<Vec<f64>>, TmError> {
    let (seq1, seq2) = parse_input(seq1, seq2)?;
    let s1 = seq1.as_bytes();
    let s2 = seq2.as_bytes();
    let emap = dna();
    let n = s1.len();

    let arr_gc = gc_cache_bytes(s1);
    let mut arr_dh = vec![vec![0.0f64; n]; n];
    let mut arr_ds = vec![vec![0.0f64; n]; n];
    let mut arr_tm = vec![vec![f64::INFINITY; n]; n];

    // fill in the diagonal
    for i in 0..n {
        if i == n - 1 {
            arr_dh[i][i] = arr_dh[i - 1][i - 1];
            arr_ds[i][i] = arr_ds[i - 1][i - 1];
            continue;
        }

        let pair = format!(
            "{}{}/{}{}",
            s1[i] as char,
            s1[i + 1] as char,
            s2[i] as char,
            s2[i + 1] as char
        );
        let (dh, ds) = emap
            .nn
            .get(&pair)
            .copied()
            .unwrap_or_else(|| emap.internal_mm[&pair]);

        arr_dh[i][i] = dh;
        arr_ds[i][i] = ds;
    }

    // fill in the tm array
    for i in 0..n {
        for j in (i + 1)..n {
            arr_dh[i][j] = arr_dh[i][j - 1] + arr_dh[j][j];
            arr_ds[i][j] = arr_ds[i][j - 1] + arr_ds[j][j];
            arr_tm[i][j] = calc_tm(arr_dh[i][j], arr_ds[i][j], pcr, arr_gc[i][j], j - i + 1);
        }
    }

    Ok(arr_tm)
}

/// GC ratio of each (i, j) range in the sequence.
pub fn gc_cache(seq: &str) -> Vec<Vec<f64>> {
    gc_cache_bytes(seq.as_bytes())
}

fn gc_cache_bytes(seq: &[u8]) -> Vec<Vec<f64>> {
    let n = seq.len();
    let mut arr_gc = vec![vec![f64::INFINITY; n]; n];

    // fill in the diagonal
    for i in 0..n {
        if i == n - 1 {
            arr_gc[i][i] = arr_gc[i - 1][i - 1];
            continue;
        }

        arr_gc[i][i] = if seq[i] == b'G' || seq[i] == b'C' {
            1.0
        } else {
            0.0
        };

        if i == n - 2 && arr_gc[i][i] == 0.0 {
            arr_gc[i][i] = if seq[i + 1] == b'G' || seq[i + 1] == b'C' {
                1.0
            } else {
                0.0
            };
        }
    }

    // fill in the upper right of the array
    for i in 0..n {
        for j in (i + 1)..n {
            arr_gc[i][j] = arr_gc[i][j - 1] + arr_gc[j][j];
        }
    }

    // convert to ratios
    for i in 0..n {
        for j in i..n {
            arr_gc[i][j] = pyround(arr_gc[i][j] / (j - i + 1) as f64, 1);
        }
    }

    arr_gc
}

fn parse_input(seq1: &str, seq2: &str) -> Result<(String, String), TmError> {
    let seq1 = seq1.to_uppercase();
    let seq2 = if seq2.is_empty() {
        let emap = dna();
        seq1.bytes()
            .map(|c| emap.complement[&c] as char)
            .collect::<String>()
    } else {
        seq2.to_string()
    };

    if seq1.len() != seq2.len() {
        return Err(TmError(format!(
            "Length mismatch between seq1 {} and seq2 {}",
            seq1.len(),
            seq2.len()
        )));
    }

    if seq1.len() < 2 {
        return Err(TmError(format!(
            "Sequence, {}bp, is too short for tm calculation",
            seq1.len()
        )));
    }

    Ok((seq1, seq2))
}

#[allow(non_snake_case)]
fn calc_tm(dh: f64, ds: f64, pcr: bool, gc: f64, seq_len: usize) -> f64 {
    // adjust salt based on mode
    let (seq1_conc, seq2_conc, Na, K, Tris, Mg, dNTPs): (
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
    ) = if pcr {
        (250.0, 0.0, 0.0, 50.0, 2.0, 1.5, 0.2)
    } else {
        (25.0, 25.0, 50.0, 0.0, 0.0, 0.0, 0.0)
    };

    // salt correction for deltaS (Owczarzy et al., 2008)
    let Mon = Na + K + Tris / 2.0;
    let mut mg = Mg * 1e-3;
    let mon = Mon * 1e-3;

    // coefficients to a multi-variate from the paper
    let (mut a, b, c, mut d, e, f, mut g) = (3.92, -0.911, 6.26, 1.42, -48.2, 52.5, 8.31);

    if dNTPs > 0.0 {
        let dntps = dNTPs * 1e-3;
        let ka = 3e4; // dissociation constant for Mg:dNTP
        mg = (-(ka * dntps - ka * mg + 1.0)
            + ((ka * dntps - ka * mg + 1.0).powi(2) + 4.0 * ka * mg).sqrt())
            / (2.0 * ka);
    }
    if Mon > 0.0 {
        let r = mg.sqrt() / mon;
        if r < 0.22 {
            return (4.29 * gc / 100.0 - 3.95) * 1e-5 * mon.ln() + 9.4e-6 * mon.ln().powi(2);
        } else if r < 6.0 {
            a = 3.92 * (0.843 - 0.352 * mon.sqrt() * mon.ln());
            d = 1.42 * (1.279 - 4.03e-3 * mon.ln() - 8.03e-3 * mon.ln().powi(2));
            g = 8.31 * (0.486 - 0.258 * mon.ln() + 5.25e-3 * mon.ln().powi(3));
        }
    }
    let corr = (a
        + b * mg.ln()
        + (gc / 100.0) * (c + d * mg.ln())
        + (1.0 / (2.0 * (seq_len as f64 - 1.0))) * (e + f * mg.ln() + g * mg.ln().powi(2)))
        * 1e-5;

    // tm with concentration consideration
    let k = (seq1_conc - (seq2_conc / 2.0)) * 1e-9;
    let r_const = 1.9872;
    let mut est = (dh * 1000.0) / (ds + r_const * k.ln()) - 273.15;

    // add in salt correction
    est = 1.0 / (1.0 / (est + 273.15) + corr) - 273.1;

    pyround(est, 1)
}

fn gc_ratio(seq: &[u8]) -> f64 {
    let count = seq.iter().filter(|&&b| b == b'G' || b == b'C').count();
    count as f64 / seq.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calc_tm() {
        let cases = [
            ("GGGACCGCCT", 51.9),
            ("CCATTGCTACC", 42.7),
            ("GCAGTGGATGTGAGA", 55.1),
            ("CTGGTCTGGATCTGAGAACTTCAGG", 67.7),
            ("CTTAAGATATGAGAACTTCAACTAATGTGT", 59.7),
            ("AGTCTGGTCTGGATCTGAGAACTTCAGGCT", 71.6),
        ];
        for (seq, actual) in cases {
            let calc = tm(seq, "", true).unwrap();
            assert!((calc - actual).abs() <= 7.0, "{}: got {}", seq, calc);
        }
    }
}
