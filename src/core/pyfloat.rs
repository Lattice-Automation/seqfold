//! Python-compatible float helpers.
//!
//! `seqfold`'s exact outputs depend on Python's `round()` (round-half-to-even
//! on the true binary value) and `str(float)` (shortest round-tripping repr,
//! always with a decimal point). Rust's `{:.*}` formatting is likewise
//! correctly rounded with ties-to-even, and `{}` is a shortest round-trip, so
//! we build both helpers on top of the standard formatter and only patch the
//! couple of spots where the textual conventions differ.

/// Equivalent of Python's `round(x, ndigits)` for the small `ndigits` values
/// (0..=2) that seqfold uses.
pub fn pyround(x: f64, ndigits: usize) -> f64 {
    if !x.is_finite() {
        return x;
    }
    // Rust's float formatting is correctly rounded, ties-to-even, matching the
    // dtoa-based rounding CPython performs.
    format!("{:.*}", ndigits, x)
        .parse::<f64>()
        .expect("formatted float must parse")
}

/// Equivalent of Python 3's `str(float)` / `repr(float)`.
pub fn pyfloat_str(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x < 0.0 { "-inf" } else { "inf" }.to_string();
    }

    // Rust's `{}` already yields the shortest string that round-trips, the same
    // policy as CPython's repr. The only divergences for the magnitudes seqfold
    // produces are: Rust prints integers without a trailing ".0", and uses a
    // different exponent format. seqfold values are small, so only the trailing
    // ".0" case can occur in practice; handle the exponent case too for safety.
    let s = format!("{}", x);
    if s.contains('.') || s.contains('e') || s.contains('E') {
        s
    } else {
        format!("{}.0", s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_matches_python() {
        assert_eq!(pyround(-13.4499999, 2), -13.45);
        assert_eq!(pyround(2.675, 2), 2.67); // binary repr < 2.675
        assert_eq!(pyround(0.5, 0), 0.0); // ties to even
        assert_eq!(pyround(1.5, 0), 2.0); // ties to even
        assert_eq!(pyround(2.5, 0), 2.0); // ties to even
        assert_eq!(pyround(-0.25, 1), -0.2); // ties to even
    }

    #[test]
    fn float_str_matches_python() {
        assert_eq!(pyfloat_str(-1.8), "-1.8");
        assert_eq!(pyfloat_str(3.1), "3.1");
        assert_eq!(pyfloat_str(-0.2), "-0.2");
        assert_eq!(pyfloat_str(45.0), "45.0");
        assert_eq!(pyfloat_str(1600.0), "1600.0");
        assert_eq!(pyfloat_str(-13.4), "-13.4");
    }
}
