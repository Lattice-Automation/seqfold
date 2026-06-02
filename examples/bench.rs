//! Wall-clock benchmark over real, structured nucleic-acid sequences
//! (see examples/known_structures.fasta). Random sequences barely fold, so they
//! under-represent the work; these have genuine hairpins, stems, and
//! multi-branch junctions.
//!
//! Run: `cargo run --release --example bench`
//! Single-thread:  `RAYON_NUM_THREADS=1 cargo run --release --example bench`

use std::time::Instant;

use seqfold::core::fold::dg;

/// (name, sequence) — a representative subset of the known-structure set.
const SEQS: &[(&str, &str)] = &[
    ("5S rRNA E. coli (120nt, 3-way junction)",
     "UGCCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGGCAU"),
    ("structured RNA (109nt)",
     "GUGAAAGUGUACCUAGGGUUCCAGCCUAUUUGUAGGUGUUCGGACCGAGCGGUACAGGUAUAUUAUAUACCACACCUUAGGGACAAAAGCCCGAGAGGAUGGUUUCACG"),
    ("structured RNA (103nt)",
     "GGAAGUGUACCUAGGGAUCCACCUCGAGAGAGGAAGGACCAAGCGGUACAGGCCUACUUCGGUAGGUUACACCGUGGGGAUAAAAGACCCGUGGCAAGUUUCG"),
    ("structured DNA oligo (99nt)",
     "TGTCAGAAGTTTCCAAATGGCCAGCAATCAACCCATTCCATTGGGGATACAATGGTACAGTTTCGCATATTGTCGGTGAAAATGGTTCCATTAAACTCC"),
    ("yeast tRNA-Phe (76nt, 4-way junction)",
     "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA"),
    ("tRNA (76nt)",
     "GGGGGCAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCGCGCUCCCACCA"),
];

fn main() {
    let iters = 5;
    let _ = dg(SEQS[0].1, 37.0).unwrap(); // warm up the energy tables

    println!("per-sequence (best of {}):", iters);
    let mut total = 0.0f64;
    for (name, seq) in SEQS {
        let mut best = f64::INFINITY;
        let mut d = 0.0;
        for _ in 0..iters {
            let start = Instant::now();
            d = dg(seq, 37.0).unwrap();
            best = best.min(start.elapsed().as_secs_f64());
        }
        total += best;
        println!("  {:>8.2} ms  dg={:>7}  {}", best * 1000.0, d, name);
    }
    println!("sum of best: {:.4}s", total);
}
