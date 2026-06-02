//! Quick wall-clock benchmark for the folding engine.
//! Run: `cargo run --release --example bench`

use std::time::Instant;

use seqfold::core::fold::dg;

const SEQS: &[&str] = &[
    // ~200bp DNA (the profiling sequence)
    "GAAATAGACGCCAAGTTCAATCCGTACTCCGACGTACGATGGAACAGTGTGGATGTGACGAGCTTCATTTATACCCTTCGCGCGCCGGACCGGGGTCCGCAAGGCGCGGCGGTGCACAAGCAATTGACAACTAACCACCGTGTATTCGTTATGGCACCAGGGAGTTTAAGCCGAGTCAATGGAGCTCGCAATACAGAGTT",
    // 118bp RNA
    "UGCCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGGCAU",
    // 98bp DNA
    "TGTCAGAAGTTTCCAAATGGCCAGCAATCAACCCATTCCATTGGGGATACAATGGTACAGTTTCGCATATTGTCGGTGAAAATGGTTCCATTAAACTCC",
    // 76bp DNA
    "GGGGGCATAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGTCTGCGGTTCGATCCCGCGCGCTCCCACCA",
];

fn main() {
    let iters = 5;
    let mut total = 0.0f64;
    // warm up energy tables
    let _ = dg(SEQS[3], 37.0).unwrap();

    for _ in 0..iters {
        let start = Instant::now();
        let mut checksum = 0.0;
        for s in SEQS {
            checksum += dg(s, 37.0).unwrap();
        }
        let el = start.elapsed().as_secs_f64();
        total += el;
        println!("iter: {:.4}s  checksum={}", el, checksum);
    }
    println!("avg over {} iters: {:.4}s", iters, total / iters as f64);
}
