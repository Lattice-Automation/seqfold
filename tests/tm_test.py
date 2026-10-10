"""Test Tm calculation"""

import unittest

from seqfold import tm, tm_cache, gc_cache


class TestTm(unittest.TestCase):
    """Test tm functions"""

    def test_calc_tm(self):
        """Oligo tm calculation."""

        # values are from Table 1 of IDT's:
        # Owczarzy et al. (2008), Biochemistry 4 7: 5336-5353
        # with a 1.5mM Mg concentration which looks typical according to NEB
        experimental_tms = {
            "GGGACCGCCT": 51.9,
            "CCATTGCTACC": 42.7,
            "GCAGTGGATGTGAGA": 55.1,
            "CTGGTCTGGATCTGAGAACTTCAGG": 67.7,
            "CTTAAGATATGAGAACTTCAACTAATGTGT": 59.7,
            "AGTCTGGTCTGGATCTGAGAACTTCAGGCT": 71.6,
        }

        for seq, actual in experimental_tms.items():
            calc = tm(seq)
            self.assertAlmostEqual(calc, actual, delta=7)  # within 7 deg tm difference

    def test_tm_cache(self):
        """Create a cache for tms over ranges in the sequence."""

        seq = "AGTCTGGTCTGGATCTGAGAACTTCAGGCT"
        n = len(seq)

        cache = tm_cache(seq)

        self.assertAlmostEqual(cache[0][n - 1], 71.6, delta=3)
        self.assertLess(cache[3][n - 3], cache[0][n - 1])
        self.assertEqual(float("inf"), cache[5][5], "inf for invalid subsequences")
        self.assertEqual(float("inf"), cache[5][1], "inf for invalid subsequences")

    def test_tm_non_pcr(self):
        """Non-PCR salt conditions return a Tm, not the salt correction.

        Regression test for https://github.com/Lattice-Automation/seqfold/issues/25
        where `tm(..., pcr=False)` returned the raw correction term (~0.0002).
        """

        seq = "ACGAGCGAGCGGAGCGGGCG"
        non_pcr = tm(seq, pcr=False)
        self.assertGreater(non_pcr, 0)
        # 250 nM oligo with Mg melts higher than 25 nM oligo in Na alone
        self.assertGreater(tm(seq, pcr=True), non_pcr)

    def test_gc_cache(self):
        """Create a cache of GC ratios from i to j."""

        seq = "GGATTACCCAGATAGATAGAT"
        ranges = [(0, len(seq) - 1), (5, 9), (3, 15)]

        cache = gc_cache(seq)

        for s, e in ranges:
            est = cache[s][e]
            ss = seq[s : e + 1]
            gc_count = ss.count("G") + ss.count("C")
            gc_actual = float(gc_count) / (e - s + 1)

            self.assertAlmostEqual(gc_actual, est, delta=0.02)

