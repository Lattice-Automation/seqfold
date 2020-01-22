"""Test Tm calculation"""

import unittest

from seqfold import tm, tm_cache


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
            self.assertAlmostEqual(calc, actual, delta=3)  # within 3 deg tm difference

    def test_tm_cache(self):
        """Create a cache for tms over ranges in the sequence."""

        seq = "AGTCTGGTCTGGATCTGAGAACTTCAGGCT"
        n = len(seq)

        cache = tm_cache(seq)

        self.assertAlmostEqual(cache[0][n - 1], 71.6, delta=3)
        self.assertLess(cache[3][n - 3], 71.6)
        self.assertEqual(float("inf"), cache[5][5], "inf for invalid subsequences")
        self.assertEqual(float("inf"), cache[5][1], "inf for invalid subsequences")
