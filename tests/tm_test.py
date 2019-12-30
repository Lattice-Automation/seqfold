"""Test Tm calculation"""

import unittest

from seqfold.tm import calc_tm

class TestTm(unittest.TestCase):
    """Test tm functions"""

    def test_calc_tm(self):
        """Test oligo tm calculation."""

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
            calc = calc_tm(seq)
            self.assertAlmostEqual(calc, actual, delta=3)  # within 3 deg tm difference