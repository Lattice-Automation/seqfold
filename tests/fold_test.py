"""Test oligo functions."""

from os import path
import unittest

from seqfold.dna import DNA_ENERGIES
from seqfold.fold import calc_tm, calc_dg, _bulge, _pair, _hairpin, _internal_loop


DIR = path.dirname(path.realpath(__file__))


class TestOligos(unittest.TestCase):
    """Test oligo functions"""

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

    def test_fold(self):
        """Test fold function."""

        # it should throw is a non-sense sequence is provided
        with self.assertRaises(RuntimeError):
            calc_dg("EASFEASFAST", 37.0)

    def test_fold_dna(self):
        """Test DNA folding to find min energy secondary structure."""

        # unafold's estimates for free energy estimates of DNA oligos
        unafold_dgs = {
            "TAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGT": -6.85,  # W(i,j) bifurcation
            "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC": -10.94,  # three branched structure
            "GGGGGCATAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGTCTGCGGTTCGATCCCGCGCGCTCCCACCA": -15.50,
            "TGAGACGGAAGGGGATGATTGTCCCCTTCCGTCTCA": -18.10,
            "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA": -3.65,  # unafold == -3.65
            "TGTCAGAAGTTTCCAAATGGCCAGCAATCAACCCATTCCATTGGGGATACAATGGTACAGTTTCGCATATTGTCGGTGAAAATGGTTCCATTAAACTCC": -9.35,
        }

        # writing results to examples for comparison
        results = {}

        for seq, unafold_est in unafold_dgs.items():
            dg = calc_dg(seq, temp=37.0)

            # accepting a 25% difference
            delta = abs(0.25 * unafold_est)
            self.assertAlmostEqual(dg, unafold_est, delta=delta)

            # save result
            results[seq] = (dg, unafold_est)

        # save results to examples
        with open(path.join(DIR, "..", "examples", "dna.csv"), "w") as ex:
            ex.write("seqfold,UNAFold,seq\n")

            for seq, (sf, uf) in results.items():
                ex.write(",".join([str(round(sf, 2)), str(uf), seq]) + "\n")

    def test_bulge(self):
        """Test delta G calc of a bulge."""

        # mock bulge of CAT on one side and AG on other
        # from pg 429 of SantaLucia, 2004
        pair = "CT/GA"
        seq = "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA"  # nonsense sequence

        pair_dg = _bulge(pair, seq, 5, 7, "CTC", 310.15, DNA_ENERGIES)
        self.assertAlmostEqual(3.22, pair_dg, delta=0.4)

        # self.assertEqual(_bulge("CT/GA", seq, 5, 22, 310.15), 0.0)

    def test_pair(self):
        """Test delta G of pairs with and without mismatches."""

        pairs = [
            ("CT/GA", -1.28),
            ("GG/CC", -2.1),
            ("TC/AG", -1.3),
            ("GT/CG", -0.59),
            ("TC/GG", 0.08),
        ]
        seq = "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA"

        for pair, dg_actual in pairs:
            dg_est = _pair(pair, seq, 5, 27, 310.15, DNA_ENERGIES)
            self.assertAlmostEqual(dg_est, dg_actual, delta=0.02)

    def test_hairpin(self):
        """Test delta G of a hairpin structure."""

        # hairpin = "CCTTGG"
        seq = "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA"
        i = 11
        j = 16
        temp = 310.15
        hairpin_dg = _hairpin(seq, i, j, temp, DNA_ENERGIES)
        # this differs from Unafold
        self.assertAlmostEqual(hairpin_dg, 4.3, delta=1.0)

        # from page 428 of SantaLucia, 2004
        # hairpin = "CGCAAG"
        seq = "ACCCGCAAGCCCTCCTTCCTTGGATCAAGGGGCTCAA"
        i = 3
        j = 8
        hairpin_dg = _hairpin(seq, i, j, temp, DNA_ENERGIES)
        self.assertAlmostEqual(0.67, hairpin_dg, delta=0.1)

    def test_internal_loop(self):
        """Test internal loop."""

        seq = "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA"
        i = 6
        j = 21
        left = "TCCTT"
        right = "ATCAA"
        temp = 310.15
        temp_est = 3.5

        loop_temp = _internal_loop(seq, i, j, left, right, temp, DNA_ENERGIES)

        self.assertAlmostEqual(loop_temp, temp_est, delta=0.1)
