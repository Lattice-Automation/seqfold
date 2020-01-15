"""Test DNA/RNA folding."""

from os import path
from time import time
import unittest

from seqfold import calc_dg
from seqfold.dna import DNA_ENERGIES
from seqfold.fold import _bulge, _stack, _hairpin, _internal_loop, _pair


DIR = path.dirname(path.realpath(__file__))


class TestFold(unittest.TestCase):
    """Test folding functions"""

    def test_fold(self):
        """Test fold function."""

        # it should throw if a non-sense sequence is provided
        with self.assertRaises(RuntimeError):
            calc_dg("EASFEASFAST", 37.0)

        # Both U and T, mix of RNA and DNA
        with self.assertRaises(RuntimeError):
            calc_dg("ATGCATGACGATUU", 37.0)

        # should not throw
        calc_dg("ATGGATTTAGATAGAT")

    def test_fold_dna(self):
        """DNA folding to find min energy secondary structure."""

        # unafold's estimates for free energy estimates of DNA oligos
        unafold_dgs = {
            "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC": -10.94,  # three branched structure
            "CGCAGGGAUACCCGCG": -3.8,
            "TAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGT": -6.85,
            "GGGGGCATAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGTCTGCGGTTCGATCCCGCGCGCTCCCACCA": -15.50,
            "TGAGACGGAAGGGGATGATTGTCCCCTTCCGTCTCA": -18.10,
            "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA": -3.65,
            "TGTCAGAAGTTTCCAAATGGCCAGCAATCAACCCATTCCATTGGGGATACAATGGTACAGTTTCGCATATTGTCGGTGAAAATGGTTCCATTAAACTCC": -9.35,
        }

        # writing results to examples for comparison
        results = {}

        for seq, ufold in unafold_dgs.items():
            dg = calc_dg(seq, temp=37.0)

            # accepting a 50% difference
            delta = abs(0.5 * min(dg, ufold))
            self.assertAlmostEqual(dg, ufold, delta=delta)

            # save result
            results[seq] = (dg, ufold)

        # save results to examples
        with open(path.join(DIR, "..", "examples", "dna.csv"), "w") as ex:
            ex.write("seqfold,UNAFold,seq\n")

            for seq, (sf, uf) in results.items():
                ex.write(",".join([str(round(sf, 2)), str(uf), seq]) + "\n")

    def test_fold_rna(self):
        """RNA folding to find min energy secondary structure."""

        # unafold's estimates for free energy estimates of RNA oligos
        # most tests available at https://github.com/jaswindersingh2/SPOT-RNA/blob/master/sample_inputs/batch_seq.fasta
        unafold_dgs = {
            "ACCCCCUCCUUCCUUGGAUCAAGGGGCUCAA": -9.5,
            "AAGGGGUUGGUCGCCUCGACUAAGCGGCUUGGAAUUCC": -10.1,
            "UUGGAGUACACAACCUGUACACUCUUUC": -4.3,
            "UGCCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGGCAU": -54.9,
            "AGGGAAAAUCCC": -3.3,
            "GCUUACGAGCAAGUUAAGCAAC": -4.6,
            "UGGGAGGUCGUCUAACGGUAGGACGGCGGACUCUGGAUCCGCUGGUGGAGGUUCGAGUCCUCCCCUCCCAGCCA": -32.8,
            "GGGCGAUGAGGCCCGCCCAAACUGCCCUGAAAAGGGCUGAUGGCCUCUACUG": -20.7,
            "GUUCUUAUCAAGAGAAGCAGAGGGACUGGCCCGACGAAGCUUCAGCAACCGGUGUAAUGGCGAAAGCCAUGACCAAGGUGCUAAAUCCAGCAAGCUCGAACAGCUUGGAAGAUAAGAACA": -46.2,
        }

        # writing results to examples for comparison
        results = {}

        for seq, ufold in unafold_dgs.items():
            dg = calc_dg(seq, temp=37.0)

            # accepting a 5% difference
            delta = abs(0.5 * min(dg, ufold))
            self.assertAlmostEqual(dg, ufold, delta=delta)

            # save result
            results[seq] = (dg, ufold)

        # save results to examples
        with open(path.join(DIR, "..", "examples", "rna.csv"), "w") as ex:
            ex.write("seqfold,UNAFold,seq\n")

            for seq, (sf, uf) in results.items():
                ex.write(",".join([str(round(sf, 2)), str(uf), seq]) + "\n")

    def test_pair(self):
        """Create a pair for stack checking."""

        seq = "ATGGAATAGTG"

        self.assertEqual(_pair(seq, 0, 1, 9, 10), "AT/TG")

    def test_bulge(self):
        """Calc delta G calc of a bulge."""

        # mock bulge of CAT on one side and AG on other
        # from pg 429 of SantaLucia, 2004
        seq = "ACCCCCATCCTTCCTTGAGTCAAGGGGCTCAA"

        pair_dg = _bulge(seq, 5, 7, 18, 17, 310.15, DNA_ENERGIES)
        self.assertAlmostEqual(3.22, pair_dg, delta=0.4)

    def test_hairpin(self):
        """Calc delta G of a hairpin structure."""

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
        """Calc dg of an internal loop."""

        seq = "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA"
        i = 6
        j = 21
        temp = 310.15

        loop_temp = _internal_loop(seq, i, i + 4, j, j - 4, temp, DNA_ENERGIES)

        self.assertAlmostEqual(loop_temp, 3.5, delta=0.1)
