"""Test DNA/RNA folding."""

from os import path
from time import time
import unittest

from seqfold import dg, dg_cache, fold, Cache, Struct
from seqfold.dna import DNA_ENERGIES
from seqfold.rna import RNA_ENERGIES
from seqfold.fold import (
    STRUCT_DEFAULT,
    _traceback,
    _bulge,
    _stack,
    _hairpin,
    _internal_loop,
    _pair,
    _w,
    Structs,
)


class TestFold(unittest.TestCase):
    """Test folding functions"""

    def test_fold(self):
        """Fold function."""

        # it should throw if a nonsense sequence is provided
        with self.assertRaises(RuntimeError):
            dg("EASFEASFAST", 37.0)

        # Both U and T, mix of RNA and DNA
        with self.assertRaises(RuntimeError):
            dg("ATGCATGACGATUU", 37.0)

        # should not throw
        dg("ATGGATTTAGATAGAT")

    def test_fold_cache(self):
        """Gather a cache of the folded structure."""

        seq = "ATGGATTTAGATAGAT"
        cache = dg_cache(seq)
        seq_dg = dg(seq)

        self.assertAlmostEqual(seq_dg, cache[0][len(seq) - 1], delta=1)

    def test_fold_dna(self):
        """DNA folding to find min energy secondary structure."""

        # unafold's estimates for free energy estimates of DNA oligos
        unafold_dgs = {
            "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC": -10.94,  # three branched structure
            "GGGAGGTCGCTCCAGCTGGGAGGAGCGTTGGGGGTATATACCCCCAACACCGGTACTGATCCGGTGACCTCCC": -23.4,  # four branched structure
            "CGCAGGGAUACCCGCG": -3.8,
            "TAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGT": -6.85,
            "GGGGGCATAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGTCTGCGGTTCGATCCCGCGCGCTCCCACCA": -15.50,
            "TGAGACGGAAGGGGATGATTGTCCCCTTCCGTCTCA": -18.10,
            "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA": -3.65,
        }

        for seq, ufold in unafold_dgs.items():
            d = dg(seq, temp=37.0)

            # accepting a 50% difference
            delta = abs(0.5 * min(d, ufold))
            self.assertAlmostEqual(d, ufold, delta=delta)

    def test_fold_rna(self):
        """RNA folding to find min energy secondary structure."""

        # unafold's estimates for free energy estimates of RNA oligos
        # most tests available at https://github.com/jaswindersingh2/SPOT-RNA/blob/master/sample_inputs/batch_seq.fasta
        unafold_dgs = {
            "ACCCCCUCCUUCCUUGGAUCAAGGGGCUCAA": -9.5,
            "AAGGGGUUGGUCGCCUCGACUAAGCGGCUUGGAAUUCC": -10.1,
            "UUGGAGUACACAACCUGUACACUCUUUC": -4.3,
            "AGGGAAAAUCCC": -3.3,
            "GCUUACGAGCAAGUUAAGCAAC": -4.6,
            "UGGGAGGUCGUCUAACGGUAGGACGGCGGACUCUGGAUCCGCUGGUGGAGGUUCGAGUCCUCCCCUCCCAGCCA": -32.8,
            "GGGCGAUGAGGCCCGCCCAAACUGCCCUGAAAAGGGCUGAUGGCCUCUACUG": -20.7,
            "GGGGGCAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCGCGCUCCCACCA": -31.4,
        }

        for seq, ufold in unafold_dgs.items():
            d = dg(seq, temp=37.0)

            # accepting a 5% difference
            delta = abs(0.5 * min(d, ufold))
            self.assertAlmostEqual(d, ufold, delta=delta)

    def test_multibranch(self):
        """Fold a multibranch structure."""

        seq = "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC"  # three branch

        structs = fold(seq)
        self.assertTrue(
            any("BIFURCATION" in s.desc and (7, 41) in s.ij for s in structs)
        )

    def test_pair(self):
        """Create a pair for stack checking."""

        seq = "ATGGAATAGTG"
        self.assertEqual(_pair(seq, 0, 1, 9, 10), "AT/TG")

    def test_stack(self):
        """Calc delta G of a stack."""

        seq = "GCUCAGCUGGGAGAGC"
        temp = 310.15

        self.assertAlmostEqual(
            _stack(seq, 1, 2, 14, 13, temp, RNA_ENERGIES), -2.1, delta=0.1
        )

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

        seq = "CUUUGCACG"
        i = 0
        j = 8
        hairpin_dg = _hairpin(seq, i, j, temp, RNA_ENERGIES)
        self.assertAlmostEqual(4.5, hairpin_dg, delta=0.2)

    def test_internal_loop(self):
        """Calc dg of an internal loop."""

        seq = "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA"
        i = 6
        j = 21
        temp = 310.15
        dg = _internal_loop(seq, i, i + 4, j, j - 4, temp, DNA_ENERGIES)
        self.assertAlmostEqual(dg, 3.5, delta=0.1)

    def test_w(self):
        """Calculate _w over some range."""

        seq = "GCUCAGCUGGGAGAGC"
        i = 0
        j = 15
        temp = 310.15
        v_cache = []
        w_cache = []
        for _ in range(len(seq)):
            v_cache.append([STRUCT_DEFAULT] * len(seq))
            w_cache.append([STRUCT_DEFAULT] * len(seq))
        struct = _w(seq, i, j, temp, v_cache, w_cache, RNA_ENERGIES)
        self.assertAlmostEqual(struct.e, -3.8, delta=0.2)

        seq = "CCUGCUUUGCACGCAGG"
        i = 0
        j = 16
        temp = 310.15
        v_cache = []
        w_cache = []
        for _ in range(len(seq)):
            v_cache.append([STRUCT_DEFAULT] * len(seq))
            w_cache.append([STRUCT_DEFAULT] * len(seq))
        struct = _w(seq, i, j, temp, v_cache, w_cache, RNA_ENERGIES)
        self.assertAlmostEqual(struct.e, -6.4, delta=0.2)

        seq = "GCGGUUCGAUCCCGC"
        i = 0
        j = 14
        v_cache = []
        w_cache = []
        for _ in range(len(seq)):
            v_cache.append([STRUCT_DEFAULT] * len(seq))
            w_cache.append([STRUCT_DEFAULT] * len(seq))
        struct = _w(seq, i, j, temp, v_cache, w_cache, RNA_ENERGIES)
        self.assertAlmostEqual(struct.e, -4.2, delta=0.2)

    def _debug(self, cache: Structs):
        """Log the contents of a Cache."""

        rows = []
        for row in cache:
            rows.append(
                ",".join(str(s.ij).replace(",", "-") if s else "." for s in row)
            )
        rows.append("")
        for row in cache:
            rows.append(",".join(str(s.e) if s else "." for s in row))
        rows.append("")
        for row in cache:
            rows.append(",".join(str(s.desc) if s else "." for s in row))
        print("\n".join(rows))

