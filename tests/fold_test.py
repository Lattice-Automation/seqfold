"""Test DNA/RNA folding (public API).

The private-internal unit tests (``_pair``, ``_stack``, ``_bulge``,
``_hairpin``, ``_internal_loop``, ``_w``) were ported to Rust and run via
``cargo test``. This module covers the public Python surface.
"""

import unittest

from seqfold import dg, dg_cache, dot_bracket, fold


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

            # accepting a 60% difference
            delta = abs(0.6 * min(d, ufold))
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
            "CAGCGCGGCGGGCGGGAGUCCGGCGCGCCCUCCAUCCCCGGCGGCGUCGGCAAGGAGUAG": -18.26,
        }

        for seq, ufold in unafold_dgs.items():
            d = dg(seq, temp=37.0)

            # accepting a 30% difference
            delta = abs(0.3 * min(d, ufold))
            self.assertAlmostEqual(d, ufold, delta=delta)

    def test_dot_bracket(self):
        """Get the dot bracket notation for a folded structure."""

        seq = "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC"
        structs = fold(seq)

        self.assertEqual(
            "((((((((.((((......))))..((((.......)))).))))))))",
            dot_bracket(seq, structs),
        )

        seq = "ACGCTCACCGTGCCCAGTGAGCGA"
        structs = fold(seq)
        self.assertEqual(len(seq), len(dot_bracket(seq, structs)))

    def test_multibranch(self):
        """Fold a multibranch structure."""

        seq = "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC"  # three branch

        structs = fold(seq)
        self.assertTrue(
            any("BIFURCATION" in s.desc and (7, 41) in s.ij for s in structs)
        )

        seq = "CAGCGCGGCGGGCGGGAGUCCGGCGCGCCCUCCAUCCCCGGCGGCGUCGGCAAGGAGUAG"

        structs = fold(seq)
        self.assertTrue(
            any("BIFURCATION" in s.desc and (2, 56) in s.ij for s in structs)
        )
