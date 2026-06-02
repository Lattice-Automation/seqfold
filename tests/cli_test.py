"""Test the `seqfold` command-line interface end-to-end.

These run the actual console entry point (falling back to `python -m
seqfold.main`) as a subprocess, so they also exercise argument parsing and the
console-script wiring after a real install. Kept to a representative subset.
"""

import shutil
import subprocess
import sys
import unittest


def _cli():
    """The command prefix for invoking the CLI."""
    exe = shutil.which("seqfold")
    if exe:
        return [exe]
    return [sys.executable, "-m", "seqfold.main"]


def run_cli(*args):
    """Run the CLI and return (stdout, returncode)."""
    proc = subprocess.run(
        _cli() + list(args),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    return proc.stdout, proc.returncode


class TestCli(unittest.TestCase):
    """Test the command-line interface."""

    SEQ = "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC"

    def test_cli_dg(self):
        """Bare invocation prints the rounded free energy."""

        out, code = run_cli(self.SEQ)
        self.assertEqual(0, code)
        self.assertEqual("-13.4", out.strip())

    def test_cli_celcius(self):
        """The --celcius flag changes the temperature."""

        out, code = run_cli(self.SEQ, "--celcius", "32")
        self.assertEqual(0, code)
        self.assertEqual("-15.3", out.strip())

    def test_cli_dot_bracket(self):
        """The -d flag prints the sequence and a dot-bracket structure."""

        seq = "ACCCCCUCCUUCCUUGGAUCAAGGGGCUCAA"
        out, code = run_cli(seq, "-d")
        self.assertEqual(0, code)

        lines = out.splitlines()
        self.assertEqual(seq, lines[0])
        # second line is the dot-bracket, same length, only ()/. characters
        self.assertEqual(len(seq), len(lines[1]))
        self.assertTrue(set(lines[1]) <= set("()."))
        # last line is the free energy
        self.assertEqual("-9.5", lines[-1])

    def test_cli_sub_structures(self):
        """The -r flag prints a substructure table header."""

        out, code = run_cli(self.SEQ, "-r")
        self.assertEqual(0, code)
        self.assertIn("description", out.splitlines()[0])
        self.assertIn("BIFURCATION", out)
        self.assertEqual("-13.4", out.splitlines()[-1])

    def test_cli_version(self):
        """--version prints the program name and version."""

        out, code = run_cli("--version")
        self.assertEqual(0, code)
        self.assertTrue(out.strip().startswith("seqfold "))


if __name__ == "__main__":
    unittest.main()
