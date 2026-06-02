#!/usr/bin/env python3
"""Smoke test for an installed seqfold (Rust-backed).

Verifies the public API and CLI against known-good values. Useful for checking
any install: after `maturin develop`, after `pip install seqfold`, or against a
wheel pulled from a CI run. Exits non-zero if anything is wrong.

    python scripts/smoke_test.py
"""
import shutil
import subprocess
import sys

import seqfold
from seqfold import fold, dg, dg_cache, dot_bracket, tm, tm_cache, gc_cache

fails = []


def check(name, got, want):
    ok = got == want
    print(f"  {'ok ' if ok else 'FAIL'}  {name}: {got!r}" + ("" if ok else f" != {want!r}"))
    if not ok:
        fails.append(name)


print(f"seqfold {seqfold.__version__}  ({seqfold._core.__file__.split('/')[-1]})")

# --- folding: free energy (dg) on known sequences ---
DNA = "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC"
check("dg DNA @37C", dg(DNA), -13.4)
check("dg DNA @32C", dg(DNA, temp=32.0), -15.3)
check(
    "dg E.coli 5S rRNA",
    dg(
        "UGCCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCC"
        "GUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGGCAU"
    ),
    -51.1,
)
check(
    "dg yeast tRNA-Phe",
    dg("GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA"),
    -21.9,
)

# --- structure: dot-bracket + Struct objects ---
structs = fold(DNA)
check(
    "dot_bracket",
    dot_bracket(DNA, structs),
    "((((((((.((((......))))..((((.......)))).))))))))",
)
check("sum(struct.e) == dg", round(sum(s.e for s in structs), 2), -13.4)
check("has a BIFURCATION", any("BIFURCATION" in s.desc for s in structs), True)

# --- dg_cache is an n x n matrix ---
cache = dg_cache(DNA)
check("dg_cache shape", (len(cache), len(cache[0])), (len(DNA), len(DNA)))

# --- Tm ---
check("tm", tm("AGTCTGGTCTGGATCTGAGAACTTCAGGCT"), 68.7)
check("gc_cache full-span", gc_cache("GGGGCCCC")[0][7], 1.0)
check("tm_cache invalid is inf", tm_cache("AGTCTGGTCTGG")[5][1], float("inf"))

# --- error handling: bad sequence must raise ---
try:
    dg("EASFEASF")
    check("bad seq raises", False, True)
except RuntimeError:
    check("bad seq raises RuntimeError", True, True)

# --- CLI end-to-end ---
cli = [shutil.which("seqfold")] if shutil.which("seqfold") else [sys.executable, "-m", "seqfold.main"]
out = subprocess.run(cli + [DNA, "--celcius", "32"], capture_output=True, text=True)
check("CLI dg @32C", out.stdout.strip(), "-15.3")
ver = subprocess.run(cli + ["--version"], capture_output=True, text=True)
check("CLI --version prefix", ver.stdout.strip().startswith("seqfold "), True)

print()
if fails:
    print(f"FAILED: {len(fails)} check(s): {', '.join(fails)}")
    sys.exit(1)
print("ALL CHECKS PASSED")
