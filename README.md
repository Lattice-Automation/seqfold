# seqfold

[![DOI](https://zenodo.org/badge/224018980.svg)](https://zenodo.org/badge/latestdoi/224018980)

Predict the minimum free energy structure of nucleic acids.

`seqfold` is an implementation of the `Zuker, 1981` dynamic programming algorithm, the basis for [UNAFold](http://unafold.rna.albany.edu/?q=DINAMelt/software)/[mfold](https://www.ibridgenetwork.org/#!/profiles/1045554571442/innovations/1/), with energy functions from `SantaLucia, 2004` (DNA) and `Turner, 2009` (RNA).

## Installation

### pypy3 (strongly recommended)

```bash
pypy3 -m ensurepip
pypy3 -m pip install seqfold
```

For a 200bp sequence (on my laptop), [pypy3](https://doc.pypy.org/en/latest/index.html) takes 2.5 seconds versus 15 seconds for CPython.

### Default pip

```bash
pip install seqfold
```

## Usage

### Python

```python
from seqfold import fold, dg, dg_cache, dot_bracket

# just returns minimum free energy
dg("GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC", temp = 37.0)  # -13.4

# `fold` returns a list of `seqfold.Struct` from the minimum free energy structure
structs = fold("GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC")
print(sum(s.e for s in structs))  # -13.4, same as dg()
for struct in structs:
    print(struct)  # prints the i, j, ddg, and description of each structure

# `dg_cache` returns a 2D array where each (i,j) combination returns the MFE from i to j inclusive
cache = dg_cache("GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC")

# `dot_bracket` returns a dot_bracket representation of the folding
print(dot_bracket(structs))  # ((((((((.((((......))))..((((.......)))).))))))))
```

### CLI

```txt
usage: seqfold [-h] [-t FLOAT] [-d] [-r] [--version] SEQ

Predict the minimum free energy (kcal/mol) of a nucleic acid sequence

positional arguments:
  SEQ                   nucleic acid sequence to fold

optional arguments:
  -h, --help            show this help message and exit
  -t FLOAT, --celcius FLOAT
                        temperature in Celsius
  -d, --dot-bracket     write a dot-bracket of the MFE folding to stdout
  -r, --sub-structures  write each substructure of the MFE folding to stdout
  --version             show program's version number and exit
```

#### Examples

```bash
$ seqfold GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC --celcius 32
-15.3
```

```bash
$ seqfold GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC --celcius 32 --dot-bracket --sub-structures
GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC
((((((((.((((......))))..((((.......)))).))))))))
   i    j    ddg  description
   0   48   -1.9  STACK:GG/CC
   1   47   -1.9  STACK:GG/CC
   2   46   -1.4  STACK:GA/CT
   3   45   -1.4  STACK:AG/TC
   4   44   -1.9  STACK:GG/CC
   5   43   -1.6  STACK:GT/CA
   6   42   -1.4  STACK:TC/AG
   7   41   -0.5  BIFURCATION:4n/3h
   9   22   -1.1  STACK:TT/AA
  10   21   -0.7  STACK:TA/AT
  11   20   -1.6  STACK:AC/TG
  12   19    3.0  HAIRPIN:CA/GG
  25   39   -1.9  STACK:CC/GG
  26   38   -2.3  STACK:CG/GC
  27   37   -1.9  STACK:GG/CC
  28   36    3.2  HAIRPIN:GT/CT
-15.3
```

### Notes

- The type of nucleic acid, DNA or RNA, is inferred from the input sequence.
- `seqfold` is case-insensitive with the input sequence.
- The default temperature is 37 degrees Celsius for both the Python and CLI interface.

## Motivation

Secondary structure prediction is used for making [PCR primers](https://academic.oup.com/nar/article/40/15/e115/1223759), designing [oligos for MAGE](https://pubs.acs.org/doi/abs/10.1021/acssynbio.5b00219), and tuning [RBS expression rates](https://www.sciencedirect.com/science/article/pii/B9780123851208000024).

While [UNAFold](http://unafold.rna.albany.edu/?q=DINAMelt/software) and [mfold](https://www.ibridgenetwork.org/#!/profiles/1045554571442/innovations/1/) are the most widely used applications for nucleic acid secondary structure prediction, their format and license are restrictive. `seqfold` is meant to be an open-source, minimalist alternative for predicting minimum free energy secondary structure.

|              | seqfold               | mfold                                                                                  | UNAFold                                                                                          |
| ------------ | --------------------- | -------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------ |
| License      | MIT                   | [Academic Non-commercial](http://unafold.rna.albany.edu/download/Academic_License.txt) | [\$200-36,000](https://www.ibridgenetwork.org/#!/profiles/1045554571442/innovations/1/products/) |
| OS           | Linux, MacOS, Windows | Linux, MacOS                                                                           | Linux, MacOS, Windows                                                                            |
| Format       | python, CLI python    | CLI binary                                                                             | CLI binary                                                                                       |
| Dependencies | none                  | (mfold_util)                                                                           | Perl, (gnuplot, glut/OpenGL)                                                                     |
| Graphical    | no                    | yes (output)                                                                           | yes (output)                                                                                     |
| Heterodimers | no                    | yes                                                                                    | yes                                                                                              |
| Constraints  | no                    | yes                                                                                    | yes                                                                                              |

## Citations

Papers, and how they helped in developing `seqfold`, are listed below.

### Nussinov, 1980

> Nussinov, Ruth, and Ann B. Jacobson. "Fast algorithm for predicting the secondary structure of single-stranded RNA." Proceedings of the National Academy of Sciences 77.11 (1980): 6309-6313.

Framework for the dynamic programming approach. It has a conceptually helpful "Maximal Matching" example that demonstrates the approach on a simple sequence with only matched or unmatched bp.

### Zuker, 1981

> Zuker, Michael, and Patrick Stiegler. "Optimal computer folding of large RNA sequences using thermodynamics and auxiliary information." Nucleic acids research 9.1 (1981): 133-148.

The most cited paper in this space. Extends further than `Nussinov, 1980` with a nearest neighbor approach to energies and a consideration of each of stack, bulge, internal loop, and hairpin. Their data structure and traceback method are both more intuitive than `Nussinov, 1980`.

### Jaeger, 1989

> Jaeger, John A., Douglas H. Turner, and Michael Zuker. "Improved predictions of secondary structures for RNA." Proceedings of the National Academy of Sciences 86.20 (1989): 7706-7710.

Zuker and colleagues expand on the 1981 paper to incorporate penalties for multibranched loops and dangling ends.

### SantaLucia, 2004

> SantaLucia Jr, John, and Donald Hicks. "The thermodynamics of DNA structural motifs." Annu. Rev. Biophys. Biomol. Struct. 33 (2004): 415-440.

The paper from which almost every DNA energy function in `seqfold` comes from (with the exception of multibranch loops). Provides neighbor entropies and enthalpies for stacks, mismatching stacks, terminal stacks, and dangling stacks. Ditto for bulges, internal loops, and hairpins.

### Turner, 2009

> Turner, Douglas H., and David H. Mathews. "NNDB: the nearest neighbor parameter database for predicting stability of nucleic acid secondary structure." Nucleic acids research 38.suppl_1 (2009): D280-D282.

Source of RNA nearest neighbor change in entropy and enthalpy parameter data. In `/data`.

### Ward, 2017

> Ward, M., Datta, A., Wise, M., & Mathews, D. H. (2017). Advanced multi-loop algorithms for RNA secondary structure prediction reveal that the simplest model is best. Nucleic acids research, 45(14), 8541-8550.

An investigation of energy functions for multibranch loops that validates the simple linear approach employed by `Jaeger, 1989` that keeps runtime within `O(n³)`.
