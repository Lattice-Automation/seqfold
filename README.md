# seqfold

Predict the minimum free energy structure of nucleic acids.

`seqfold` is an implementation of the `Zuker, 1981` dynamic programming algorithm, that is the basis for [UNAFold](http://unafold.rna.albany.edu/?q=DINAMelt/software)/[mfold](https://www.ibridgenetwork.org/#!/profiles/1045554571442/innovations/1/), plus energy functions from `SantaLucia, 2004`.

```python
from seqfold import fold

# a triple-branch DNA structure
fold("GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC")  # -10.9
```

## Motivation

Knowing about nucleic acid secondary structure is essential in synbio for selecting [primers for PCR](https://academic.oup.com/nar/article/40/15/e115/1223759), designing [oligos for MAGE](https://pubs.acs.org/doi/abs/10.1021/acssynbio.5b00219), and tuning [RBS expression rates](https://www.sciencedirect.com/science/article/pii/B9780123851208000024).

While [UNAFold](http://unafold.rna.albany.edu/?q=DINAMelt/software) and [mfold](https://www.ibridgenetwork.org/#!/profiles/1045554571442/innovations/1/) are the most widely used applications for nucleic acid secondary structure prediction, they have several limitations. Those and others are noted below.

|              | seqfold             | mfold                                                                                  | UNAFold                                                                                          |
| ------------ | ------------------- | -------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------ |
| License      | MIT                 | [Academic Non-commercial](http://unafold.rna.albany.edu/download/Academic_License.txt) | [\$200-36,000](https://www.ibridgenetwork.org/#!/profiles/1045554571442/innovations/1/products/) |
| OS           | Linux,MacOS,Windows | Linux,MacOS                                                                            | Linux,MacOS,Windows                                                                              |
| Format       | python,CLI python   | CLI binary                                                                             | CLI binary                                                                                       |
| Dependencies | none                | (mfold_util)                                                                           | Perl,(gnuplot,GD Library,glut/OpenGL)                                                            |
| Graphical    | no                  | [yes](http://unafold.rna.albany.edu/?q=mfold/documentation) (output)                   | yes (output)                                                                                     |
| Heterodimers | no                  | [yes](http://unafold.rna.albany.edu/?q=mfold/documentation)                            | yes                                                                                              |
| Constraints  | no                  | [yes](http://unafold.rna.albany.edu/?q=mfold/constraints)                              | yes                                                                                              |

## Citations

That papers and others that were used to develop this library are below. Each paper is listed along with how it relates to `seqfold`.

### Nussinov, 1980

Framework for the dynamic programming approach. It has a conceptually helpful "Maximal Matching" example that demonstrates the approach on a simple sequence with only matched or unmatched bp.

> Nussinov, Ruth, and Ann B. Jacobson. "Fast algorithm for predicting the secondary structure of single-stranded RNA." Proceedings of the National Academy of Sciences 77.11 (1980): 6309-6313.

### Zuker, 1981

The most cited paper in this space. Extends further than `Nussinov, 1980` with a nearest neighbor approach to energies and a consideration of each of stack, bulge, internal loop, and hairpin. Their data structure and traceback method are both more intuitive than `Nussinov, 1980`.

> Zuker, Michael, and Patrick Stiegler. "Optimal computer folding of large RNA sequences using thermodynamics and auxiliary information." Nucleic acids research 9.1 (1981): 133-148.

### Jaeger, 1989

Zuker and colleagues expand on the 1981 paper to incorporate penalties for multibranched loops and dangling ends.

> Jaeger, John A., Douglas H. Turner, and Michael Zuker. "Improved predictions of secondary structures for RNA." Proceedings of the National Academy of Sciences 86.20 (1989): 7706-7710.

### SantaLucia, 2004

The paper from which almost every DNA energy function in `seqfold` comes from (with the exception of multibranch loops). Provides neighbor entropies and enthalpies for stacks, mismatching stacks, terminal stacks, and dangling stacks. Ditto for bulges, internal loops, and hairpins.

> SantaLucia Jr, John, and Donald Hicks. "The thermodynamics of DNA structural motifs." Annu. Rev. Biophys. Biomol. Struct. 33 (2004): 415-440.

### Ward, 2017

An investigation of energy functions for multibranch loops that validates the simple linear approach employed by `Jaeger, 1989` that keeps runtime at `O(n³)`.

> Ward, M., Datta, A., Wise, M., & Mathews, D. H. (2017). Advanced multi-loop algorithms for RNA secondary structure prediction reveal that the simplest model is best. Nucleic acids research, 45(14), 8541-8550.
