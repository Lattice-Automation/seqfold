# seqfold

Find the minimum free energy structure of nucleic acids.

## Sources

This tool is an implementation of the `Zuker, 1981` dynamic programming algorithm that is the basis for Mfold/UNAFold — with energy functions from `SantaLucia, 2004`. That papers and others that were used to develop this library are below. Each paper is listed along with how it relates to `seqfold`.

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

> Advanced multi-loop algorithms for RNA secondary structure prediction reveal that the simplest model is best
