from os import path

from seqfold.fold import calc_dg, fold

import seqfold

print(seqfold.__file__)

DIR = path.dirname(path.realpath(__file__))

ufold_dna = {
    "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC": -10.94,  # three branched structure
    "GGGAGGTCGCTCCAGCTGGGAGGAGCGTTGGGGGTATATACCCCCAACACCGGTACTGATCCGGTGACCTCCC": -23.4,  # four branched structure
    "CGCAGGGAUACCCGCG": -3.8,
    "TAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGT": -6.85,
    "GGGGGCATAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGTCTGCGGTTCGATCCCGCGCGCTCCCACCA": -15.50,
    "TGAGACGGAAGGGGATGATTGTCCCCTTCCGTCTCA": -18.10,
    "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA": -3.65,
    "TGTCAGAAGTTTCCAAATGGCCAGCAATCAACCCATTCCATTGGGGATACAATGGTACAGTTTCGCATATTGTCGGTGAAAATGGTTCCATTAAACTCC": -9.35,
}

# writing results to examples for comparison
DNA_RESULTS = {}

for seq, ufold in ufold_dna.items():
    ss = fold(seq, temp=37.0)
    dg = round(sum(s.e for s in ss), 2)
    DNA_RESULTS[seq] = (dg, ufold)

# save DNA_RESULTS to examples
with open(path.join(DIR, "..", "examples", "dna.csv"), "w") as ex:
    ex.write("seqfold,UNAFold,seq\n")

    for seq, (sf, uf) in DNA_RESULTS.items():
        ex.write(",".join([str(sf), str(uf), seq]) + "\n")


ufold_rna = {
    "ACCCCCUCCUUCCUUGGAUCAAGGGGCUCAA": -9.5,
    "AAGGGGUUGGUCGCCUCGACUAAGCGGCUUGGAAUUCC": -10.1,
    "UUGGAGUACACAACCUGUACACUCUUUC": -4.3,
    "UGCCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGGCAU": -54.9,
    "AGGGAAAAUCCC": -3.3,
    "GCUUACGAGCAAGUUAAGCAAC": -4.6,
    "UGGGAGGUCGUCUAACGGUAGGACGGCGGACUCUGGAUCCGCUGGUGGAGGUUCGAGUCCUCCCCUCCCAGCCA": -32.8,
    "GGGCGAUGAGGCCCGCCCAAACUGCCCUGAAAAGGGCUGAUGGCCUCUACUG": -20.7,
    "GGGGGCAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCGCGCUCCCACCA": -31.4,
}


# writing results to examples for comparison
RNA_RESULTS = {}

for seq, ufold in ufold_rna.items():
    dg = calc_dg(seq, temp=37.0)
    RNA_RESULTS[seq] = (dg, ufold)

# save RNA_RESULTS to examples
with open(path.join(DIR, "..", "examples", "rna.csv"), "w") as ex:
    ex.write("seqfold,UNAFold,seq\n")

    for seq, (sf, uf) in RNA_RESULTS.items():
        ex.write(",".join([str(round(sf, 2)), str(uf), seq]) + "\n")

