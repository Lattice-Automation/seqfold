"""Functions for oligos' Tm and delta G"""

import argparse
import math
import sys
from typing import Dict, List, Tuple, Any

from . import __version__


DNA_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

DNA_NN = {
    "init": (0.2, -5.7),
    "init_G/C": (0.0, 0.0),
    "init_A/T": (2.2, 6.9),
    "sym": (0, -1.4),
    "AA/TT": (-7.6, -21.3),
    "AT/TA": (-7.2, -20.4),
    "TA/AT": (-7.2, -20.4),
    "CA/GT": (-8.5, -22.7),
    "GT/CA": (-8.4, -22.4),
    "CT/GA": (-7.8, -21.0),
    "GA/CT": (-8.2, -22.2),
    "CG/GC": (-10.6, -27.2),
    "GC/CG": (-9.8, -24.4),
    "GG/CC": (-8.0, -19.0),
}
"""
SantaLucia (1998)
A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics
"""

DNA_INTERNAL_MM = {
    "AG/TT": (1.0, 0.9),
    "AT/TG": (-2.5, -8.3),
    "CG/GT": (-4.1, -11.7),
    "CT/GG": (-2.8, -8.0),
    "GG/CT": (3.3, 10.4),
    "GG/TT": (5.8, 16.3),
    "GT/CG": (-4.4, -12.3),
    "GT/TG": (4.1, 9.5),
    "TG/AT": (-0.1, -1.7),
    "TG/GT": (-1.4, -6.2),
    "TT/AG": (-1.3, -5.3),
    "AA/TG": (-0.6, -2.3),
    "AG/TA": (-0.7, -2.3),
    "CA/GG": (-0.7, -2.3),
    "CG/GA": (-4.0, -13.2),
    "GA/CG": (-0.6, -1.0),
    "GG/CA": (0.5, 3.2),
    "TA/AG": (0.7, 0.7),
    "TG/AA": (3.0, 7.4),
    "AC/TT": (0.7, 0.2),
    "AT/TC": (-1.2, -6.2),
    "CC/GT": (-0.8, -4.5),
    "CT/GC": (-1.5, -6.1),
    "GC/CT": (2.3, 5.4),
    "GT/CC": (5.2, 13.5),
    "TC/AT": (1.2, 0.7),
    "TT/AC": (1.0, 0.7),
    "AA/TC": (2.3, 4.6),
    "AC/TA": (5.3, 14.6),
    "CA/GC": (1.9, 3.7),
    "CC/GA": (0.6, -0.6),
    "GA/CC": (5.2, 14.2),
    "GC/CA": (-0.7, -3.8),
    "TA/AC": (3.4, 8.0),
    "TC/AA": (7.6, 20.2),
    "AA/TA": (1.2, 1.7),
    "CA/GA": (-0.9, -4.2),
    "GA/CA": (-2.9, -9.8),
    "TA/AA": (4.7, 12.9),
    "AC/TC": (0.0, -4.4),
    "CC/GC": (-1.5, -7.2),
    "GC/CC": (3.6, 8.9),
    "TC/AC": (6.1, 16.4),
    "AG/TG": (-3.1, -9.5),
    "CG/GG": (-4.9, -15.3),
    "GG/CG": (-6.0, -15.8),
    "TG/AG": (1.6, 3.6),
    "AT/TT": (-2.7, -10.8),
    "CT/GT": (-5.0, -15.8),
    "GT/CT": (-2.2, -8.4),
    "TT/AT": (0.2, -1.5),
}
"""
Internal mismatch table (DNA)
Allawi & SantaLucia (1997), Biochemistry 36: 10581-10594
Allawi & SantaLucia (1998), Biochemistry 37: 9435-9444
Allawi & SantaLucia (1998), Biochemistry 37: 2170-2179 * 
Allawi & SantaLucia (1998), Nucl Acids Res 26: 2694-2701 *
Peyret et al. (1999), Biochemistry 38: 3468-3477 *
"""

DNA_TERMINAL_MM = {
    "AA/TA": (-3.1, -7.8),
    "TA/AA": (-2.5, -6.3),
    "CA/GA": (-4.3, -10.7),
    "GA/CA": (-8.0, -22.5),
    "AC/TC": (-0.1, 0.5),
    "TC/AC": (-0.7, -1.3),
    "CC/GC": (-2.1, -5.1),
    "GC/CC": (-3.9, -10.6),
    "AG/TG": (-1.1, -2.1),
    "TG/AG": (-1.1, -2.7),
    "CG/GG": (-3.8, -9.5),
    "GG/CG": (-0.7, -19.2),
    "AT/TT": (-2.4, -6.5),
    "TT/AT": (-3.2, -8.9),
    "CT/GT": (-6.1, -16.9),
    "GT/CT": (-7.4, -21.2),
    "AA/TC": (-1.6, -4.0),
    "AC/TA": (-1.8, -3.8),
    "CA/GC": (-2.6, -5.9),
    "CC/GA": (-2.7, -6.0),
    "GA/CC": (-5.0, -13.8),
    "GC/CA": (-3.2, -7.1),
    "TA/AC": (-2.3, -5.9),
    "TC/AA": (-2.7, -7.0),
    "AC/TT": (-0.9, -1.7),
    "AT/TC": (-2.3, -6.3),
    "CC/GT": (-3.2, -8.0),
    "CT/GC": (-3.9, -10.6),
    "GC/CT": (-4.9, -13.5),
    "GT/CC": (-3.0, -7.8),
    "TC/AT": (-2.5, -6.3),
    "TT/AC": (-0.7, -1.2),
    "AA/TG": (-1.9, -4.4),
    "AG/TA": (-2.5, -5.9),
    "CA/GG": (-3.9, -9.6),
    "CG/GA": (-6.0, -15.5),
    "GA/CG": (-4.3, -11.1),
    "GG/CA": (-4.6, -11.4),
    "TA/AG": (-2.0, -4.7),
    "TG/AA": (-2.4, -5.8),
    "AG/TT": (-3.2, -8.7),
    "AT/TG": (-3.5, -9.4),
    "CG/GT": (-3.8, -9.0),
    "CT/GG": (-6.6, -18.7),
    "GG/CT": (-5.7, -15.9),
    "GT/CG": (-5.9, -16.1),
    "TG/AT": (-3.9, -10.5),
    "TT/AG": (-3.6, -9.8),
}
"""
Terminal mismatch table (DNA)
SantaLucia & Peyret (2001) Patent Application WO 01/94611
"""

DNA_DE = {
    "AA/.T": (0.2, 2.3),
    "AC/.G": (-6.3, -17.1),
    "AG/.C": (-3.7, -10.0),
    "AT/.A": (-2.9, -7.6),
    "CA/.T": (0.6, 3.3),
    "CC/.G": (-4.4, -12.6),
    "CG/.C": (-4.0, -11.9),
    "CT/.A": (-4.1, -13.0),
    "GA/.T": (-1.1, -1.6),
    "GC/.G": (-5.1, -14.0),
    "GG/.C": (-3.9, -10.9),
    "GT/.A": (-4.2, -15.0),
    "TA/.T": (-6.9, -20.0),
    "TC/.G": (-4.0, -10.9),
    "TG/.C": (-4.9, -13.8),
    "TT/.A": (-0.2, -0.5),
    ".A/AT": (-0.7, -0.8),
    ".C/AG": (-2.1, -3.9),
    ".G/AC": (-5.9, -16.5),
    ".T/AA": (-0.5, -1.1),
    ".A/CT": (4.4, 14.9),
    ".C/CG": (-0.2, -0.1),
    ".G/CC": (-2.6, -7.4),
    ".T/CA": (4.7, 14.2),
    ".A/GT": (-1.6, -3.6),
    ".C/GG": (-3.9, -11.2),
    ".G/GC": (-3.2, -10.4),
    ".T/GA": (-4.1, -13.1),
    ".A/TT": (2.9, 10.4),
    ".C/TG": (-4.4, -13.1),
    ".G/TC": (-5.2, -15.0),
    ".T/TA": (-3.8, -12.6),
}
"""DNA dangling ends

Bommarito et al. (2000), Nucl Acids Res 28: 1929-1934
"""

# the energies are the same for each loop stack in the
# reverse complementary direction
DNA_NN.update({k[::-1]: v for k, v in DNA_NN.items()})
DNA_INTERNAL_MM.update(
    {k[::-1]: v for k, v in DNA_INTERNAL_MM.items() if k[::-1] not in DNA_INTERNAL_MM}
)
DNA_TERMINAL_MM.update(
    {k[::-1]: v for k, v in DNA_TERMINAL_MM.items() if k[::-1] not in DNA_TERMINAL_MM}
)
DNA_DE.update({k[::-1]: v for k, v in DNA_DE.items() if k[::-1] not in DNA_DE})

DNA_TRI_TETRA_LOOPS = {
    "AGAAT": (-1.5, 0.0),
    "AGCAT": (-1.5, 0.0),
    "AGGAT": (-1.5, 0.0),
    "AGTAT": (-1.5, 0.0),
    "CGAAG": (-2.0, 0.0),
    "CGCAG": (-2.0, 0.0),
    "CGGAG": (-2.0, 0.0),
    "CGTAG": (-2.0, 0.0),
    "GGAAC": (-2.0, 0.0),
    "GGCAC": (-2.0, 0.0),
    "GGGAC": (-2.0, 0.0),
    "GGTAC": (-2.0, 0.0),
    "TGAAA": (-1.5, 0.0),
    "TGCAA": (-1.5, 0.0),
    "TGGAA": (-1.5, 0.0),
    "TGTAA": (-1.5, 0.0),
    "AAAAAT": (0.5, 0.6),
    "AAAACT": (0.7, -1.6),
    "AAACAT": (1.0, -1.6),
    "ACTTGT": (0.0, -4.2),
    "AGAAAT": (-1.1, -1.6),
    "AGAGAT": (-1.1, -1.6),
    "AGATAT": (-1.5, -1.6),
    "AGCAAT": (-1.6, -1.6),
    "AGCGAT": (-1.1, -1.6),
    "AGCTTT": (0.2, -1.6),
    "AGGAAT": (-1.1, -1.6),
    "AGGGAT": (-1.1, -1.6),
    "AGGGGT": (0.5, -0.6),
    "AGTAAT": (-1.6, -1.6),
    "AGTGAT": (-1.1, -1.6),
    "AGTTCT": (0.8, -1.6),
    "ATTCGT": (-0.2, -1.6),
    "ATTTGT": (0.0, -1.6),
    "ATTTTT": (-0.5, -1.6),
    "CAAAAG": (0.5, 1.3),
    "CAAACG": (0.7, 0.0),
    "CAACAG": (1.0, 0.0),
    "CAACCG": (0.0, 0.0),
    "CCTTGG": (0.0, -2.6),
    "CGAAAG": (-1.1, 0.0),
    "CGAGAG": (-1.1, 0.0),
    "CGATAG": (-1.5, 0.0),
    "CGCAAG": (-1.6, 0.0),
    "CGCGAG": (-1.1, 0.0),
    "CGCTTG": (0.2, 0.0),
    "CGGAAG": (-1.1, 0.0),
    "CGGGAG": (-1.0, 0.0),
    "CGGGGG": (0.5, 1.0),
    "CGTAAG": (-1.6, 0.0),
    "CGTGAG": (-1.1, 0.0),
    "CGTTCG": (0.8, 0.0),
    "CTTCGG": (-0.2, 0.0),
    "CTTTGG": (0.0, 0.0),
    "CTTTTG": (-0.5, 0.0),
    "GAAAAC": (0.5, 3.2),
    "GAAACC": (0.7, 0.0),
    "GAACAC": (1.0, 0.0),
    "GCTTGC": (0.0, -2.6),
    "GGAAAC": (-1.1, 0.0),
    "GGAGAC": (-1.1, 0.0),
    "GGATAC": (-1.6, 0.0),
    "GGCAAC": (-1.6, 0.0),
    "GGCGAC": (-1.1, 0.0),
    "GGCTTC": (0.2, 0.0),
    "GGGAAC": (-1.1, 0.0),
    "GGGGAC": (-1.1, 0.0),
    "GGGGGC": (0.5, 1.0),
    "GGTAAC": (-1.6, 0.0),
    "GGTGAC": (-1.1, 0.0),
    "GGTTCC": (0.8, 0.0),
    "GTTCGC": (-0.2, 0.0),
    "GTTTGC": (0.0, 0.0),
    "GTTTTC": (-0.5, 0.0),
    "GAAAAT": (0.5, 3.2),
    "GAAACT": (1.0, 0.0),
    "GAACAT": (1.0, 0.0),
    "GCTTGT": (0.0, -1.6),
    "GGAAAT": (-1.1, 0.0),
    "GGAGAT": (-1.1, 0.0),
    "GGATAT": (-1.6, 0.0),
    "GGCAAT": (-1.6, 0.0),
    "GGCGAT": (-1.1, 0.0),
    "GGCTTT": (-0.1, 0.0),
    "GGGAAT": (-1.1, 0.0),
    "GGGGAT": (-1.1, 0.0),
    "GGGGGT": (0.5, 1.0),
    "GGTAAT": (-1.6, 0.0),
    "GGTGAT": (-1.1, 0.0),
    "GTATAT": (-0.5, 0.0),
    "GTTCGT": (-0.4, 0.0),
    "GTTTGT": (-0.4, 0.0),
    "GTTTTT": (-0.5, 0.0),
    "TAAAAA": (0.5, -0.3),
    "TAAACA": (0.7, -1.6),
    "TAACAA": (1.0, -1.6),
    "TCTTGA": (0.0, -4.2),
    "TGAAAA": (-1.1, -1.6),
    "TGAGAA": (-1.1, -1.6),
    "TGATAA": (-1.6, -1.6),
    "TGCAAA": (-1.6, -1.6),
    "TGCGAA": (-1.1, -1.6),
    "TGCTTA": (0.2, -1.6),
    "TGGAAA": (-1.1, -1.6),
    "TGGGAA": (-1.1, -1.6),
    "TGGGGA": (0.5, -0.6),
    "TGTAAA": (-1.6, -1.6),
    "TGTGAA": (-1.1, -1.6),
    "TGTTCA": (0.8, -1.6),
    "TTTCGA": (-0.2, -1.6),
    "TTTTGA": (0.0, -1.6),
    "TTTTTA": (-0.5, -1.6),
    "TAAAAG": (0.5, 1.6),
    "TAAACG": (1.0, -1.6),
    "TAACAG": (1.0, -1.6),
    "TCTTGG": (0.0, -3.2),
    "TGAAAG": (-1.0, -1.6),
    "TGAGAG": (-1.0, -1.6),
    "TGATAG": (-1.5, -1.6),
    "TGCAAG": (-1.5, -1.6),
    "TGCGAG": (-1.0, -1.6),
    "TGCTTG": (-0.1, -1.6),
    "TGGAAG": (-1.0, -1.6),
    "TGGGAG": (-1.0, -1.6),
    "TGGGGG": (0.5, -0.6),
    "TGTAAG": (-1.5, -1.6),
    "TGTGAG": (-1.0, -1.6),
    "TTTCGG": (-0.4, -1.6),
    "TTTTAG": (-1.0, -1.6),
    "TTTTGG": (-0.4, -1.6),
    "TTTTTG": (-0.5, -1.6),
}
"""Experimental delta H and delta S for tri/tetra loops

Supplemental Material: Annu.Rev.Biophs.Biomol.Struct.33:415-40
doi: 10.1146/annurev.biophys.32.110601.141800
The Thermodynamics of DNA Structural Motifs
SantaLucia and Hicks, 2004

delta S was computed using delta G and delta H and is in cal / (K x mol)
(versus delta H in kcal / mol)
"""

DNA_MIN_HAIRPIN_LEN = 3
"""Cannot have extremely sharp angles in DNA.
This limit is from Nussinov, et al. (1980)
"""


DNA_INTERNAL_LOOPS: Dict[int, Tuple[float, float]] = {
    1: (0, 0),
    2: (0, 0),
    3: (0, -10.3),
    4: (0, -11.6),
    5: (0, -12.9),
    6: (0, -14.2),
    7: (0, -14.8),
    8: (0, -15.5),
    9: (0, -15.8),
    10: (0, -15.8),
    11: (0, -16.1),
    12: (0, -16.8),
    13: (0, -16.4),
    14: (0, -17.4),
    15: (0, -17.7),
    16: (0, -18.1),
    17: (0, -18.4),
    18: (0, -18.7),
    19: (0, -18.7),
    20: (0, -19.0),
    21: (0, -19.0),
    22: (0, -19.3),
    23: (0, -19.7),
    24: (0, -20.0),
    25: (0, -20.3),
    26: (0, -20.3),
    27: (0, -20.6),
    28: (0, -21.0),
    29: (0, -21.0),
    30: (0, -21.3),
}
"""Enthalpy and entropy increments for length dependence of internal loops

Were calculated from delta G Table 4 of SantaLucia, 2004:

Annu.Rev.Biophs.Biomol.Struct.33:415-40
doi: 10.1146/annurev.biophys.32.110601.141800
The Thermodynamics of DNA Structural Motifs
SantaLucia and Hicks, 2004

Additional loop sizes are accounted for with the Jacobson-Stockmayer
entry extrapolation formula in paper:
delta G (loop-n) = delta G (loop-x) + 2.44 x R x 310.15 x ln(n / x)

Additional correction is applied for asymmetric loops in paper:
delta G (asymmetry) = |length A - length B| x 0.3 (kcal / mol)
where A and B are lengths of both sides of loop
"""

DNA_BULGE_LOOPS: Dict[int, Tuple[float, float]] = {
    1: (0, -12.9),
    2: (0, -9.4),
    3: (0, -10.0),
    4: (0, -10.3),
    5: (0, -10.6),
    6: (0, -11.3),
    7: (0, -11.9),
    8: (0, -12.6),
    9: (0, -13.2),
    10: (0, -13.9),
    11: (0, -14.2),
    12: (0, -14.5),
    13: (0, -14.8),
    14: (0, -15.5),
    15: (0, -15.8),
    16: (0, -16.1),
    17: (0, -16.4),
    18: (0, -16.8),
    19: (0, -16.8),
    20: (0, -17.1),
    21: (0, -17.4),
    22: (0, -17.4),
    23: (0, -17.7),
    24: (0, -17.7),
    25: (0, -18.1),
    26: (0, -18.1),
    27: (0, -18.4),
    28: (0, -18.7),
    29: (0, -18.7),
    30: (0, -19.0),
}
"""Enthalpy and entropy increments for length depedence of bulge loops

Were calculated from delta G Table 4 of SantaLucia, 2004:

Annu.Rev.Biophs.Biomol.Struct.33:415-40
doi: 10.1146/annurev.biophys.32.110601.141800
The Thermodynamics of DNA Structural Motifs
SantaLucia and Hicks, 2004

For bulge loops of size 1, the intervening NN energy is used.
Closing AT penalty is applied on both sides
"""

DNA_HAIRPIN_LOOPS: Dict[int, Tuple[float, float]] = {
    1: (0, 0.0),
    2: (0, 0.0),
    3: (0, -11.3),
    4: (0, -11.3),
    5: (0, -10.6),
    6: (0, -12.9),
    7: (0, -13.5),
    8: (0, -13.9),
    9: (0, -14.5),
    10: (0, -14.8),
    11: (0, -15.5),
    12: (0, -16.1),
    13: (0, -16.1),
    14: (0, -16.4),
    15: (0, -16.8),
    16: (0, -17.1),
    17: (0, -17.4),
    18: (0, -17.7),
    19: (0, -18.1),
    20: (0, -18.4),
    21: (0, -18.7),
    22: (0, -18.7),
    23: (0, -19.0),
    24: (0, -19.3),
    25: (0, -19.7),
    26: (0, -19.7),
    27: (0, -19.7),
    28: (0, -20.0),
    29: (0, -20.0),
    30: (0, -20.3),
}
"""Enthalpy and entropy increments for length depedence of hairpin loops

Were calculated from delta G Table 4 of SantaLucia, 2004:

Annu.Rev.Biophs.Biomol.Struct.33:415-40
doi: 10.1146/annurev.biophys.32.110601.141800
The Thermodynamics of DNA Structural Motifs
SantaLucia and Hicks, 2004

For hairpins of length 3 and 4, the entropy values are looked up
in the DNA_TRI_TETRA_LOOPS Dict

From formula 8-9 of the paper:
An additional 1.6 delta entropy penalty if the hairpin is closed by AT
"""


def calc_tm(seq1: str, seq2: str = "", pcr: bool = True) -> float:
    """Calculate the annealing temperature between seq1 and seq2.

    If seq2 is not provided, it's exact complement is used.
    In otherwords, it's assumed to be an exact match. This tm
    calculate does not account for pseudoknots or anything other
    than an exact, unpadded alignment between seq1 and seq2.

    This is largley influenced by Bio.SeqUtils.MeltingTemp with
    some different defaults. Here, the reaction mixture is assumed to
    be PCR and concentrations for Mg, Tris, K, and dNTPs are included
    that match a typical PCR reaction according to Thermo and NEB. Additionally,
    the salt correction formula from IDT's Owczarzy et al. (2008) is used.

    NEB: https://www.neb.com/tools-and-resources/usage-guidelines/guidelines-for-pcr-optimization-with-taq-dna-polymerase
    ThermoFisher: https://www.thermofisher.com/order/catalog/product/18067017?SID=srch-srp-18067017

    NOTE: Sequences are assumed not to be symmetrical. Oligo not binding to self.

    Args:
        seq1: The seq whose tm is calculated

    Keyword Args:
        seq2: The seq that seq1 anneals to in 3' -> 5' direction
        pcr: Whether tm is being calculated for the oligo is in a
            PCR reaction mixture. If so, ion and Tris concentrations
            that match a typical NEB/ThermoFisher PCR mixture are used

    Returns:
        float: The estimated tm as a float
    """

    if not seq2:
        seq1 = seq1.upper()
        seq2 = "".join([DNA_COMPLEMENT[c] for c in seq1])

    if len(seq1) != len(seq2):
        raise ValueError(
            f"length mismatch between seq1 {len(seq1)} and seq2 {len(seq2)}"
        )

    # sum enthalpy and entropy. Enthaply is first value of each tuple and
    # entropy is the second value of each tuple in:
    # SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct 33: 415-440

    # start with initiation enthalpy and entropy
    delta_h, delta_s = DNA_NN["init"]

    # add in initial A/T and initial G/Cs
    init = seq1[0] + seq1[-1]
    init_at = init.count("A") + init.count("T")
    init_gc = init.count("G") + init.count("C")
    init_at_h, init_at_s = DNA_NN["init_A/T"]
    init_gc_h, init_gc_s = DNA_NN["init_G/C"]
    delta_h += init_at * init_at_h + init_gc * init_gc_h
    delta_s += init_at * init_at_s + init_gc * init_gc_s

    # work through each nearest neighbor pair
    for i in range(len(seq1) - 1):
        pair1 = seq1[i] + seq1[i + 1]
        pair2 = seq2[i] + seq2[i + 1]
        pair = pair1 + "/" + pair2

        # assuming internal neighbor pair
        pair_delta_h, pair_delta_s = 0.0, 0.0
        if pair in DNA_NN:
            pair_delta_h, pair_delta_s = DNA_NN[pair]
        elif pair in DNA_INTERNAL_MM:
            pair_delta_h, pair_delta_s = DNA_INTERNAL_MM[pair]

        # overwrite if it's a terminal pair
        if i in (0, len(seq1) - 2):
            if pair in DNA_TERMINAL_MM:
                pair_delta_h, pair_delta_s = DNA_TERMINAL_MM[pair]

        delta_h += pair_delta_h
        delta_s += pair_delta_s

    # adjust salt based on mode
    if pcr:
        seq1_conc = 250.0
        seq2_conc = 0.0
        Na = 0
        K = 50
        Tris = 2  # see Thermo
        Mg = 1.5  # see NEB
        dNTPs = 0.2  # see NEB
    else:
        seq1_conc = 25.0
        seq2_conc = 25.0
        Na = 50
        K = 0
        Tris = 0
        Mg = 0
        dNTPs = 0

    # salt correction for deltaS
    # copied-pasted from Bio.SeqUtils' use of a decision tree by:
    # Owczarzy et al. (2008), Biochemistry 4 7: 5336-5353
    Mon = Na + K + Tris / 2.0  # monovalent ions
    mg = Mg * 1e-3  # Lowercase ions (mg, mon, dntps) are molar
    mon = Mon * 1e-3

    # coefficients to a multi-variate from the paper
    a, b, c, d, e, f, g = 3.92, -0.911, 6.26, 1.42, -48.2, 52.5, 8.31

    if dNTPs > 0:
        dntps = dNTPs * 1e-3
        ka = 3e4  # Dissociation constant for Mg:dNTP
        # Free Mg2+ calculation:
        mg = (
            -(ka * dntps - ka * mg + 1.0)
            + math.sqrt((ka * dntps - ka * mg + 1.0) ** 2 + 4.0 * ka * mg)
        ) / (2.0 * ka)
    if Mon > 0:
        R = math.sqrt(mg) / mon
        if R < 0.22:
            corr = (4.29 * _gc(seq1) / 100 - 3.95) * 1e-5 * math.log(
                mon
            ) + 9.40e-6 * math.log(mon) ** 2
            return corr
        elif R < 6.0:
            a = 3.92 * (0.843 - 0.352 * math.sqrt(mon) * math.log(mon))
            d = 1.42 * (1.279 - 4.03e-3 * math.log(mon) - 8.03e-3 * math.log(mon) ** 2)
            g = 8.31 * (0.486 - 0.258 * math.log(mon) + 5.25e-3 * math.log(mon) ** 3)
    corr = (
        a
        + b * math.log(mg)
        + (_gc(seq1) / 100) * (c + d * math.log(mg))
        + (1 / (2.0 * (len(seq1) - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
    ) * 1e-5

    # tm with concentration consideration
    k = (seq1_conc - (seq2_conc / 2.0)) * 1e-9
    R = 1.9872
    tm = (delta_h * 1000.0) / (delta_s + R * math.log(k)) - 273.15

    # add in salt correction
    tm = 1 / (1 / (tm + 273.15) + corr) - 273.15

    return tm


Cache = List[List[Any]]
"""A map from i, j tuple to a value."""


def fold(seq: str, temp: float = 37.0) -> List[Tuple[int, int, str, float]]:
    """Fold the DNA sequence and return lowest free energy score.

    Based on the approach described in:
    Zuker and Stiegler, 1981
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC326673/pdf/nar00394-0137.pdf

    Args:
        seq: The sequence to fold

    Keyword Args:
        temp: The temperature the fold takes place in, in Celcius

    Returns:
        List[Tuple[int, int, str, float]]: A list of structures. Stacks, bulges, hairpins, etc.
            The first two ints are the index on the left and right of the structure (the
            closing indexes). Then a description of the structure followed by the delta delta G
            of the structure: the added free energy change to the overall secondary structure
    """

    seq = seq.upper()
    temp = temp + 273.15  # kelvin

    v_cache: Any = []
    w_cache: Any = []
    for _ in range(len(seq)):
        v_cache.append([(None, "", [])] * len(seq))
        w_cache.append([(None, "", [])] * len(seq))

    # for increasing fragment length; smaller fragments first
    for f_len in range(4, len(seq)):
        # for increasing start index
        for i in range(len(seq) - f_len):
            # fill
            _v(seq, i, i + f_len, temp, v_cache, w_cache)
            _w(seq, i, i + f_len, temp, v_cache, w_cache)

    # gather the min energy structure over the full sequence
    min_e, _, _ = _w(seq, 0, len(seq) - 1, temp, v_cache, w_cache)
    min_e = round(min_e, 2)

    # get the structure out of the cache
    # _debug(v_cache, w_cache)
    structs = _traceback(0, len(seq) - 1, v_cache, w_cache)
    total_e = sum(s[-1] for s in structs)
    assert abs(total_e - min_e) < 0.2, f"{total_e} != {min_e}"

    return structs


def calc_dg(seq: str, temp: float) -> float:
    """Fold the sequence and return just the delta G of the structure
    
    Args:
        seq: The sequence to fold
        temp: The temperature to fold at
    
    Returns:
        float: The minimum free energy of the folded sequence
    """

    structs = fold(seq, temp)

    return round(sum([s[-1] for s in structs]), 2)


def _v(
    seq: str, i: int, j: int, temp: float, v_cache: Cache, w_cache: Cache
) -> Tuple[float, str, List[Tuple[int, int]]]:
    """Find, store and return the minimum free energy of the structure between i and j

    If i and j don't bp, store and return INF.
    See: Figure 2B of Zuker, 1981

    Args:
        seq: The sequence being folded
        i: The start index
        j: The end index (inclusive)
        temp: The temperature in Kelvin
        v_cache: Free energy cache for if i and j bp. INF otherwise
        w_cache: Free energy cache for lowest energy structure from i to j. 0 otherwise

    Returns:
        float: The minimum energy folding structure possible between i and j on seq
    """

    if v_cache[i][j][0] is not None:
        return v_cache[i][j]

    # the ends must basepair for V(i,j)
    if DNA_COMPLEMENT[seq[i]] != seq[j]:
        v_cache[i][j] = (math.inf, "", [])
        return v_cache[i][j]

    # E1 = FH(i, j); hairpin
    pair = seq[i] + seq[i + 1] + "/" + seq[j] + seq[j - 1]
    e1, e1_type = _hairpin(seq, i, j, temp), "HAIRPIN:" + pair
    e1_ij: List[Tuple[int, int]] = []
    if j - i == 4:  # small hairpin; 4bp
        v_cache[i][j] = (e1, e1_type, [])
        w_cache[i][j] = (e1, e1_type, [])
        return v_cache[i][j]

    # E2 = min{FL(i, j, i', j') + V(i', j')}, i<i'<j'<j
    # stacking region or bulge or interior loop; Figure 2A(2)
    # j-i=d>4; various pairs i',j' for j'-i'<d
    e2 = math.inf
    e2_type = ""
    e2_ij: List[Tuple[int, int]] = []
    for i_1 in range(i + 1, j - 4):
        for j_1 in range(i_1 + 4, j):
            pair = seq[i] + seq[i_1] + "/" + seq[j] + seq[j_1]
            pair_left = seq[i] + seq[i + 1] + "/" + seq[j] + seq[j - 1]
            pair_right = seq[i_1 - 1] + seq[i_1] + "/" + seq[j_1 + 1] + seq[j_1]
            pair_outer = pair_left in DNA_NN or pair_right in DNA_NN

            stack = i_1 == i + 1 and j_1 == j - 1 and pair in DNA_NN
            bulge_left = i_1 > i + 1 and pair in DNA_NN
            bulge_right = j_1 < j - 1 and pair in DNA_NN

            loop_left = seq[i : i_1 + 1]
            loop_right = seq[j_1 : j + 1]

            e2_test, e2_test_type = math.inf, ""
            if stack:
                # it's a neighboring/stacking pair in a helix
                e2_test = _pair(pair, seq, i, j, temp)
                e2_test_type = "STACK:" + pair

                if i > 0 and j == len(seq) - 1 or i == 0 and j < len(seq) - 1:
                    # there's a dangling end
                    e2_test_type = "STACK_DE:" + pair
            elif bulge_left and bulge_right and not pair_outer:
                # it's an interior loop
                e2_test = _internal_loop(seq, i, j, loop_left, loop_right, temp)
                e2_test_type = "INTERIOR_LOOP"

                if len(loop_left) == 3 and len(loop_right) == 3:
                    # technically an interior loop of 1. really 1bp mismatch
                    e2_test_type = "STACK:" + loop_left + "/" + loop_right[::-1]
            elif bulge_left and not bulge_right:
                # it's a bulge on the left side
                e2_test = _bulge(pair, seq, i, j, loop_left, temp)
                e2_test_type = "BULGE:" + str(i) + "-" + str(i_1)
            elif bulge_right and not bulge_left:
                # it's a bulge on the right side
                e2_test = _bulge(pair, seq, i, j, loop_right, temp)
                e2_test_type = "BULGE:" + str(j_1) + "-" + str(j)
            else:
                # it's basically a hairpin, only outside bp match
                continue

            # add V(i', j')
            e2_test += _v(seq, i_1, j_1, temp, v_cache, w_cache)[0]
            if e2_test < e2:
                e2, e2_type, e2_ij = e2_test, e2_test_type, [(i_1, j_1)]

    # E3 = min{W(i+1,i') + W(i'+1,j-1)}, i+1<i'<j-2
    e3, e3_type = math.inf, "BIFURCATION"
    e3_ij: List[Tuple[int, int]] = []
    for i_1 in range(i + 2, j - 2):
        e3_test, unpaired, helixes = _multi_branch(
            seq, i + 1, i_1, j - 1, temp, v_cache, w_cache, True
        )

        if e3_test < e3:
            e3 = e3_test
            e3_ij = [(i + 1, i_1), (i_1 + 1, j - 1)]
            e3_type = "BIFURCATION:" + str(unpaired) + "n/" + str(helixes) + "h"

    e = min(
        [(e1, e1_type, e1_ij), (e2, e2_type, e2_ij), (e3, e3_type, e3_ij)],
        key=lambda x: x[0],
    )
    v_cache[i][j] = e
    return e


def _w(
    seq: str, i: int, j: int, temp: float, v_cache: Cache, w_cache: Cache
) -> Tuple[float, str, List[Tuple[int, int]]]:
    """Find and return the lowest free energy structure in Sij subsequence

    Figure 2B

    Args:
        seq: The sequence being folded
        i: The start index
        j: The end index (inclusive)
        temp: The temperature in Kelvin
        v_cache: Free energy cache for if i and j bp
        w_cache: Free energy cache for lowest energy structure from i to j. 0 otherwise

    Returns:
        float: The free energy for the subsequence from i to j
    """

    if w_cache[i][j][0] is not None:
        return w_cache[i][j]

    if j - i < 4:
        w_cache[i][j] = (math.inf, "", [])
        return w_cache[i][j]

    w1 = _w(seq, i + 1, j, temp, v_cache, w_cache)
    w2 = _w(seq, i, j - 1, temp, v_cache, w_cache)
    w3 = _v(seq, i, j, temp, v_cache, w_cache)

    w4, w4_type = math.inf, "BIFURCATION"
    w4_ij: List[Tuple[int, int]] = []
    for i_1 in range(i + 1, j - 1):
        w4_test, unpaired, helixes = _multi_branch(
            seq, i, i_1, j, temp, v_cache, w_cache, False
        )

        if w4_test < w4:
            w4 = w4_test
            w4_ij = [(i, i_1), (i_1 + 1, j)]
            w4_type = "BIFURCATION:" + str(unpaired) + "n/" + str(helixes) + "h"

    w = min([w1, w2, w3, (w4, w4_type, w4_ij)], key=lambda x: x[0])
    w_cache[i][j] = w
    return w


def _d_g(d_h: float, d_s: float, temp: float) -> float:
    """Find the free energy given delta h, s and temp

    Args:
        d_h: The enthalpy increment in kcal / mol
        d_s: The entropy increment in cal / mol
        temp: The temperature in Kelvin

    Returns:
        The free energy increment in kcal / (mol x K)
    """

    return d_h - temp * (d_s / 1000.0)


def _j_s(query_len: int, known_len: int, d_g_x: float, temp: float) -> float:
    """Estimate the free energy of length query_len based on one of length known_len.

    The Jacobson-Stockmayer entry extrapolation formula is used
    for bulges, hairpins, etc that fall outside the 30nt upper limit
    for pre-calculated free-energies. See SantaLucia and Hicks (2014).

    Args:
        query_len: Length of element without known free energy value
        known_len: Length of element with known free energy value (d_g_x)
        d_g_x: The free energy of the element known_len
        temp: Temperature in Kelvin

    Returns:
        float: The free energy for a structure of length query_len
    """

    gas_constant = 1.9872e-3
    return d_g_x + 2.44 * gas_constant * temp * math.log(query_len / float(known_len))


def _gc(seq: str) -> float:
    """Return the GC ratio of a sequence."""

    seq = seq.upper()
    return float(seq.count("G") + seq.count("C")) / float(len(seq))


def _pair(pair: str, seq: str, i: int, j: int, temp: float) -> float:
    """Get the free energy for a pair.

    Using the indexes i and j, check whether it's at the end of
    the sequence or internal. Then check whether it's a match
    or mismatch, and return.

    Two edge-cases are terminal mismatches and dangling ends.
    The energy of a dangling end is added to the energy of a pair
    where i XOR j is at the sequence's end.

    Args:
        pair: The pair sequence, ex: (AG/TC)
        seq: The full folding sequence
        i: The start index on left side of the pair/stack
        j: The end index on right side of the pair/stack
        temp: Temperature in Kelvin

    Returns:
        float: The free energy of the NN pairing
    """

    if i > 0 and j < len(seq) - 1:
        # it's internal
        d_h, d_s = DNA_NN[pair] if pair in DNA_NN else DNA_INTERNAL_MM[pair]
        return _d_g(d_h, d_s, temp)

    if i == 0 and j == len(seq) - 1:
        # it's terminal
        d_h, d_s = DNA_NN[pair] if pair in DNA_NN else DNA_TERMINAL_MM[pair]
        return _d_g(d_h, d_s, temp)

    if i > 0 and j == len(seq) - 1:
        # it's dangling on left
        d_h, d_s = DNA_NN[pair] if pair in DNA_NN else DNA_TERMINAL_MM[pair]
        d_g = _d_g(d_h, d_s, temp)

        pair_de = seq[i - 1] + seq[i] + "/." + seq[j]
        d_h, d_s = DNA_DE[pair_de]
        return d_g + _d_g(d_h, d_s, temp)

    if i == 0 and j < len(seq) - 1:
        # it's dangling on right
        d_h, d_s = DNA_NN[pair] if pair in DNA_NN else DNA_TERMINAL_MM[pair]
        d_g = _d_g(d_h, d_s, temp)

        pair_de = "." + seq[i] + "/" + seq[j + 1] + seq[j]
        d_h, d_s = DNA_DE[pair_de]
        return d_g + _d_g(d_h, d_s, temp)

    raise RuntimeError


def _hairpin(seq: str, i: int, j: int, temp: float) -> float:
    """Calculate the free energy of a hairpin.

    Args:
        seq: The sequence we're folding
        i: The index of start of hairpin
        j: The index of end of hairpin
        temp: Temperature in Kelvin

    Returns:
        float: The free energy increment from the hairpin structure
    """

    if j - i < 4:
        return math.inf

    hairpin = seq[i : j + 1]
    hairpin_len = len(hairpin) - 2
    pair = hairpin[0] + hairpin[1] + "/" + hairpin[-1] + hairpin[-2]

    if DNA_COMPLEMENT[hairpin[0]] != hairpin[-1]:
        # not known terminal pair, nothing to close "hairpin"
        raise RuntimeError()

    d_g = 0.0
    if hairpin in DNA_TRI_TETRA_LOOPS:
        # it's a pre-known hairpin with known value
        d_h, d_s = DNA_TRI_TETRA_LOOPS[hairpin]
        d_g = _d_g(d_h, d_s, temp)

    # add penalty based on size
    if hairpin_len in DNA_HAIRPIN_LOOPS:
        d_h, d_s = DNA_HAIRPIN_LOOPS[hairpin_len]
        d_g += _d_g(d_h, d_s, temp)
    else:
        # it's too large, extrapolate
        hairpin_max_len = max(DNA_HAIRPIN_LOOPS.keys())
        d_h, d_s = DNA_HAIRPIN_LOOPS[hairpin_max_len]
        d_g_inc = _d_g(d_h, d_s, temp)
        d_g += _j_s(hairpin_len, hairpin_max_len, d_g_inc, temp)

    # add penalty for a terminal mismatch
    if hairpin_len > 3 and pair in DNA_TERMINAL_MM:
        if pair in DNA_TERMINAL_MM:
            d_h, d_s = DNA_TERMINAL_MM[pair]
            d_g += _d_g(d_h, d_s, temp)

    # add penalty if length 3 and AT closing, formula 8 from SantaLucia, 2004
    if hairpin_len == 3 and (pair.startswith("A") or pair.startswith("T")):
        d_g += 0.5  # convert to entropy

    return d_g


def _bulge(pair: str, seq: str, i: int, j: int, bulge: str, temp: float) -> float:
    """Calculate the free energy associated with a bulge.

    Args:
        pair: The NN pair outside the bulge
        seq: The full folding DNA sequence
        i: The start index of the bulge
        j: The end index of the bulge
        loop: The sequence of the bulge
        temp: Temperature in Kelvin

    Returns:
        float: The increment in free energy from the bulge
    """

    loop_len = len(bulge) - 2  # bulge seq includes edges
    if loop_len <= 0:
        return math.inf

    # add penalty based on size
    if loop_len in DNA_BULGE_LOOPS:
        d_h, d_s = DNA_BULGE_LOOPS[loop_len]
        d_g = _d_g(d_h, d_s, temp)
    else:
        # it's too large for pre-calculated list, extrapolate
        loop_len_max = max(DNA_BULGE_LOOPS.keys())
        d_h, d_s = DNA_BULGE_LOOPS[loop_len_max]
        d_g = _d_g(d_h, d_s, temp)
        d_g = _j_s(loop_len, loop_len_max, d_g, temp)

    if loop_len == 1:
        # if len 1, include the delta G of intervening NN (SantaLucia 2004)
        d_g += _pair(pair, seq, i, j, temp)

    # penalize AT terminal bonds
    if pair.count("A"):
        d_g += 0.5

    return d_g


def _internal_loop(
    seq: str, i: int, j: int, left: str, right: str, temp: float
) -> float:
    """Calculate the free energy of an internal loop.

    The first and last bp of both left and right sequences
    are not themselves parts of the loop, but are the terminal
    bp on either side of it. They are needed for when there's
    a single internal looping bp (where just the mismatching
    free energies are used)

    Note that both left and right sequences are in 5' to 3' direction

    This is adapted from the "Internal Loops" section of SantaLucia/Hicks, 2004

    Args:
        left: The sequence on the left side
        right: The sequence on the right side
        temp: Temperature in Kelvin

    Returns:
        float: The free energy associated with the internal loop
    """

    pair_left_mm = left[:2] + "/" + right[-2::][::-1]
    pair_right_mm = left[-2:] + "/" + right[:2][::-1]
    loop_left = len(left) - 2
    loop_right = len(right) - 2
    loop_len = loop_left + loop_right

    # single bp mismatch, sum up the two single mismatch pairs
    if loop_left == 1 and loop_right == 1:
        if pair_left_mm in DNA_NN or pair_right_mm in DNA_NN:
            raise RuntimeError()

        return _pair(pair_left_mm, seq, i, j, temp) + _pair(
            pair_right_mm, seq, i + 1, j - 1, temp
        )

    # apply a penalty based on loop size
    if loop_len in DNA_INTERNAL_LOOPS:
        d_h, d_s = DNA_INTERNAL_LOOPS[loop_len]
        d_g = _d_g(d_h, d_s, temp)
    else:
        # it's too large an internal loop, extrapolate
        loop_max_len = max(DNA_INTERNAL_LOOPS.keys())
        d_h, d_s = DNA_INTERNAL_LOOPS[loop_max_len]
        d_g = _d_g(d_h, d_s, temp)
        d_g = _j_s(loop_len, loop_max_len, d_g, temp)

    # apply an asymmetry penalty
    loop_asymmetry = abs(loop_left - loop_right)
    d_g += 0.3 * loop_asymmetry

    # apply penalty based on the mismatching pairs on either side the loop
    d_h, d_s = DNA_TERMINAL_MM[pair_left_mm]
    d_g += _d_g(d_h, d_s, temp)

    d_h, d_s = DNA_TERMINAL_MM[pair_right_mm]
    d_g += _d_g(d_h, d_s, temp)

    return d_g


def _multi_branch(
    seq: str,
    i: int,
    k: int,
    j: int,
    temp: float,
    v_cache: Cache,
    w_cache: Cache,
    helix: bool = False,
) -> Tuple[float, int, int]:
    """Calculate a multi-branch energy penalty using linear formula.

    From Jaeger, Turner, and Zuker, 1989.
    Found to be better than logarithmic in Ward, et al. 2017

    Args:
        seq: The sequence being folded
        i: The left starting index
        k: The mid-point in the search
        j: The right ending index
        temp: Folding temp
        v_cache: Cache of energies where V(i,j) bond
        w_cache: Cache of min energy of substructures between W(i,j)
        helix: Whether this multibranch is enclosed by another helix

    Returns:
        float: The energy of the multi-branch penalty
    """

    e_left, e_ltype, e_lpairs = _w(seq, i, k, temp, v_cache, w_cache)
    e_right, e3_rtype, e_rpairs = _w(seq, k + 1, j, temp, v_cache, w_cache)

    # at least three multi-loops here; Fig 2A
    helixes = 3 if helix else 2
    if "BIFURCATION" in e_ltype:  # TODO: recurse to go higher
        helixes += 1
    if "BIFURCATION" in e3_rtype:
        helixes += 1

    # add up the unpaired bp count
    unpaired = 0
    i2_1 = 0
    if e_lpairs:
        i2, i2_1 = e_lpairs[0]
        unpaired += i2 - i - 1

    if e_rpairs:
        i3, i3_1 = e_rpairs[0]
        unpaired += i3 - k - 2
        unpaired += j - i3_1 - 1
        if i2_1:
            unpaired += i3 - k - 2

    # penalty for unmatched bp and multi-branch
    e_multibranch = 4.6 + 0.4 * unpaired + 0.1 * helixes

    # energy of min-energy neighbors
    e = e_multibranch + e_left + e_right

    return e, unpaired, helixes


def _traceback(
    i: int, j: int, v_cache: Cache, w_cache: Cache
) -> List[Tuple[int, int, str, float]]:
    """Traceback thru the V(i,j) and W(i,j) caches to find the structure

    For each step, get to the lowest energy W(i,j) within that block
    Store the structure in W(i,j)
    Inc i and j
    If the next structure is viable according to V(i,j), store as well
    Repeat

    Args:
        i: The leftmost index to start searching in
        j: The rightmost index to start searching in
        v_cache: Energies/structures where i and j bond
        w_cache: Energies/sub-structures between or with i and j

    Returns:
        A list of pairs (i, j) and the type of structure enclosed
    """

    # move i,j down-left to start coordinates
    _, desc, ij = w_cache[i][j]
    if "HAIRPIN" not in desc:
        while w_cache[i + 1][j][2] == ij:
            i += 1
        while w_cache[i][j - 1][2] == ij:
            j -= 1

    structs: List[Tuple[int, int, str, float]] = []
    while True:
        e, desc, ij = v_cache[i][j]

        # it's a hairpin, end of structure
        if not ij:
            # set the energy of everything relative to the hairpin
            structs.append((i, j, desc, e))
            return _trackback_energy(structs)

        # it's a stack, bulge, etc
        # there's another single structure beyond this
        if len(ij) == 1:
            structs.append((i, j, desc, e))
            i, j = ij[0]
            continue

        # it's a bifurcation
        (i1, j1), (i2, j2) = ij
        structs.append((i, j, desc, e))

        # next structure might not exist; Figure 2A
        while v_cache[i1][j1][0] == math.inf:
            i1 += 1

        while v_cache[i2][j2][0] == math.inf:
            i2 += 1

        structs = _trackback_energy(structs)
        traceback_left = _traceback(i1, j1, v_cache, w_cache)
        traceback_right = _traceback(i2, j2, v_cache, w_cache)

        e_sum = 0.0
        if traceback_left:
            left_i, left_j, _, _ = traceback_left[0]
            e_sum += w_cache[left_i][left_j][0]
        if traceback_right:
            right_i, right_j, _, _ = traceback_right[0]
            e_sum += w_cache[right_i][right_j][0]

        if structs:
            last_i, last_j, last_desc, last_e = structs[-1]
            structs[-1] = (last_i, last_j, last_desc, round(last_e - e_sum, 2))

        return structs + traceback_left + traceback_right

    return _trackback_energy(structs)


def _trackback_energy(
    structs: List[Tuple[int, int, str, float]]
) -> List[Tuple[int, int, str, float]]:
    """Add energy to each structure, based on how it's W(i,j) differs from the one after
    
    Args:
        structs: The structures for whom energy is being calculated
    
    Returns:
        List[Tuple[int, int, str, float]]: Structures in the folded DNA with energy
    """

    structs_e: List[Tuple[int, int, str, float]] = []
    for index, struct in enumerate(structs):
        i, j, desc, e = struct
        e_next = 0.0
        if index < len(structs) - 1:
            e_next = structs[index + 1][3]
        structs_e.append((i, j, desc, round(e - e_next, 2)))
    return structs_e


def _debug(v_cache, w_cache):
    """Temporary _debug function for logging energies

    Args:
        v_cache: V(i, j)
        w_cache: W(i,j)
    """

    print("\n")
    print(",".join([str(n) for n in range(len(w_cache))]))
    for i, row in enumerate(w_cache):
        print(
            ",".join([str(round(r, 1)) if r is not None else "." for r, _, _ in row])
            + ","
            + str(i)
        )

    print("\n")
    print(",".join([str(n) for n in range(len(v_cache))]))
    for i, row in enumerate(v_cache):
        print(
            ",".join([str(round(r, 1)) if r is not None else "." for r, _, _ in row])
            + ","
            + str(i)
        )

    print("\n")
    print(",".join([str(n) for n in range(len(w_cache))]))
    for i, row in enumerate(w_cache):
        print(",".join([str(ij).replace(", ", "-") for _, _, ij in row]) + "," + str(i))

    print("\n")
    print(",".join([str(n) for n in range(len(v_cache))]))
    for i, row in enumerate(v_cache):
        print(",".join([str(ij).replace(", ", "-") for _, _, ij in row]) + "," + str(i))

    print("\n")
    print(",".join([str(n) for n in range(len(w_cache))]))
    for i, row in enumerate(w_cache):
        print(",".join([t if t else "." for _, t, _ in row]) + "," + str(i))

    print("\n")
    print(",".join([str(n) for n in range(len(v_cache))]))
    for i, row in enumerate(v_cache):
        print(",".join([t if t else "." for _, t, _ in row]) + "," + str(i))


def parse_args(args):
    """Parse command line parameters

    Created and based on example from pyscaffold:
    https://github.com/pyscaffold/pyscaffold/blob/master/src/pyscaffold/templates/skeleton.template
    
    Args:
        args ([str]): List of parameters as strings
    
    Returns:
        :obj:`argparse.Namespace`: command line parameters namespace
    """

    parser = argparse.ArgumentParser(
        description="Predict minimum free energy structure of a nucleic acid sequence"
    )

    parser.add_argument(
        "seq", type=str, metavar="SEQ", help="nucleic acid sequence to fold",
    )
    parser.add_argument(
        "-t",
        dest="temp",
        type=float,
        metavar="FLOAT",
        default=37.0,
        help="temperature to fold at (Celcius)",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="seqfold {ver}".format(ver=__version__),
    )

    return parser.parse_args(args)


def run():
    """Entry point for console_scripts.

    Simply fold the DNA and report the minimum free energy
    """

    parsed_args = parse_args(sys.argv[1:])
    dg = calc_dg(parsed_args.seq, temp=parsed_args.temp)
    print(dg)


if __name__ == "__main__":
    run()
