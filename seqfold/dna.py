"""DNA enthalpy and entropy change parameters."""

from typing import Any, Dict, Tuple

from .types import Comp, MultiBranch, BpEnergy, LoopEnergy, Energies


DNA_COMPLEMENT: Comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

DNA_MULTIBRANCH: MultiBranch = (2.6, 0.2, 0.2, 2.0)
"""a, b, c, d in a linear multi-branch energy change function.

Inferred from:
Supplemental Material: Annu.Rev.Biophs.Biomol.Struct.33:415-40
doi: 10.1146/annurev.biophys.32.110601.141800
The Termodynamicso f DNA Structural Motifs
SantaLucia and Hicks, 2004
"""

DNA_NN: BpEnergy = {
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

DNA_INTERNAL_MM: BpEnergy = {
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

DNA_TERMINAL_MM: BpEnergy = {
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

DNA_DE: BpEnergy = {
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

DNA_TRI_TETRA_LOOPS: BpEnergy = {
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

DNA_INTERNAL_LOOPS: LoopEnergy = {
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

DNA_BULGE_LOOPS: LoopEnergy = {
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

DNA_HAIRPIN_LOOPS: LoopEnergy = {
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


DNA_ENERGIES = Energies(
    DNA_BULGE_LOOPS,
    DNA_COMPLEMENT,
    DNA_DE,
    DNA_HAIRPIN_LOOPS,
    DNA_MULTIBRANCH,
    DNA_INTERNAL_LOOPS,
    DNA_INTERNAL_MM,
    DNA_NN,
    DNA_TERMINAL_MM,
    DNA_TRI_TETRA_LOOPS,
)

