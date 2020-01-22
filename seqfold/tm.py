"""Calculate the tm of a DNA sequence"""

import math
from typing import List, Tuple

from .dna import DNA_COMPLEMENT, DNA_ENERGIES

TmCache = List[List[float]]


def tm(seq1: str, seq2: str = "", pcr: bool = True) -> float:
    """Calculate the annealing temperature between seq1 and seq2.

    If seq2 is not provided, its exact complement is used.
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

    seq1, seq2 = _parse_input(seq1, seq2)

    # sum enthalpy and entropy. Enthaply is first value of each tuple and
    # entropy is the second value of each tuple in:
    # SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct 33: 415-440

    # start with initiation enthalpy and entropy
    dh, ds = DNA_ENERGIES.NN["init"]

    # add in initial A/T and initial G/Cs
    init = seq1[0] + seq1[-1]
    init_at = init.count("A") + init.count("T")
    init_gc = init.count("G") + init.count("C")
    init_at_h, init_at_s = DNA_ENERGIES.NN["init_A/T"]
    init_gc_h, init_gc_s = DNA_ENERGIES.NN["init_G/C"]
    dh += init_at * init_at_h + init_gc * init_gc_h
    ds += init_at * init_at_s + init_gc * init_gc_s

    # work through each nearest neighbor pair
    for i in range(len(seq1) - 1):
        pair = f"{seq1[i]}{seq1[i + 1]}/{seq2[i]}{seq2[i + 1]}"

        # assuming internal neighbor pair
        pair_dh, pair_ds = 0.0, 0.0
        if pair in DNA_ENERGIES.NN:
            pair_dh, pair_ds = DNA_ENERGIES.NN[pair]
        elif pair in DNA_ENERGIES.INTERNAL_MM:
            pair_dh, pair_ds = DNA_ENERGIES.INTERNAL_MM[pair]

        # overwrite if it's a terminal pair
        if i in (0, len(seq1) - 2):
            if pair in DNA_ENERGIES.TERMINAL_MM:
                pair_dh, pair_ds = DNA_ENERGIES.TERMINAL_MM[pair]

        dh += pair_dh
        ds += pair_ds

    gc = _gc(seq1)
    return _calc_tm(dh, ds, pcr, gc, len(seq1))


def tm_cache(seq1: str, seq2: str = "", pcr: bool = True) -> TmCache:
    """Return a TmCache where each (i, j) returns the Tm for that subspan.

    1. Build up the 2D matrixes for the tm calculation:
        - dh
        - ds
        - gc count
    2. Fill in each cell, (i, j), with the estimated tm for that range from i to j

    Args:
        seq1: The seq whose tm is calculated

    Keyword Args:
        seq2: The seq that seq1 anneals to in 3' -> 5' direction
        pcr: Whether tm is being calculated for the oligo is in a
            PCR reaction mixture. If so, ion and Tris concentrations
            that match a typical NEB/ThermoFisher PCR mixture are used

    Returns:
        TmCache: Where querying the cache with (i, j) returns the tm of the
            subsequence starting with i and ending with j, inclusive
    """

    seq1, seq2 = _parse_input(seq1, seq2)
    n = len(seq1)  # using nearest neighbors, -1

    arr_dh, arr_ds, arr_gc, arr_tm = [], [], [], []
    for _ in range(n):
        arr_dh.append([0.0] * n)
        arr_ds.append([0.0] * n)
        arr_gc.append([0.0] * n)
        arr_tm.append([math.inf] * n)

    # fill in the diagonal
    for i in range(n):
        if i == n - 1:  # hackish
            arr_dh[i][i] = arr_dh[i - 1][i - 1]
            arr_ds[i][i] = arr_ds[i - 1][i - 1]
            arr_gc[i][i] = arr_gc[i - 1][i - 1]
            continue

        pair = f"{seq1[i]}{seq1[i + 1]}/{seq2[i]}{seq2[i + 1]}"
        dh, ds = (
            DNA_ENERGIES.NN[pair]
            if pair in DNA_ENERGIES.NN
            else DNA_ENERGIES.INTERNAL_MM[pair]
        )

        arr_dh[i][i] = dh
        arr_ds[i][i] = ds
        arr_gc[i][i] = 1.0 if seq1[i] in "GC" else 0.0

        if i == n - 2 and not arr_gc[i][i]:  # don't ignore last pair
            arr_gc[i][i] = 1.0 if seq1[i + 1] in "GC" else 0.0

    # fill in the tm array
    for i in range(n):
        for j in range(i + 1, n):
            arr_dh[i][j] = arr_dh[i][j - 1] + arr_dh[j][j]
            arr_ds[i][j] = arr_ds[i][j - 1] + arr_ds[j][j]
            arr_gc[i][j] = arr_gc[i][j - 1] + arr_gc[j][j]
            arr_tm[i][j] = _calc_tm(
                arr_dh[i][j], arr_ds[i][j], pcr, arr_gc[i][j] / (j - i + 1), j - i + 1
            )

    return arr_tm


def _parse_input(seq1: str, seq2: str = "") -> Tuple[str, str]:
    """Parse and prepare the input sequences. Throw if there's an issue.

    Args:
        seq1: The main sequence whose tm is being calculated
    
    Keyword Args:
        seq2: The second sequence that seq2 is annealing to
    
    Returns:
        (str, str): The sequences to use for tm calculation
    """

    if "SeqRecord" in str(type(seq1)):
        seq1 = str(seq1.seq)  # type: ignore

    if "SeqRecord" in str(type(seq2)):
        seq2 = str(seq2.seq)  # type: ignore

    seq1 = seq1.upper()
    if not seq2:
        seq2 = "".join([DNA_COMPLEMENT[c] for c in seq1])

    if len(seq1) != len(seq2):
        raise ValueError(
            f"Length mismatch between seq1 {len(seq1)} and seq2 {len(seq2)}"
        )

    if len(seq1) < 2:
        raise ValueError(f"Sequence, {len(seq1)}bp, is too short for tm calculation")

    return seq1, seq2


def _calc_tm(dh: float, ds: float, pcr: bool, gc: float, seq_len: int) -> float:
    """Apply the correction formula to estimate Tm

    Args:
        dh: Accumulated entropy
        ds: Accumulated enthalpy
        pcr: Whether this is for PCR or not
        gc: The GC% of the sequence
        seq_len: The length of the sequence

    Returns:
        float: The estimated tm
    """

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
            return (4.29 * gc / 100 - 3.95) * 1e-5 * math.log(mon) + 9.4e-6 * math.log(
                mon
            ) ** 2
        elif R < 6.0:
            a = 3.92 * (0.843 - 0.352 * math.sqrt(mon) * math.log(mon))
            d = 1.42 * (1.279 - 4.03e-3 * math.log(mon) - 8.03e-3 * math.log(mon) ** 2)
            g = 8.31 * (0.486 - 0.258 * math.log(mon) + 5.25e-3 * math.log(mon) ** 3)
    corr = (
        a
        + b * math.log(mg)
        + (gc / 100) * (c + d * math.log(mg))
        + (1 / (2.0 * (seq_len - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
    ) * 1e-5

    # tm with concentration consideration
    k = (seq1_conc - (seq2_conc / 2.0)) * 1e-9
    R = 1.9872
    est = (dh * 1000.0) / (ds + R * math.log(k)) - 273.15

    # add in salt correction
    est = 1 / (1 / (est + 273.15) + corr) - 273.1

    return round(est, 1)


def _gc(seq: str) -> float:
    """Return the GC ratio of a sequence."""

    return float(seq.count("G") + seq.count("C")) / float(len(seq))
