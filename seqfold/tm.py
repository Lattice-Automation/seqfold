"""Calculate the tm of a DNA sequence"""

import math

from .dna import DNA_COMPLEMENT, DNA_ENERGIES


def calc_tm(seq1: str, seq2: str = "", pcr: bool = True) -> float:
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
    delta_h, delta_s = DNA_ENERGIES.NN["init"]

    # add in initial A/T and initial G/Cs
    init = seq1[0] + seq1[-1]
    init_at = init.count("A") + init.count("T")
    init_gc = init.count("G") + init.count("C")
    init_at_h, init_at_s = DNA_ENERGIES.NN["init_A/T"]
    init_gc_h, init_gc_s = DNA_ENERGIES.NN["init_G/C"]
    delta_h += init_at * init_at_h + init_gc * init_gc_h
    delta_s += init_at * init_at_s + init_gc * init_gc_s

    # work through each nearest neighbor pair
    for i in range(len(seq1) - 1):
        pair1 = seq1[i] + seq1[i + 1]
        pair2 = seq2[i] + seq2[i + 1]
        pair = pair1 + "/" + pair2

        # assuming internal neighbor pair
        pair_delta_h, pair_delta_s = 0.0, 0.0
        if pair in DNA_ENERGIES.NN:
            pair_delta_h, pair_delta_s = DNA_ENERGIES.NN[pair]
        elif pair in DNA_ENERGIES.INTERNAL_MM:
            pair_delta_h, pair_delta_s = DNA_ENERGIES.INTERNAL_MM[pair]

        # overwrite if it's a terminal pair
        if i in (0, len(seq1) - 2):
            if pair in DNA_ENERGIES.TERMINAL_MM:
                pair_delta_h, pair_delta_s = DNA_ENERGIES.TERMINAL_MM[pair]

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


def _gc(seq: str) -> float:
    """Return the GC ratio of a sequence."""

    seq = seq.upper()
    return float(seq.count("G") + seq.count("C")) / float(len(seq))
