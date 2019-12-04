"""Predict nucleic acid secondary structure"""

import argparse
import math
import sys
from typing import Any, Dict, List, Tuple

from . import __version__
from .dna import DNA_COMPLEMENT, DNA_ENERGIES
from .rna import RNA_ENERGIES


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
    delta_h, delta_s = DNA_ENERGIES["NN"]["init"]

    # add in initial A/T and initial G/Cs
    init = seq1[0] + seq1[-1]
    init_at = init.count("A") + init.count("T")
    init_gc = init.count("G") + init.count("C")
    init_at_h, init_at_s = DNA_ENERGIES["NN"]["init_A/T"]
    init_gc_h, init_gc_s = DNA_ENERGIES["NN"]["init_G/C"]
    delta_h += init_at * init_at_h + init_gc * init_gc_h
    delta_s += init_at * init_at_s + init_gc * init_gc_s

    # work through each nearest neighbor pair
    for i in range(len(seq1) - 1):
        pair1 = seq1[i] + seq1[i + 1]
        pair2 = seq2[i] + seq2[i + 1]
        pair = pair1 + "/" + pair2

        # assuming internal neighbor pair
        pair_delta_h, pair_delta_s = 0.0, 0.0
        if pair in DNA_ENERGIES["NN"]:
            pair_delta_h, pair_delta_s = DNA_ENERGIES["NN"][pair]
        elif pair in DNA_ENERGIES["INTERNAL_MM"]:
            pair_delta_h, pair_delta_s = DNA_ENERGIES["INTERNAL_MM"][pair]

        # overwrite if it's a terminal pair
        if i in (0, len(seq1) - 2):
            if pair in DNA_ENERGIES["TERMINAL_MM"]:
                pair_delta_h, pair_delta_s = DNA_ENERGIES["TERMINAL_MM"][pair]

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


def calc_dg(seq: str, temp: float, log_structs: bool = False) -> float:
    """Fold the sequence and return just the delta G of the structure
    
    Args:
        seq: The sequence to fold
        temp: The temperature to fold at

    Optional Args:
        log_structs: Whether to log the structure as well
    
    Returns:
        float: The minimum free energy of the folded sequence
    """

    structs = fold(seq, temp)

    if log_structs:
        for (i, j, desc, ddg) in structs:
            print(i, j, desc, ddg)

    return round(sum([s[-1] for s in structs]), 2)


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

    # figure out whether it's DNA or RNA, choose energy map
    dna = True
    bps = set(seq)
    if not bps.difference("ATGU"):
        dna = False
    elif set("ATGC").difference(bps):
        diff = str(set("ATGC").difference(bps))
        raise RuntimeError(f"Unknown bp: ${diff}.\n\tOnly DNA/RNA foldable")
    emap = DNA_ENERGIES if dna else RNA_ENERGIES

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
            _v(seq, i, i + f_len, temp, v_cache, w_cache, emap)

    # gather the min energy structure over the full sequence
    min_e, _, _ = _w(seq, 0, len(seq) - 1, temp, v_cache, w_cache, emap)
    min_e = round(min_e, 2)

    # get the structure out of the cache
    # _debug(v_cache, w_cache)
    structs = _traceback(0, len(seq) - 1, v_cache, w_cache)
    total_e = sum(s[-1] for s in structs)
    assert abs(total_e - min_e) < 0.2, f"{total_e} != {min_e}"

    return structs


def _v(
    seq: str, i: int, j: int, temp: float, v_cache: Cache, w_cache: Cache, emap: Any
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
        emap: Energy map for DNA/RNA

    Returns:
        float: The minimum energy folding structure possible between i and j on seq
    """

    if v_cache[i][j][0] is not None:
        return v_cache[i][j]

    # the ends must basepair for V(i,j)
    if emap["COMPLEMENT"][seq[i]] != seq[j]:
        v_cache[i][j] = (math.inf, "", [])
        return v_cache[i][j]

    # E1 = FH(i, j); hairpin
    pair = seq[i] + seq[i + 1] + "/" + seq[j] + seq[j - 1]
    e1, e1_type = _hairpin(seq, i, j, temp, emap), "HAIRPIN:" + pair
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
            pair_outer = pair_left in emap["NN"] or pair_right in emap["NN"]

            stack = i_1 == i + 1 and j_1 == j - 1 and pair in emap["NN"]
            bulge_left = i_1 > i + 1 and pair in emap["NN"]
            bulge_right = j_1 < j - 1 and pair in emap["NN"]

            loop_left = seq[i : i_1 + 1]
            loop_right = seq[j_1 : j + 1]

            e2_test, e2_test_type = math.inf, ""
            if stack:
                # it's a neighboring/stacking pair in a helix
                e2_test = _pair(pair, seq, i, j, temp, emap)
                e2_test_type = "STACK:" + pair

                if i > 0 and j == len(seq) - 1 or i == 0 and j < len(seq) - 1:
                    # there's a dangling end
                    e2_test_type = "STACK_DE:" + pair
            elif bulge_left and bulge_right and not pair_outer:
                # it's an interior loop
                e2_test = _internal_loop(seq, i, j, loop_left, loop_right, temp, emap)
                e2_test_type = "INTERIOR_LOOP"

                if len(loop_left) == 3 and len(loop_right) == 3:
                    # technically an interior loop of 1. really 1bp mismatch
                    e2_test_type = "STACK:" + loop_left + "/" + loop_right[::-1]
            elif bulge_left and not bulge_right:
                # it's a bulge on the left side
                e2_test = _bulge(pair, seq, i, j, loop_left, temp, emap)
                e2_test_type = "BULGE:" + str(i) + "-" + str(i_1)
            elif bulge_right and not bulge_left:
                # it's a bulge on the right side
                e2_test = _bulge(pair, seq, i, j, loop_right, temp, emap)
                e2_test_type = "BULGE:" + str(j_1) + "-" + str(j)
            else:
                # it's basically a hairpin, only outside bp match
                continue

            # add V(i', j')
            e2_test += _v(seq, i_1, j_1, temp, v_cache, w_cache, emap)[0]
            if e2_test < e2:
                e2, e2_type, e2_ij = e2_test, e2_test_type, [(i_1, j_1)]

    # E3 = min{W(i+1,i') + W(i'+1,j-1)}, i+1<i'<j-2
    e3, e3_type = math.inf, "BIFURCATION"
    e3_ij: List[Tuple[int, int]] = []
    for i_1 in range(i + 2, j - 2):
        e3_test, unpaired, helixes = _multi_branch(
            seq, i + 1, i_1, j - 1, temp, v_cache, w_cache, emap, True
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
    seq: str, i: int, j: int, temp: float, v_cache: Cache, w_cache: Cache, emap: Any
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

    w1 = _w(seq, i + 1, j, temp, v_cache, w_cache, emap)
    w2 = _w(seq, i, j - 1, temp, v_cache, w_cache, emap)
    w3 = _v(seq, i, j, temp, v_cache, w_cache, emap)

    w4, w4_type = math.inf, "BIFURCATION"
    w4_ij: List[Tuple[int, int]] = []
    for i_1 in range(i + 1, j - 1):
        w4_test, unpaired, helixes = _multi_branch(
            seq, i, i_1, j, temp, v_cache, w_cache, emap, False,
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


def _pair(
    pair: str, seq: str, i: int, j: int, temp: float, emap: Dict[str, Any]
) -> float:
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
        d_h, d_s = emap["NN"][pair] if pair in emap["NN"] else emap["INTERNAL_MM"][pair]
        return _d_g(d_h, d_s, temp)

    if i == 0 and j == len(seq) - 1:
        # it's terminal
        d_h, d_s = emap["NN"][pair] if pair in emap["NN"] else emap["TERMINAL_MM"][pair]
        return _d_g(d_h, d_s, temp)

    if i > 0 and j == len(seq) - 1:
        # it's dangling on left
        d_h, d_s = emap["NN"][pair] if pair in emap["NN"] else emap["TERMINAL_MM"][pair]
        d_g = _d_g(d_h, d_s, temp)

        pair_de = seq[i - 1] + seq[i] + "/." + seq[j]
        d_h, d_s = emap["DE"][pair_de]
        return d_g + _d_g(d_h, d_s, temp)

    if i == 0 and j < len(seq) - 1:
        # it's dangling on right
        d_h, d_s = emap["NN"][pair] if pair in emap["NN"] else emap["TERMINAL_MM"][pair]
        d_g = _d_g(d_h, d_s, temp)

        pair_de = "." + seq[i] + "/" + seq[j + 1] + seq[j]
        d_h, d_s = emap["DE"][pair_de]
        return d_g + _d_g(d_h, d_s, temp)

    raise RuntimeError


def _hairpin(seq: str, i: int, j: int, temp: float, emap: Dict[str, Any]) -> float:
    """Calculate the free energy of a hairpin.

    Args:
        seq: The sequence we're folding
        i: The index of start of hairpin
        j: The index of end of hairpin
        temp: Temperature in Kelvin
        emap: Map of energies

    Returns:
        float: The free energy increment from the hairpin structure
    """

    if j - i < 4:
        return math.inf

    hairpin = seq[i : j + 1]
    hairpin_len = len(hairpin) - 2
    pair = hairpin[0] + hairpin[1] + "/" + hairpin[-1] + hairpin[-2]

    if emap["COMPLEMENT"][hairpin[0]] != hairpin[-1]:
        # not known terminal pair, nothing to close "hairpin"
        raise RuntimeError()

    d_g = 0.0
    if "TRI_TETRA_LOOPS" in emap and hairpin in emap["TRI_TETRA_LOOPS"]:
        # it's a pre-known hairpin with known value
        d_h, d_s = emap["TRI_TETRA_LOOPS"][hairpin]
        d_g = _d_g(d_h, d_s, temp)

    # add penalty based on size
    if hairpin_len in emap["HAIRPIN_LOOPS"]:
        d_h, d_s = emap["HAIRPIN_LOOPS"][hairpin_len]
        d_g += _d_g(d_h, d_s, temp)
    else:
        # it's too large, extrapolate
        hairpin_max_len = max(emap["HAIRPIN_LOOPS"].keys())
        d_h, d_s = emap["HAIRPIN_LOOPS"][hairpin_max_len]
        d_g_inc = _d_g(d_h, d_s, temp)
        d_g += _j_s(hairpin_len, hairpin_max_len, d_g_inc, temp)

    # add penalty for a terminal mismatch
    if hairpin_len > 3 and pair in emap["TERMINAL_MM"]:
        if pair in emap["TERMINAL_MM"]:
            d_h, d_s = emap["TERMINAL_MM"][pair]
            d_g += _d_g(d_h, d_s, temp)

    # add penalty if length 3 and AT closing, formula 8 from SantaLucia, 2004
    if hairpin_len == 3 and (pair.startswith("A") or pair.startswith("T")):
        d_g += 0.5  # convert to entropy

    return d_g


def _bulge(
    pair: str, seq: str, i: int, j: int, bulge: str, temp: float, emap: Dict[str, Any]
) -> float:
    """Calculate the free energy associated with a bulge.

    Args:
        pair: The NN pair outside the bulge
        seq: The full folding DNA sequence
        i: The start index of the bulge
        j: The end index of the bulge
        loop: The sequence of the bulge
        temp: Temperature in Kelvin
        emap: Map to DNA/RNA energies

    Returns:
        float: The increment in free energy from the bulge
    """

    loop_len = len(bulge) - 2  # bulge seq includes edges
    if loop_len <= 0:
        return math.inf

    # add penalty based on size
    if loop_len in emap["BULGE_LOOPS"]:
        d_h, d_s = emap["BULGE_LOOPS"][loop_len]
        d_g = _d_g(d_h, d_s, temp)
    else:
        # it's too large for pre-calculated list, extrapolate
        loop_len_max = max(emap["BULGE_LOOPS"].keys())
        d_h, d_s = emap["BULGE_LOOPS"][loop_len_max]
        d_g = _d_g(d_h, d_s, temp)
        d_g = _j_s(loop_len, loop_len_max, d_g, temp)

    if loop_len == 1:
        # if len 1, include the delta G of intervening NN (SantaLucia 2004)
        d_g += _pair(pair, seq, i, j, temp, emap)

    # penalize AT terminal bonds
    if pair.count("A"):
        d_g += 0.5

    return d_g


def _internal_loop(
    seq: str, i: int, j: int, left: str, right: str, temp: float, emap: Dict[str, Any]
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
        seq: The sequence we're folding
        i: The index of the start of structure on left side
        j: The index of the end of structure on right side
        left: The sequence on the left side
        right: The sequence on the right side
        temp: Temperature in Kelvin
        emap: Dictionary mapping to energies for DNA/RNA

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
        if pair_left_mm in emap["NN"] or pair_right_mm in emap["NN"]:
            raise RuntimeError()

        return _pair(pair_left_mm, seq, i, j, temp, emap) + _pair(
            pair_right_mm, seq, i + 1, j - 1, temp, emap
        )

    # apply a penalty based on loop size
    if loop_len in emap["INTERNAL_LOOPS"]:
        d_h, d_s = emap["INTERNAL_LOOPS"][loop_len]
        d_g = _d_g(d_h, d_s, temp)
    else:
        # it's too large an internal loop, extrapolate
        loop_max_len = max(emap["INTERNAL_LOOPS"].keys())
        d_h, d_s = emap["INTERNAL_LOOPS"][loop_max_len]
        d_g = _d_g(d_h, d_s, temp)
        d_g = _j_s(loop_len, loop_max_len, d_g, temp)

    # apply an asymmetry penalty
    loop_asymmetry = abs(loop_left - loop_right)
    d_g += 0.3 * loop_asymmetry

    # apply penalty based on the mismatching pairs on either side the loop
    d_h, d_s = emap["TERMINAL_MM"][pair_left_mm]
    d_g += _d_g(d_h, d_s, temp)

    d_h, d_s = emap["TERMINAL_MM"][pair_right_mm]
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
    emap: Dict[str, Any],
    helix: bool = False,
) -> Tuple[float, int, int]:
    """Calculate a multi-branch energy penalty using a linear formula.

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
        emap: Map to DNA/RNA energies

    Returns:
        float: The energy of the multi-branch penalty
    """

    e_left, e_ltype, e_lpairs = _w(seq, i, k, temp, v_cache, w_cache, emap)
    e_right, e3_rtype, e_rpairs = _w(seq, k + 1, j, temp, v_cache, w_cache, emap)

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
    a, b, c = emap["MULTIBRANCH"]
    e_multibranch = a + b * unpaired + c * helixes

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
        description="Predict the minimum free energy (kcal/mol) of a nucleic acid sequence"
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
        "--verbose",
        action="store_true",
        help="log each structure (i,j,description,ddG)",
    )
    parser.add_argument(
        "--version", action="version", version="seqfold {ver}".format(ver=__version__),
    )

    return parser.parse_args(args)


def run():
    """Entry point for console_scripts.

    Simply fold the DNA and report the minimum free energy
    """

    args = parse_args(sys.argv[1:])
    dg = calc_dg(args.seq, temp=args.temp, log_structs=args.verbose)
    print(dg)


if __name__ == "__main__":
    run()
