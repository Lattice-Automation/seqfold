"""Predict nucleic acid secondary structure"""

import math
from typing import Any, Dict, Optional, List, Tuple

from .dna import DNA_ENERGIES
from .rna import RNA_ENERGIES
from .types import Energies


class Struct:
    """A single structure with a free energy, description, and inward children."""

    fmt = "{:>4} {:>4} {:>6}  {:<15}"

    def __init__(
        self, e: float = -math.inf, desc: str = "", ij: List[Tuple[int, int]] = []
    ):
        self.e: float = e
        self.desc: str = desc
        self.ij: List[Tuple[int, int]] = ij

    def __eq__(self, other) -> bool:
        return self.e == other.e and self.ij == other.ij

    def __str__(self) -> str:
        i = self.ij[0][0] if self.ij else ""
        j = self.ij[0][1] if self.ij else ""
        e = str(self.e)
        return self.fmt.format(i, j, e, self.desc)

    def __bool__(self) -> bool:
        return self.e != math.inf and self.e != -math.inf

    def with_ij(self, ij: List[Tuple[int, int]]):
        return Struct(self.e, self.desc, ij)


STRUCT_DEFAULT = Struct(-math.inf)
STRUCT_NULL = Struct(math.inf)


Cache = List[List[Struct]]
"""A map from i, j tuple to a value."""


def fold_cache(seq: str, temp: float = 37.0) -> Tuple[Cache, Cache]:
    """Fold a nucleic acid sequence and return the w_cache.

    The Cache is useful for gathering many possible energies
    between a series of (i,j) combinations.

    Args:
        seq: The sequence to fold
    
    Keyword Args:
        temp: The temperature to fold at
    
    Returns:
        (Cache, Cache): The w_cache and the v_cache for traversal later
    """

    # if it's a SeqRecord, gather just the seq
    if "SeqRecord" in str(type(seq)):
        seq = str(seq.seq)  # type: ignore

    seq = seq.upper()
    temp = temp + 273.15  # kelvin

    # figure out whether it's DNA or RNA, choose energy map
    dna = True
    bps = set(seq)
    if "U" in bps and "T" in bps:
        raise RuntimeError(
            "Both T and U in sequence. Provide one or the other for DNA OR RNA."
        )
    if all(bp in "AUCG" for bp in bps):
        dna = False
    elif any(bp not in "ATGC" for bp in bps):
        diff = [bp for bp in bps if bp not in "ATUGC"]
        raise RuntimeError(f"Unknown bp: {diff}. Only DNA/RNA foldable")
    emap = DNA_ENERGIES if dna else RNA_ENERGIES

    n = len(seq)
    v_cache: Cache = []
    w_cache: Cache = []
    for _ in range(n):
        v_cache.append([STRUCT_DEFAULT] * n)
        w_cache.append([STRUCT_DEFAULT] * n)

    # fill the cache
    _w(seq, 0, n - 1, temp, v_cache, w_cache, emap)

    return v_cache, w_cache


def fold(seq: str, temp: float = 37.0) -> List[Struct]:
    """Fold the DNA sequence and return lowest free energy score.

    Based on the approach described in:
    Zuker and Stiegler, 1981
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC326673/pdf/nar00394-0137.pdf

    If the sequence is 50 or more bp long, "isolated" matching bp
    are ignored in V(i,j). This is based on an approach described in:
    Mathews, Sabina, Zuker and Turner, 1999
    https://www.ncbi.nlm.nih.gov/pubmed/10329189

    Args:
        seq: The sequence to fold

    Keyword Args:
        temp: The temperature the fold takes place in, in Celcius

    Returns:
        List[Struct]: A list of structures. Stacks, bulges, hairpins, etc.
    """

    v_cache, w_cache = fold_cache(seq, temp)
    n = len(seq)

    # get the minimum free energy structure out of the cache
    return traceback(0, n - 1, v_cache, w_cache)


def calc_dg(seq: str, temp: float = 37.0) -> float:
    """Fold the sequence and return just the delta G of the structure

    Args:
        seq: The sequence to fold
    
    Keyword Args:
        temp: The temperature to fold at

    Returns:
        float: The minimum free energy of the folded sequence
    """

    structs = fold(seq, temp)
    return round(sum(s.e for s in structs), 2)


def _w(
    seq: str,
    i: int,
    j: int,
    temp: float,
    v_cache: Cache,
    w_cache: Cache,
    emap: Energies,
) -> Struct:
    """Find and return the lowest free energy structure in Sij subsequence

    Figure 2B in Zuker and Stiegler, 1981

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

    if w_cache[i][j] != STRUCT_DEFAULT:
        return w_cache[i][j]

    if j - i < 4:
        w_cache[i][j] = STRUCT_NULL
        return w_cache[i][j]

    w1 = _w(seq, i + 1, j, temp, v_cache, w_cache, emap)
    w2 = _w(seq, i, j - 1, temp, v_cache, w_cache, emap)
    w3 = _v(seq, i, j, temp, v_cache, w_cache, emap)

    w4 = STRUCT_NULL
    for k in range(i + 1, j - 1):
        w4_test = _multi_branch(seq, i, k, j, temp, v_cache, w_cache, emap, False)

        if w4_test and w4_test.e < w4.e:
            w4 = w4_test

    w = _min_struct(w1, w2, w3, w4)
    w_cache[i][j] = w
    return w


def _v(
    seq: str,
    i: int,
    j: int,
    temp: float,
    v_cache: Cache,
    w_cache: Cache,
    emap: Energies,
) -> Struct:
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

    if v_cache[i][j] != STRUCT_DEFAULT:
        return v_cache[i][j]

    # the ends must basepair for V(i,j)
    if emap.COMPLEMENT[seq[i]] != seq[j]:
        v_cache[i][j] = STRUCT_NULL
        return v_cache[i][j]
    # if the basepair is isolated, and the seq large, penalize at 1,600 kcal/mol
    # heuristic for speeding this up
    # from https://www.ncbi.nlm.nih.gov/pubmed/10329189
    isolated_outer = True
    if i and j < len(seq) - 1:
        isolated_outer = emap.COMPLEMENT[seq[i - 1]] != seq[j + 1]
    isolated_inner = emap.COMPLEMENT[seq[i + 1]] != seq[j - 1]

    if isolated_outer and isolated_inner:
        v_cache[i][j] = Struct(1600)
        return v_cache[i][j]

    # E1 = FH(i, j); hairpin
    pair = _pair(seq, i, i + 1, j, j - 1)
    e1 = Struct(_hairpin(seq, i, j, temp, emap), "HAIRPIN:" + pair)
    if j - i == 4:  # small hairpin; 4bp
        v_cache[i][j] = e1
        w_cache[i][j] = e1
        return v_cache[i][j]

    # E2 = min{FL(i, j, i', j') + V(i', j')}, i<i'<j'<j
    # stacking region or bulge or interior loop; Figure 2A(2)
    # j-i=d>4; various pairs i',j' for j'-i'<d
    n = len(seq)
    e2 = Struct(math.inf)
    for i1 in range(i + 1, j - 4):
        for j1 in range(i1 + 4, j):
            # i1 and j1 must match
            if emap.COMPLEMENT[seq[i1]] != seq[j1]:
                continue

            pair = _pair(seq, i, i1, j, j1)
            pair_left = _pair(seq, i, i + 1, j, j - 1)
            pair_right = _pair(seq, i1 - 1, i1, j1 + 1, j1)
            pair_inner = pair_left in emap.NN or pair_right in emap.NN

            stack = i1 == i + 1 and j1 == j - 1
            bulge_left = i1 > i + 1
            bulge_right = j1 < j - 1

            e2_test, e2_test_type = math.inf, ""
            if stack:
                # it's a neighboring/stacking pair in a helix
                e2_test = _stack(seq, i, i1, j, j1, temp, emap)
                e2_test_type = f"STACK:{pair}"

                if i > 0 and j == n - 1 or i == 0 and j < n - 1:
                    # there's a dangling end
                    e2_test_type = f"STACK_DE:{pair}"
            elif bulge_left and bulge_right and not pair_inner:
                # it's an interior loop
                e2_test = _internal_loop(seq, i, i1, j, j1, temp, emap)
                e2_test_type = f"INTERIOR_LOOP:{str(i1 - i)}/{str(j - j1)}"

                if i1 - i == 2 and j - j1 == 2:
                    loop_left = seq[i : i1 + 1]
                    loop_right = seq[j1 : j + 1]
                    # technically an interior loop of 1. really 1bp mismatch
                    e2_test_type = f"STACK:{loop_left}/{loop_right[::-1]}"
            elif bulge_left and not bulge_right:
                # it's a bulge on the left side
                e2_test = _bulge(seq, i, i1, j, j1, temp, emap)
                e2_test_type = f"BULGE:{str(i1 - i)}"
            elif not bulge_left and bulge_right:
                # it's a bulge on the right side
                e2_test = _bulge(seq, i, i1, j, j1, temp, emap)
                e2_test_type = f"BULGE:{str(j - j1)}"
            else:
                # it's basically a hairpin, only outside bp match
                continue

            # add V(i', j')
            e2_test += _v(seq, i1, j1, temp, v_cache, w_cache, emap).e
            if e2_test != -math.inf and e2_test < e2.e:
                e2 = Struct(e2_test, e2_test_type, [(i1, j1)])

    # E3 = min{W(i+1,i') + W(i'+1,j-1)}, i+1<i'<j-2
    e3 = STRUCT_NULL
    if not isolated_outer or not i or j == len(seq) - 1:
        for k in range(i + 1, j - 1):
            e3_test = _multi_branch(seq, i, k, j, temp, v_cache, w_cache, emap, True)

            if e3_test and e3_test.e < e3.e:
                e3 = e3_test

    e = _min_struct(e1, e2, e3)
    v_cache[i][j] = e
    return e


def _pair(s: str, i: int, i1: int, j: int, j1: int) -> str:
    """Return a stack representation, a key for the NN maps

    Args:
        s: Sequence being folded
        i: leftmost index
        i1: index to right of i
        j: rightmost index
        j1: index to left of j
    
    Returns:
        str: string representation of the pair
    """

    return (
        (s[i] if i >= 0 else ".")
        + (s[i1] if i1 >= 0 else ".")
        + "/"
        + (s[j] if j >= 0 else ".")
        + (s[j1] if j1 >= 0 else ".")
    )


def _min_struct(*structs: Struct) -> Struct:
    """Return the struct with the lowest free energy that isn't -inf (undef)

    Args:
        structs: Structures being compared
    
    Returns:
        struct: The min free energy structure
    """

    s: Struct = STRUCT_NULL
    for struct in structs:
        if struct.e != -math.inf and struct.e < s.e:
            s = struct
    return s


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
    for pre-calculated free-energies. See SantaLucia and Hicks (2004).

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


def _stack(
    seq: str, i: int, i1: int, j: int, j1: int, temp: float, emap: Energies
) -> float:
    """Get the free energy for a stack.

    Using the indexes i and j, check whether it's at the end of
    the sequence or internal. Then check whether it's a match
    or mismatch, and return.

    Two edge-cases are terminal mismatches and dangling ends.
    The energy of a dangling end is added to the energy of a pair
    where i XOR j is at the sequence's end.

    Args:
        seq: The full folding sequence
        i: The start index on left side of the pair/stack
        i1: The index to the right of i
        j: The end index on right side of the pair/stack
        j1: The index to the left of j
        temp: Temperature in Kelvin

    Returns:
        float: The free energy of the NN pairing
    """

    if any(x >= len(seq) for x in [i, i1, j, j1]):
        return 0.0

    pair = _pair(seq, i, i1, j, j1)
    if any(x == -1 for x in [i, i1, j, j1]):
        # it's a dangling end
        d_h, d_s = emap.DE[pair]
        return _d_g(d_h, d_s, temp)

    if i > 0 and j < len(seq) - 1:
        # it's internal
        d_h, d_s = emap.NN[pair] if pair in emap.NN else emap.INTERNAL_MM[pair]
        return _d_g(d_h, d_s, temp)

    if i == 0 and j == len(seq) - 1:
        # it's terminal
        d_h, d_s = emap.NN[pair] if pair in emap.NN else emap.TERMINAL_MM[pair]
        return _d_g(d_h, d_s, temp)

    if i > 0 and j == len(seq) - 1:
        # it's dangling on left
        d_h, d_s = emap.NN[pair] if pair in emap.NN else emap.TERMINAL_MM[pair]
        d_g = _d_g(d_h, d_s, temp)

        pair_de = seq[i - 1] + seq[i] + "/." + seq[j]
        if pair_de in emap.DE:
            d_h, d_s = emap.DE[pair_de]
            d_g += _d_g(d_h, d_s, temp)
        return d_g

    if i == 0 and j < len(seq) - 1:
        # it's dangling on right
        d_h, d_s = emap.NN[pair] if pair in emap.NN else emap.TERMINAL_MM[pair]
        d_g = _d_g(d_h, d_s, temp)

        pair_de = "." + seq[i] + "/" + seq[j + 1] + seq[j]
        if pair_de in emap.DE:
            d_h, d_s = emap.DE[pair_de]
            d_g += _d_g(d_h, d_s, temp)
        return d_g

    return 0


def _hairpin(seq: str, i: int, j: int, temp: float, emap: Energies) -> float:
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
    pair = _pair(seq, i, i + 1, j, j - 1)

    if emap.COMPLEMENT[hairpin[0]] != hairpin[-1]:
        # not known terminal pair, nothing to close "hairpin"
        raise RuntimeError

    d_g = 0.0
    if emap.TRI_TETRA_LOOPS and hairpin in emap.TRI_TETRA_LOOPS:
        # it's a pre-known hairpin with known value
        d_h, d_s = emap.TRI_TETRA_LOOPS[hairpin]
        d_g = _d_g(d_h, d_s, temp)

    # add penalty based on size
    if hairpin_len in emap.HAIRPIN_LOOPS:
        d_h, d_s = emap.HAIRPIN_LOOPS[hairpin_len]
        d_g += _d_g(d_h, d_s, temp)
    else:
        # it's too large, extrapolate
        d_h, d_s = emap.HAIRPIN_LOOPS[30]
        d_g_inc = _d_g(d_h, d_s, temp)
        d_g += _j_s(hairpin_len, 30, d_g_inc, temp)

    # add penalty for a terminal mismatch
    if hairpin_len > 3 and pair in emap.TERMINAL_MM:
        if pair in emap.TERMINAL_MM:
            d_h, d_s = emap.TERMINAL_MM[pair]
            d_g += _d_g(d_h, d_s, temp)

    # add penalty if length 3 and AT closing, formula 8 from SantaLucia, 2004
    if hairpin_len == 3 and (hairpin[0] == "A" or hairpin[-1] == "A"):
        d_g += 0.5  # convert to entropy

    return d_g


def _bulge(
    seq: str, i: int, i1: int, j: int, j1: int, temp: float, emap: Energies
) -> float:
    """Calculate the free energy associated with a bulge.

        seq: The full folding DNA sequence
        i: The start index of the bulge
        i1: The index to the right of i
        j: The end index of the bulge
        j1: The index to the left of j
        loop: The sequence of the bulge
        temp: Temperature in Kelvin
        emap: Map to DNA/RNA energies

    Returns:
        float: The increment in free energy from the bulge
    """

    loop_len = max(i1 - i - 1, j - j1 - 1)
    if loop_len <= 0:
        raise RuntimeError

    # add penalty based on size
    if loop_len in emap.BULGE_LOOPS:
        d_h, d_s = emap.BULGE_LOOPS[loop_len]
        d_g = _d_g(d_h, d_s, temp)
    else:
        # it's too large for pre-calculated list, extrapolate
        d_h, d_s = emap.BULGE_LOOPS[30]
        d_g = _d_g(d_h, d_s, temp)
        d_g = _j_s(loop_len, 30, d_g, temp)

    if loop_len == 1:
        # if len 1, include the delta G of intervening NN (SantaLucia 2004)
        pair = _pair(seq, i, i1, j, j1)
        assert pair in emap.NN
        d_g += _stack(seq, i, i1, j, j1, temp, emap)

    # penalize AT terminal bonds
    if any(seq[k] == "A" for k in [i, i1, j, j1]):
        d_g += 0.5

    return d_g


def _internal_loop(
    seq: str, i: int, i1: int, j: int, j1: int, temp: float, emap: Energies
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
        i1: The index to the right of i
        j: The index of the end of structure on right side
        j1: The index to the left of j
        temp: Temperature in Kelvin
        emap: Dictionary mapping to energies for DNA/RNA

    Returns:
        float: The free energy associated with the internal loop
    """

    loop_left = i1 - i - 1
    loop_right = j - j1 - 1
    loop_len = loop_left + loop_right

    if loop_left < 1 or loop_right < 1:
        raise RuntimeError

    # single bp mismatch, sum up the two single mismatch pairs
    if loop_left == 1 and loop_right == 1:
        mm_left = _stack(seq, i, i1, j, j1, temp, emap)
        mm_right = _stack(seq, i1 - 1, i1, j1 + 1, j1, temp, emap)
        return mm_left + mm_right

    # apply a penalty based on loop size
    if loop_len in emap.INTERNAL_LOOPS:
        d_h, d_s = emap.INTERNAL_LOOPS[loop_len]
        d_g = _d_g(d_h, d_s, temp)
    else:
        # it's too large an internal loop, extrapolate
        d_h, d_s = emap.INTERNAL_LOOPS[30]
        d_g = _d_g(d_h, d_s, temp)
        d_g = _j_s(loop_len, 30, d_g, temp)

    # apply an asymmetry penalty
    loop_asymmetry = abs(loop_left - loop_right)
    d_g += 0.3 * loop_asymmetry

    # apply penalty based on the mismatching pairs on either side of the loop
    pair_left_mm = _pair(seq, i, i + 1, j, j - 1)
    d_h, d_s = emap.TERMINAL_MM[pair_left_mm]
    d_g += _d_g(d_h, d_s, temp)

    pair_right_mm = _pair(seq, i1 - 1, i1, j1 + 1, j1)
    d_h, d_s = emap.TERMINAL_MM[pair_right_mm]
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
    emap: Energies,
    helix: bool = False,
) -> Struct:
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
        helix: Whether this multibranch is enclosed by a helix
        emap: Map to DNA/RNA energies
    
    Keyword Args:
        helix: Whether V(i, j) bond with one another in a helix

    Returns:
        Struct: A multi-branch structure
    """

    if helix:
        left = _w(seq, i + 1, k, temp, v_cache, w_cache, emap)
        right = _w(seq, k + 1, j - 1, temp, v_cache, w_cache, emap)
    else:
        left = _w(seq, i, k, temp, v_cache, w_cache, emap)
        right = _w(seq, k + 1, j, temp, v_cache, w_cache, emap)

    if not left or not right:
        return STRUCT_NULL

    # gather all branches of this multi-branch structure
    branches: List[Tuple[int, int]] = []

    def add_branch(s: Struct):
        if not s or not s.ij:
            return
        if len(s.ij) == 1:
            branches.append(s.ij[0])
            return
        for i1, j1 in s.ij:
            add_branch(_w(seq, i1, j1, temp, v_cache, w_cache, emap))

    add_branch(left)
    add_branch(right)

    # this isn't multi-branched
    if len(branches) < 2:
        return STRUCT_NULL

    # if there's a helix, i,j counts as well
    if helix:
        branches.append((i, j))

    # count up unpaired bp and asymmetry
    branches_count = len(branches)
    unpaired = 0
    e_sum = 0.0
    for index, (i2, j2) in enumerate(branches):
        _, j1 = branches[(index - 1) % len(branches)]
        i3, j3 = branches[(index + 1) % len(branches)]

        # add energy from unpaired bp to the right
        # of the helix as though it was a dangling end
        # if there's only one bp, it goes to whichever
        # helix (this or the next) has the more favorable energy
        unpaired_left = 0
        unpaired_right = 0
        de = 0.0
        if index == len(branches) - 1 and not helix:
            pass
        elif (i3, j3) == (i, j):
            unpaired_left = i2 - j1 - 1
            unpaired_right = j3 - j2 - 1

            if unpaired_left and unpaired_right:
                de = _stack(seq, i2 - 1, i2, j2 + 1, j2, temp, emap)
            elif unpaired_right:
                de = _stack(seq, -1, i2, j2 + 1, j2, temp, emap)
                if unpaired_right == 1:
                    de = min(_stack(seq, i3, -1, j3, j3 - 1, temp, emap), de)
        elif (i2, j2) == (i, j):
            unpaired_left = j2 - j1 - 1
            unpaired_right = i3 - i2 - 1

            if unpaired_left and unpaired_right:
                de = _stack(seq, i2 - 1, i2, j2 + 1, j2, temp, emap)
            elif unpaired_right:
                de = _stack(seq, i2, i2 + 1, j2, -1, temp, emap)
                if unpaired_right == 1:
                    de = min(_stack(seq, i3 - 1, i3, -1, j3, temp, emap), de)
        else:
            unpaired_left = i2 - j1 - 1
            unpaired_right = i3 - j2 - 1

            if unpaired_left and unpaired_right:
                de = _stack(seq, i2 - 1, i2, j2 + 1, j2, temp, emap)
            elif unpaired_right:
                de = _stack(seq, -1, i2, j2 + 1, j2, temp, emap)
                if unpaired_right == 1:
                    de = min(_stack(seq, i2 - 1, i2, j2 + 1, j2, temp, emap), de)

        e_sum += de
        unpaired += unpaired_right
        assert unpaired_right >= 0

        if (i2, j2) != (i, j):  # add energy
            e_sum += _w(seq, i2, j2, temp, v_cache, w_cache, emap).e

    assert unpaired >= 0

    # penalty for unmatched bp and multi-branch
    a, b, c, d = emap.MULTIBRANCH
    e_multibranch = a + b * len(branches) + c * unpaired

    if unpaired == 0:
        e_multibranch = a + d

    # energy of min-energy neighbors
    e = e_multibranch + e_sum

    # pointer to next structures
    if helix:
        branches.pop()

    return Struct(e, f"BIFURCATION:{str(unpaired)}n/{str(branches_count)}h", branches)


def traceback(i: int, j: int, v_cache: Cache, w_cache: Cache) -> List[Struct]:
    """Traceback thru the V(i,j) and W(i,j) caches to find the structure

    For each step, get to the lowest energy W(i,j) within that block
    Store the structure in W(i,j)
    Inc i and j
    If the next structure is viable according to V(i,j), store as well
    Repeat

    Args:
        i: The leftmost index to start searching in
        j: The rightmost index to start searching in
        v_cache: Energies where i and j bond
        w_cache: Energies/sub-structures between or with i and j

    Returns:
        A list of Structs in the final secondary structure
    """

    # move i,j down-left to start coordinates
    s = w_cache[i][j]
    if "HAIRPIN" not in s.desc:
        while w_cache[i + 1][j] == s:
            i += 1
        while w_cache[i][j - 1] == s:
            j -= 1

    structs: List[Struct] = []
    while True:
        s = v_cache[i][j]
        structs.append(s.with_ij([(i, j)]))

        # it's a hairpin, end of structure
        if not s.ij:
            # set the energy of everything relative to the hairpin
            return _trackback_energy(structs)

        # it's a stack, bulge, etc
        # there's another single structure beyond this
        if len(s.ij) == 1:
            i, j = s.ij[0]
            continue

        # it's a multibranch
        e_sum = 0.0
        structs = _trackback_energy(structs)
        branches: List[Struct] = []
        for i1, j1 in s.ij:
            tb = traceback(i1, j1, v_cache, w_cache)
            if tb and tb[0].ij:
                i2, j2 = tb[0].ij[0]
                e_sum += w_cache[i2][j2].e
                branches += tb

        last = structs[-1]
        structs[-1] = Struct(round(last.e - e_sum, 1), last.desc, list(last.ij))
        return structs + branches

    return _trackback_energy(structs)


def _trackback_energy(structs: List[Struct]) -> List[Struct]:
    """Add energy to each structure, based on how it's W(i,j) differs from the one after

    Args:
        structs: The structures for whom energy is being calculated

    Returns:
        List[Struct]: Structures in the folded DNA with energy
    """

    structs_e: List[Struct] = []
    for index, struct in enumerate(structs):
        e_next = 0.0 if index == len(structs) - 1 else structs[index + 1].e
        structs_e.append(
            Struct(round(struct.e - e_next, 1), struct.desc, list(struct.ij))
        )
    return structs_e
