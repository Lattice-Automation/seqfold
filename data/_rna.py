"""Parse the energies from Turner, 2004 to rna.py."""

from os import path

DIR = path.dirname(path.realpath(__file__))
DIR_DATA = path.join(DIR, "..", "data")

DE = path.join(DIR_DATA, "rna.dangle.txt")  # input files
DE_DH = path.join(DIR_DATA, "rna.dangle.dh.txt")

LOOP = path.join(DIR_DATA, "rna.loop.txt")
LOOP_DH = path.join(DIR_DATA, "rna.loop.dh.txt")

STACK = path.join(DIR_DATA, "rna.stack.txt")
STACK_DH = path.join(DIR_DATA, "rna.stack.dh.txt")

TSTACK = path.join(DIR_DATA, "rna.tstack.txt")
TSTACK_DH = path.join(DIR_DATA, "rna.tstack.dh.txt")

RNA_PY = path.join(DIR_DATA, "..", "seqfold", "rna.py")  # output

RNA_COMPLEMENT = {"A": "U", "U": "A", "G": "C", "C": "G", "N": "N"}

RNA_EXPORT = """
# the energies are the same for each loop stack in the
# reverse complementary direction
RNA_NN.update({k[::-1]: v for k, v in RNA_NN.items()})
RNA_INTERNAL_MM.update(
    {k[::-1]: v for k, v in RNA_INTERNAL_MM.items() if k[::-1] not in RNA_INTERNAL_MM}
)
RNA_TERMINAL_MM.update(
    {k[::-1]: v for k, v in RNA_TERMINAL_MM.items() if k[::-1] not in RNA_TERMINAL_MM}
)
RNA_DE.update({k[::-1]: v for k, v in RNA_DE.items() if k[::-1] not in RNA_DE})

RNA_ENERGIES = Energies(
    RNA_BULGE_LOOPS,
    RNA_COMPLEMENT,
    RNA_DE,
    RNA_HAIRPIN_LOOPS,
    RNA_MULTIBRANCH,
    RNA_INTERNAL_LOOPS,
    RNA_INTERNAL_MM,
    RNA_NN,
    RNA_TERMINAL_MM,
    None,
)
"""


def parse():
    """
    These energies are originally in a non-standard txt format.
    Here am parsing to a Python format usable by seqfold.
    """

    template = path.join(DIR_DATA, "rna.py.template")
    outfile = ""
    with open(template, "r") as temp_file:
        outfile = temp_file.read()

    stack_nn_map, stack_mm_map, tstack_map = parse_stack()

    de_map = parse_de()

    iloops, bloops, hloops = parse_loops()

    outfile += "RNA_NN: BpEnergy = " + str(stack_nn_map) + "\n\n"
    outfile += "RNA_INTERNAL_MM: BpEnergy = " + str(stack_mm_map) + "\n\n"
    outfile += "RNA_TERMINAL_MM: BpEnergy = " + str(tstack_map) + "\n\n"
    outfile += "RNA_DE: BpEnergy = " + str(de_map) + "\n\n"
    outfile += "RNA_INTERNAL_LOOPS: LoopEnergy = " + str(iloops) + "\n\n"
    outfile += "RNA_BULGE_LOOPS: LoopEnergy = " + str(bloops) + "\n\n"
    outfile += "RNA_HAIRPIN_LOOPS: LoopEnergy = " + str(hloops) + "\n\n"
    outfile += RNA_EXPORT

    with open(RNA_PY, "w") as out:
        out.write(outfile)


def parse_de():
    """Parse dangling ends file: DE + DE_DH."""

    # map from dangling end stack to tuple with dg
    bps = ["A", "C", "G", "U"]

    def parse_de_file(filename):
        e_map = {}

        with open(filename, "r") as f:
            lines = f.readlines()
            lines = [l.strip() for l in lines if l.strip()]

            for bi in range(len(lines) // 9):
                ci = bi * 9

                top = lines[ci + 5]
                vs = [
                    float(v) if v != "." else 0.0 for v in lines[ci + 8].strip().split()
                ]

                assert len(vs) == 16, len(vs)

                topleft = top.split()[0][0]

                for i in range(16):
                    topright = "."
                    if bi < 4:
                        topright = bps[i % 4]
                    bottomleft = bps[i // 4]
                    bottomright = "."
                    if bi >= 4:
                        bottomright = bps[i % 4]

                    stack = topleft + topright + "/" + bottomleft + bottomright

                    e_map[stack] = vs[i]

        assert "UC/A." in e_map

        return e_map

    dg_map = parse_de_file(DE)
    dh_map = parse_de_file(DE_DH)

    dh_ds_map = {}
    for k, dg in dg_map.items():
        dh = dh_map[k]
        dh_ds_map[k] = (dh, _ds(dg, dh))

    return dh_ds_map


def parse_loops():
    """Parse the loop initiation energies. Return three map: internal, bulge, hairpin"""

    def parse_loop_file(filename):
        loop_map = {}
        with open(filename, "r") as f:
            lines = f.readlines()
            lines = lines[4:]

            for i, line in enumerate(lines):
                vs = line.split()[1:]
                vs_float = [float(v) if v != "." else 0.0 for v in vs]
                loop_map[i + 1] = vs_float
        return loop_map

    dg_loop_map = parse_loop_file(LOOP)
    dh_loop_map = parse_loop_file(LOOP_DH)

    iloops = {}
    bloops = {}
    hloops = {}
    for k, dg in dg_loop_map.items():
        dh = dh_loop_map[k]

        iloops[k] = (dh[0], _ds(dg[0], dh[0]))
        bloops[k] = (dh[1], _ds(dg[1], dh[1]))
        hloops[k] = (dh[2], _ds(dg[2], dh[2]))
    return (iloops, bloops, hloops)


def parse_stack():
    """Parse a stack file with matching or mismatching bp

    Return two maps. one for matching and one for mismatching stacks
    """

    bps = "ACGU"

    def parse_stack_file(filename):
        stack_map = {}

        with open(filename, "r") as f:
            lines = [l.strip() for l in f.readlines() if l.strip()]
            lines = lines[15:]

            for k in range(4):
                bi = k * 12

                vals = []
                for i in range(4):
                    vals = [
                        float(v) if v != "." else 0.0 for v in lines[bi + 8 + i].split()
                    ]

                    assert len(vals) == 16

                    for j in range(16):
                        topleft = bps[k]
                        topright = bps[i]
                        bottomleft = bps[j // 4]
                        bottomright = bps[j % 4]

                        stack = topleft + topright + "/" + bottomleft + bottomright
                        stack_map[stack] = vals[j]

        return stack_map

    stack_dg_map = parse_stack_file(STACK)
    stack_dh_map = parse_stack_file(STACK_DH)

    stack_map = {}
    for k, dg in stack_dg_map.items():
        dh = stack_dh_map[k]
        stack_map[k] = (dh, _ds(dg, dh))

    # separate the stack map into matches and mismatches
    stack_nn_map = {}
    stack_mm_map = {}
    for k, v in stack_map.items():
        if RNA_COMPLEMENT[k[0]] == k[3] and RNA_COMPLEMENT[k[1]] == k[-1]:
            stack_nn_map[k] = v
        else:
            stack_mm_map[k] = v

    tstack_dg_map = parse_stack_file(TSTACK)
    tstack_dh_map = parse_stack_file(TSTACK_DH)

    tstack_map = {}
    for k, dg in tstack_dg_map.items():
        dh = tstack_dh_map[k]
        tstack_map[k] = (dh, _ds(dg, dh))

    return (stack_nn_map, stack_mm_map, tstack_map)


def _ds(dg, dh):
    ds = round(((dg - dh) / -310.15) * 1000.0, 1)
    if ds == -0.0:
        ds = 0.0
    return ds


if __name__ == "__main__":
    parse()
