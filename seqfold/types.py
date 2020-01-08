"""Types shared between dna.py and rna.py
"""

from typing import Dict, Optional, Tuple, Union


Comp = Dict[str, str]
MultiBranch = Tuple[float, float, float]
BpEnergy = Dict[str, Tuple[float, float]]
LoopEnergy = Dict[int, Tuple[float, float]]


class Energies:
    def __init__(
        self,
        bulge_loops: LoopEnergy,
        complement: Comp,
        de: BpEnergy,
        hairpin_loops: LoopEnergy,
        multibranch: MultiBranch,
        internal_loops: LoopEnergy,
        internal_mm: BpEnergy,
        nn: BpEnergy,
        terminal_mm: BpEnergy,
        tri_tetra_loops: Optional[BpEnergy] = None,
    ):
        self.BULGE_LOOPS = bulge_loops
        self.COMPLEMENT = complement
        self.DE = de
        self.HAIRPIN_LOOPS = hairpin_loops
        self.MULTIBRANCH = multibranch
        self.INTERNAL_LOOPS = internal_loops
        self.INTERNAL_MM = internal_mm
        self.NN = nn
        self.TERMINAL_MM = terminal_mm
        self.TRI_TETRA_LOOPS = tri_tetra_loops
