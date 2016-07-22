from .molecules import AtomicStructure

class ResiduicStructure(AtomicStructure):

    def __init__(self, *residues):
        self._residues = set(residues)
