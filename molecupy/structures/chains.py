from .molecules import AtomicStructure

class ResiduicStructure(AtomicStructure):

    def __init__(self, *residues):
        self._residues = set(residues)


    def __repr__(self):
        return "<ResiduicStructure (%i residues)>" % len(self._residues)
