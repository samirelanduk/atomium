from .molecules import AtomicStructure, Residue
from ..exceptions import NoResiduesError

class ResiduicStructure(AtomicStructure):

    def __init__(self, *residues):
        if len(residues) == 0:
            raise NoResiduesError("Cannot make a ResiduicStructure with no residues")
        for residue in residues:
            if not isinstance(residue, Residue):
                raise TypeError(
                 "Can only make ResiduicStructure with Residuess, not '%s'" % str(residue)
                )
        self._residues = set(residues)


    def __repr__(self):
        return "<ResiduicStructure (%i residues)>" % len(self._residues)


    def residues(self, include_missing=True):
        if include_missing:
            return set(self._residues)
        else:
            return set([res for res in self._residues if not res.is_missing()])
