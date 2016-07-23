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


    def __getattr__(self, attribute):
        if attribute == "_atoms":
            atoms = set()
            for residue in self.residues():
                atoms.update(residue.atoms())
            return atoms
        else:
            return self.__getattribute__(attribute)


    def residues(self, include_missing=True):
        if include_missing:
            return set(self._residues)
        else:
            return set([res for res in self._residues if not res.is_missing()])


    def add_residue(self, residue):
        if not isinstance(residue, Residue):
            raise TypeError(
             "Can only add Residues to ResiduicStructures, not '%s'" % str(residue)
            )
        self._residues.add(residue)


    def remove_residue(self, residue):
        self._residues.remove(residue)


    def get_residue_by_id(self, residue_id):
        if not isinstance(residue_id, str):
            raise TypeError(
             "Residue ID search must be by str, not '%s'" % str(residue_id)
            )
        for residue in self.residues():
            if residue.residue_id() == residue_id:
                return residue
