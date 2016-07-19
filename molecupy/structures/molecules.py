from .atoms import Atom
from ..exceptions import NoAtomsError

class AtomicStructure:

    def __init__(self, *atoms):
        if len(atoms) == 0:
            raise NoAtomsError("Cannot make an AtomicStructure with no atoms")
        for atom in atoms:
            if not isinstance(atom, Atom):
                raise TypeError(
                 "Can only make AtomicStructures with Atoms, not '%s'" % str(atom)
                )
        self._atoms = set(atoms)


    def __repr__(self):
        return "<AtomicStructure (%i atoms)>" % len(self._atoms)
