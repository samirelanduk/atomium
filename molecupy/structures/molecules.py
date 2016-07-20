from .atoms import Atom, PdbAtom
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


    def atoms(self, atom_type="pdb"):
        if not isinstance(atom_type, str):
            raise TypeError("atom_type must be str, not '%s'" % str(atom_type))
        if atom_type == "pdb":
            return set(
             [atom for atom in self._atoms if isinstance(atom, PdbAtom)]
            )
        elif atom_type == "generic":
            return set(
             [atom for atom in self._atoms if not isinstance(atom, PdbAtom)]
            )
        elif atom_type == "all":
            return set(self._atoms)
        else:
            raise ValueError("'%s' is not a valid atom_type" % atom_type)


    def add_atom(self, atom):
        if not isinstance(atom, Atom):
            raise TypeError(
             "Can only make add Atoms to AtomicStructures, not '%s'" % str(atom)
            )
        self._atoms.add(atom)


    def remove_atom(self, atom):
        if atom in self._atoms:
            self._atoms.remove(atom)


    def mass(self, atom_type="all"):
        return sum([atom.mass() for atom in self.atoms(atom_type)])
