from collections import Counter
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


    def atoms(self, atom_type="all"):
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


    def formula(self, atom_type="all", include_hydrogens=False):
        return Counter([
         atom.element() for atom in self.atoms(atom_type=atom_type)
          if include_hydrogens or atom.element().upper() != "H"
        ])


    def get_atom_by_id(self, atom_id):
        if not isinstance(atom_id, int):
            raise TypeError("Atom ID search must be by int, not '%s'" % str(atom_id))
        for atom in self.atoms():
            if atom.atom_id() == atom_id:
                return atom


    def get_atoms_by_element(self, element, atom_type="all"):
        if not isinstance(element, str):
            raise TypeError("Atom element search must be by str, not '%s'" % str(element))
        return set([
         atom for atom in self.atoms(atom_type=atom_type) if atom.element() == element
        ])


    def get_atom_by_element(self, element, atom_type="all"):
        if not isinstance(element, str):
            raise TypeError("Atom element search must be by str, not '%s'" % str(element))
        for atom in self.atoms(atom_type=atom_type):
            if atom.element() == element:
                return atom


    def get_atoms_by_name(self, atom_name, atom_type="all"):
        if not isinstance(atom_name, str):
            raise TypeError("Atom name search must be by str, not '%s'" % str(atom_name))
        return set([
         atom for atom in self.atoms(atom_type=atom_type) if atom.atom_name() == atom_name
        ])


    def get_atom_by_name(self, atom_name, atom_type="all"):
        if not isinstance(atom_name, str):
            raise TypeError("Atom name search must be by str, not '%s'" % str(atom_name))
        for atom in self.atoms(atom_type=atom_type):
            if atom.atom_name() == atom_name:
                return atom
