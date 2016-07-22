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


    def contacts_with(self, other_atomic_structure, distance=4, include_hydrogens=True):
        contacts = set()
        our_atoms = list(self.atoms(atom_type="pdb"))
        their_atoms = [
         a for a in list(other_atomic_structure.atoms(atom_type="pdb"))
          if a not in our_atoms
        ]
        if not include_hydrogens:
            our_atoms = [
             atom for atom in our_atoms if atom.element().upper() != "H"
            ]
            their_atoms = [
             atom for atom in their_atoms if atom.element().upper() != "H"
            ]
        for atom in our_atoms:
            for other_atom in their_atoms:
                if atom.distance_to(other_atom) <= distance:
                    contacts.add(frozenset((atom, other_atom)))
        return contacts


    def internal_contacts(self, distance=4, include_hydrogens=True):
        contacts = set()
        atoms = list(self.atoms(atom_type="pdb"))
        if not include_hydrogens:
            atoms = [atom for atom in atoms if atom.element().upper() != "H"]
        for index, atom in enumerate(atoms[:-1]):
            too_close_atoms = set((atom,))
            for bonded_atom in atom.bonded_atoms():
                too_close_atoms.add(bonded_atom)
            second_tier_atoms = set()
            for close_atom in too_close_atoms:
                for bonded_atom in close_atom.bonded_atoms():
                    second_tier_atoms.add(bonded_atom)
            too_close_atoms.update(second_tier_atoms)
            for other_atom in atoms[index + 1:]:
                if other_atom not in too_close_atoms and atom.distance_to(other_atom) <= distance:
                    contacts.add(frozenset((atom, other_atom)))
        return contacts


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



class SmallMolecule(AtomicStructure):

    def __init__(self, molecule_id, molecule_name, *atoms):
        if not isinstance(molecule_id, str):
            raise TypeError("'%s' is not a valid molecule_id" % str(molecule_id))
        self._molecule_id = molecule_id
        if not isinstance(molecule_name, str):
            raise TypeError("'%s' is not a valid molecule_name" % str(molecule_name))
        self._molecule_name = molecule_name
        AtomicStructure.__init__(self, *atoms)


    def __repr__(self):
        return "<SmallMolecule %s (%s)>" % (self._molecule_id, self._molecule_name)
