"""Contains classes for structures made of atoms."""

from .atoms import Atom

class AtomicStructure:
    """Represents structures made of :py:class:`.Atom`s, which tends to be
    rather a lot of things in practice. This class would not generally be
    instantiated directly, and is here to be a parent class to other, more
    specific entities.

    AtomicStructures are containers of their atoms.

    :param \*atoms: The :py:class:`.Atom`s that make up the structure.
    :raises TypeError: if non-atoms are given."""

    def __init__(self, *atoms):
        if not all(isinstance(atom, Atom) for atom in atoms):
            non_atoms = [atom for atom in atoms if not isinstance(atom, Atom)]
            raise TypeError(
             "AtomicStructures need atoms, not '{}'".format(non_atoms[0])
            )
        self._atoms = set(atoms)


    def __repr__(self):
        return "<AtomicStructure ({} atoms)>".format(len(self._atoms))


    def __contains__(self, member):
        return member in self._atoms


    def atoms(self, element=None):
        """Returns the :py:class:`.Atom` objects in the structure. You can
        filter these by element if you wish.

        :param str element: If given, only atoms whose element matches this\
        will be returned.
        :rtype: ``set``"""

        atoms = set(self._atoms)
        if element:
            atoms = set(filter(lambda a: a.element() == element, atoms))
        return atoms


    def atom(self, *args, **kwargs):
        """Returns the first :py:class:`.Atom` that matches the criteria given.

        :param str element: If given, only atoms whose element matches this\
        will be searched.
        :rtype: ``Atom``"""

        atoms = self.atoms(*args, **kwargs)
        for atom in atoms: return atom


    def add_atom(self, atom):
        """Adds an :py:class:`.Atom` to the structure.

        :param Atom atom: The atom to add.
        :raises TypeError: if the atom given is not an Atom."""

        if not isinstance(atom, Atom):
            raise TypeError("Can only add atoms, not '{}'".format(atom))
        self._atoms.add(atom)


    def remove_atom(self, atom):
        """Removes an :py:class:`.Atom` from the structure.

        :param Atom atom: The atom to remove."""

        self._atoms.remove(atom)
