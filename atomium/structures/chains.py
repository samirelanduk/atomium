"""This module contains chains and related polymer classes."""

from .molecules import Molecule, Residue, AtomicStructure
from .exceptions import SequenceConnectivityError

class ResidueSequence:
    """Base class: :py:class:`ResidueStructure`

    This is an interface class which confers upon an object the ability to
    retrieve the :py:class:`.Residue` objects it contains, ordered by residue
    connectivity. It should not be instantiated directly.

    They are iterable and indexable.

    Only classes which inherit from :py:class:`.AtomicStructure` should inherit
    this class, because it requires the :py:meth:`~.AtomicStructure.atoms`
    method."""

    @staticmethod
    def verify(sequence):
        """A static method for checking that the residues in a sequence are all
        connected together, and that there are no gaps in the sequence.

        :param ResidueSequence sequence: The sequence to check.
        :raises SequenceConnectivityError: if not properly connected.
        :returns: ``True`` if the test passes."""

        residues = set()
        for atom in sequence._atoms:
            residues.add(atom.residue)
        if residues:
            seq = [list(residues)[0]]
            while seq[-1].next is not None and seq[-1].next in residues:
                seq.append(seq[-1].next)
            while seq[-1].previous is not None and seq[-1].previous in residues:
                seq.append(seq[-1].previous)
            if set(seq) != residues:
                raise SequenceConnectivityError(
                 "{} missing from sequence {} - check connections".format(
                  residues, residues - set(seq)
                 )
                )
        return True


    def __len__(self):
        return len(self.residues())


    def __iter__(self):
        return iter(self.residues())


    def __getitem__(self, index):
        return self.residues()[index]


    @property
    def length(self):
        """Returns the length of the structure in residues.

        :rtype: ``int``"""

        return len(self)


    def residues(self, *args, **kwargs):
        """Returns the :py:class:`.Residue` objects in the structure. It can be
        given search criteria if you wish.

        :param str residue_id: Filter by residue ID.
        :param str name: Filter by name.
        :rtype: ``tuple``"""

        all_res = AtomicStructure.residues(self, *args, **kwargs)
        if len(all_res) == 0: return tuple()
        initial_res = list(all_res)[0]
        full_sequence = [initial_res]
        while full_sequence[-1].next:
            full_sequence.append(full_sequence[-1].next)
        while full_sequence[0].previous:
            full_sequence.insert(0, full_sequence[0].previous)
        return tuple(sorted(all_res, key=lambda k: full_sequence.index(k)))



class Chain(ResidueSequence, Molecule):
    """Base classes: :py:class:`.Molecule` and :py:class:`ResidueSequence`

    A Chain is a polymer of :py:class:`.Residue` objects that form a
    molecular unit. They are iterable and indexable.

    :param \*atoms: The :py:class:`.Atom` objects that make up the structure.\
    These can also be :py:class:`.AtomicStructure` objects, in which case the\
    atoms of that structure will be used in its place.
    :param str id: A unique str ID for the chain. Uniqueness is not\
    actually enforced.
    :param str name: A name for the chain.
    :raises TypeError: if non-atoms or AtomicStructures are given.
    :raises TypeError: if the chain_id is not str."""

    def __init__(self, *atoms, **kwargs):
        Molecule.__init__(self, *atoms, **kwargs)
        ResidueSequence.verify(self)
        for atom in self._atoms:
            atom._chain = self


    def __repr__(self):
        return "<Chain {}({} residues)>".format(
         self._name + " " if self._name else "",
         len(self.residues())
        )



class Site(AtomicStructure):
    """Base classes: :py:class:`.AtomicStructure`

    A Site is a set of residues which bind to a ligand molecule.

    :param \*atoms: The :py:class:`.Atom` objects that make up the structure.\
    These can also be :py:class:`.AtomicStructure` objects, in which case the\
    atoms of that structure will be used in its place.
    :param Molecule ligand: the ligand molecule for this site.
    :raises TypeError: if non-atoms or AtomicStructures are given.
    :raises TypeError: if the ligand is not a :py:class:`.Molecule`."""

    def __init__(self, *atoms, ligand=None):
        if ligand is not None and not isinstance(ligand, Molecule):
            raise TypeError("ligand {} is not a Molecule".format(ligand))
        AtomicStructure.__init__(self, *atoms)
        self._ligand = ligand


    def __repr__(self):
        if self._ligand:
            return "<'{}' Site ({} residues)>".format(
             self._ligand._id, len(self.residues())
            )
        return "<Site ({} residues)>".format(len(self.residues()))


    @property
    def ligand(self):
        """The ligand :py:class:`.Molecule` for this site.

        :raises TypeError: if the ligand given is not a :py:class:`.Molecule`.
        :rtype: ``Molecule``"""

        return self._ligand


    @ligand.setter
    def ligand(self, ligand):
        if not isinstance(ligand, Molecule):
            raise TypeError("ligand {} is not a Molecule".format(ligand))
        self._ligand = ligand


    def residues(self, *args, **kwargs):
        """Returns the :py:class:`.Residue` objects in the structure, including
        water molecules. It can be given search criteria if you wish.

        :param str residue_id: filter by residue ID.
        :param str name: filter by name.
        :rtype: ``set``"""

        residues = AtomicStructure.residues(self, *args, **kwargs)
        for atom in self.atoms():
            molecule = atom.molecule
            if molecule and not isinstance(molecule, Chain):
                residues.add(molecule)
        return residues
