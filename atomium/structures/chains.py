"""Contains chains and related polymer classes"""

from .molecules import Molecule
from .exceptions import SequenceConnectivityError

class ResidueStructure:
    """This is an interface class which confers upon an object the ability to
    retrieve the :py:class:`.Residue` objects (unordered) it contains. It should
    not be instantiated directly.

    Only classes which inherit from :py:class:`.AtomicStructure` should inherit
    this class, because it requires the :py:meth:`~AtomicStructure.atoms`
    method."""

    def residues(self, residue_id=None, name=None):
        """Returns the py:class:`.Residue` objects in the structure. It can be
        given search criteria if you wish.

        :param str residue_id: Filter by residue ID.
        :param str name: Filter by name.
        :rtype: ``set``"""

        res = set()
        for atom in self.atoms():
            res.add(atom.residue())
        if residue_id:
            res = set(filter(lambda r: r.residue_id() == residue_id, res))
        if name:
            res = set(filter(lambda r: r.name() == name, res))
        return res


    def residue(self, *args, **kwargs):
        """Returns the first py:class:`.Residue` object in the structure which
        matches the given criteria.

        :param str residue_id: Filter by residue ID.
        :param str name: Filter by name.
        :rtype: ``Residue``"""

        residues = self.residues(*args, **kwargs)
        for res in residues: return res



class ResidueSequence(ResidueStructure):
    """Base class: :py:class:`ResidueStructure`

    This is an interface class which confers upon an object the ability to
    retrieve the :py:class:`.Residue` objects it contains, ordered by residue
    connectivity. It should not be instantiated directly.

    Only classes which inherit from :py:class:`.AtomicStructure` should inherit
    this class, because it requires the :py:meth:`~AtomicStructure.atoms`
    method."""

    @staticmethod
    def verify(sequence):
        """A static method for checking that the residues in a sequence are all
        connected together, and that there are no gaps in the sequence.

        :param ResidueSequence sequence: The sequence to check.
        :raises SequenceConnectivityError: if not properly connected.
        :returns: ``True`` if the test passes."""

        residues = set()
        for atom in sequence.atoms():
            residues.add(atom.residue())
        seq = [list(residues)[0]]
        while seq[-1].next() is not None and seq[-1].next() in residues:
            seq.append(seq[-1].next())
        while seq[-1].previous() is not None and seq[-1].previous() in residues:
            seq.append(seq[-1].previous())
        if set(seq) != residues:
            raise SequenceConnectivityError(
             "{} missing from sequence {} - check connections".format(
              residues, residues - set(seq)
             )
            )
        return True


    def residues(self, *args, **kwargs):
        """Returns the py:class:`.Residue` objects in the structure. It can be
        given search criteria if you wish.

        :param str residue_id: Filter by residue ID.
        :param str name: Filter by name.
        :rtype: ``tuple``"""

        all_res = ResidueStructure.residues(self, *args, **kwargs)
        initial_res = list(all_res)[0]
        full_sequence = [initial_res]
        while full_sequence[-1].next():
            full_sequence.append(full_sequence[-1].next())
        while full_sequence[0].previous():
            full_sequence.insert(0, full_sequence[0].previous())
        return tuple(sorted(all_res, key=lambda k: full_sequence.index(k)))



class Chain(Molecule, ResidueSequence):
    """Base class: :py:class:`Molcule`

    A Chain is a polymer of :py:class:`.Residue` objects that form a
    molecular unit."""

    def __init__(self, *atoms, chain_id=None, **kwargs):
        if chain_id: kwargs["molecule_id"] = chain_id
        Molecule.__init__(self, *atoms, **kwargs)
        ResidueSequence.verify(self)
        for atom in atoms:
            atom._chain = self


    def __repr__(self):
        return "<Chain {}({} residues)>".format(
         self._name + " " if self._name else "",
         len(self.residues())
        )


    def add_atom(self, atom, *args, **kwargs):
        Molecule.add_atom(self, atom, *args, **kwargs)
        atom._chain = self


    def remove_atom(self, atom, *args, **kwargs):
        Molecule.remove_atom(self, atom, *args, **kwargs)
        atom._chain = None
