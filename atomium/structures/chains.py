"""Contains chains and related polymer classes"""

from .molecules import Molecule

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



class Chain(Molecule):
    """Base class: :py:class:`Molcule`

    A Chain is a polymer of :py:class:`.Residue` objects that form a
    molecular unit."""

    def __init__(self, *atoms, chain_id=None, **kwargs):
        if chain_id: kwargs["molecule_id"] = chain_id
        Molecule.__init__(self, *atoms, **kwargs)
        for atom in atoms:
            atom._chain = self
