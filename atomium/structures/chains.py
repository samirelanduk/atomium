"""Contains chains and related polymer classes"""

from .molecules import Molecule

class Chain(Molecule):
    """Base class: :py:class:`Molcule`

    A Chain is a polymer of :py:class:`.Residue` objects that form a
    molecular unit."""

    def __init__(self, *atoms, chain_id=None, **kwargs):
        if chain_id: kwargs["molecule_id"] = chain_id
        Molecule.__init__(self, *atoms, **kwargs)
        for atom in atoms:
            atom._chain = self
