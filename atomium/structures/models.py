"""Contains the Model class."""

from .molecules import AtomicStructure
from .chains import ResidueStructure

class ChainStructure:

    def chains(self, chain_id=None, name=None):
        """Returns the py:class:`.Chain` objects in the structure. It can be
        given search criteria if you wish.

        :param str chain_id: Filter by chain ID.
        :param str name: Filter by name.
        :rtype: ``set``"""

        chains = set()
        for atom in self.atoms():
            chains.add(atom.chain())
        if chain_id:
            chains = set(filter(lambda c: c.chain_id() == chain_id, chains))
        if name:
            chains = set(filter(lambda c: c.name() == name, chains))
        return chains


    def chain(self, *args, **kwargs):
        """Returns the first py:class:`.Chain` object in the structure which
        matches the given criteria.

        :param str chain_id: Filter by chain ID.
        :param str name: Filter by name.
        :rtype: ``Chain``"""

        chains = self.chains(*args, **kwargs)
        for chain in chains: return chain



class Model(AtomicStructure, ResidueStructure, ChainStructure):
    """Base class: :py:class:`.AtomicStructure`

    Represents molecular systems.

    :param \*atoms: The atoms that make up the model."""

    def __init__(self, *atoms):
        AtomicStructure.__init__(self, *atoms)
        for atom in self._atoms:
            atom._model = self


    def add_atom(self, atom, *args, **kwargs):
        AtomicStructure.add_atom(self, atom, *args, **kwargs)
        atom._model = self


    def remove_atom(self, atom, *args, **kwargs):
        AtomicStructure.remove_atom(self, atom, *args, **kwargs)
        atom._model = None
