"""Contains classes pertaining to complexes and multi-chain assemblies."""

from . import ResiduicStructure, Chain
from ..exceptions import DuplicateChainsError

class Complex(ResiduicStructure):
    """Base class: :py:class:`.ResiduicStructure`

    Represents complexes of multiple :py:class:`.Chain` objects.

    :param str complex_id: The complex's unique ID.
    :param str complex_name: The complex's name.
    :param \*chains: The chains to create the complex from."""

    def __init__(self, complex_id, complex_name, *chains):
        if not isinstance(complex_id, str):
            raise TypeError(
             "complex_id must be str, not '%s'" % str(complex_id)
            )
        if not isinstance(complex_name, str):
            raise TypeError(
             "complex_name must be str, not '%s'" % str(complex_name)
            )
        for chain in chains:
            if not isinstance(chain, Chain):
                raise TypeError(
                 "Can only make Complexes with Chains, not '%s'" % str(chain)
                )
            chain._complex = self
        self._complex_id = complex_id
        self._complex_name = complex_name
        self._chains = set(chains)
        self._model = None


    def __repr__(self):
        return "<Complex '%s' (%i chains)>" % (
         self._complex_name, len(self._chains)
        )


    def __getattr__(self, attribute):
        if attribute == "_residues":
            residues = set()
            for chain in self._chains:
                residues.update(chain.residues(include_missing=True))
            return residues
        elif attribute == "_atoms":
            atoms = set()
            for residue in self._residues:
                atoms.update(residue.atoms(atom_type="all"))
            return atoms
        else:
            return self.__getattribute__(attribute)


    def complex_id(self):
        """Returns the complex's ID.

        :rtype: ``str``"""

        return self._complex_id


    def complex_name(self, complex_name=None):
        """Returns or sets the complex's name.

        :param str complex_name: If given, the complex's name will be set to this.
        :rtype: ``str``"""

        if complex_name:
            if not isinstance(complex_name, str):
                raise TypeError(
                 "complex_name must be str, not '%s'" % str(complex_name)
                )
            self._complex_name = complex_name
        else:
            return self._complex_name


    def chains(self):
        """Returns the :py:class:`.Chain` objects in this complex.

        :returns: ``set`` of :py:class:`.Chain` objects"""

        return set(self._chains)


    def model(self):
        """Returns the :py:class:`.Model` that the complex inhabits.

        :rtype: ``Model``"""

        return self._model


    def add_chain(self, chain):
        """Adds a :py:class:`.Chain` to the structure.

        :param Chain chain: The chain to add."""

        if not isinstance(chain, Chain):
            raise TypeError(
             "Can only add Chains to Complexes, not '%s'" % str(chain)
            )
        if chain.chain_id() in [chain.chain_id() for chain in self.chains()]:
            raise DuplicateChainsError(
             "Cannot add chain with ID %s to %s as there is already a chain with that ID" % (
              chain.chain_id(), self
             )
            )
        self._chains.add(chain)


    def remove_chain(self, chain):
        """Adds a :py:class:`.Chain` to the structure.

        :param Chain chain: The chain to add."""

        self._chains.remove(chain)
