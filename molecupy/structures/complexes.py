from . import ResiduicStructure, Chain
from ..exceptions import DuplicateChainsError

class Complex(ResiduicStructure):

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
        return self._complex_id


    def complex_name(self, complex_name=None):
        if complex_name:
            if not isinstance(complex_name, str):
                raise TypeError(
                 "complex_name must be str, not '%s'" % str(complex_name)
                )
            self._complex_name = complex_name
        else:
            return self._complex_name


    def chains(self):
        return set(self._chains)


    def model(self):
        return self._model


    def add_chain(self, chain):
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
        self._chains.remove(chain)
