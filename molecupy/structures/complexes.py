from . import ResiduicStructure, Chain

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
        self._complex_id = complex_id
        self._complex_name = complex_name
        self._chains = set(chains)


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
