from .molecules import AtomicStructure, SmallMolecule
from .chains import Chain

class Model(AtomicStructure):

    def __init__(self):
        self._small_molecules = set()
        self._chains = set()


    def __getattr__(self, attribute):
        if attribute == "_atoms":
            return set()
        else:
            return self.__getattribute__(attribute)


    def small_molecules(self):
        return set(self._small_molecules)


    def add_small_molecule(self, small_molecule):
        if not isinstance(small_molecule, SmallMolecule):
            raise TypeError(
             "Can only add SmallMolecule to Model, not '%s'" % str(small_molecule)
            )
        self._small_molecules.add(small_molecule)
        small_molecule._model = self


    def remove_small_molecule(self, small_molecule):
        self._small_molecules.remove(small_molecule)
        small_molecule._model = None


    def get_small_molecule_by_id(self, molecule_id):
        if not isinstance(molecule_id, str):
            raise TypeError(
             "Small molecule ID search must be by str, not '%s'" % str(molecule_id)
            )
        for molecule in self.small_molecules():
            if molecule.molecule_id() == molecule_id:
                return molecule


    def get_small_molecule_by_name(self, molecule_name):
        if not isinstance(molecule_name, str):
            raise TypeError(
             "Small molecule name search must be by str, not '%s'" % str(molecule_name)
            )
        for molecule in self.small_molecules():
            if molecule.molecule_name() == molecule_name:
                return molecule


    def get_small_molecules_by_name(self, molecule_name):
        if not isinstance(molecule_name, str):
            raise TypeError(
             "Small molecule name search must be by str, not '%s'" % str(molecule_name)
            )
        return set([molecule for molecule in self.small_molecules()
         if molecule.molecule_name() == molecule_name])


    def chains(self):
        return set(self._chains)


    def add_chain(self, chain):
        if not isinstance(chain, Chain):
            raise TypeError(
             "Can only add Chain to Model, not '%s'" % str(chain)
            )
        self._chains.add(chain)
        chain._model = self


    def remove_chain(self, chain):
        self._chains.remove(chain)
        chain._model = None


    def get_chain_by_id(self, chain_id):
        if not isinstance(chain_id, str):
            raise TypeError(
             "Chain ID search must be by str, not '%s'" % str(chain_id)
            )
        for chain in self.chains():
            if chain.chain_id() == chain_id:
                return chain
