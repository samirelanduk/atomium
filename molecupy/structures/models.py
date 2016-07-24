from .molecules import AtomicStructure, SmallMolecule

class Model(AtomicStructure):

    def __init__(self):
        self._small_molecules = set()


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
