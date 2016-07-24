from .molecules import AtomicStructure

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
        self._small_molecules.add(small_molecule)
