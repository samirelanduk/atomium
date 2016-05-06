from . import atomic

class Monomer(atomic.Molecule):

    def __init__(self, monomer_id, monomer_name, *atoms):
        atomic.Molecule.__init__(self, *atoms, molecule_id=monomer_id, molecule_name=monomer_name)
        self.monomer_id = self.molecule_id
        self.monomer_name = self.molecule_name
        del self.__dict__["molecule_id"]
        del self.__dict__["molecule_name"]


    def __repr__(self):
        return "<Monomer (%s)>" % self.monomer_name
