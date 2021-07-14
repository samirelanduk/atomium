import numpy as np

class Model:

    def __init__(self, chains=None, ligands=None, waters=None):
        self.chains, self.ligands, self.waters = chains, ligands, waters
    

    def __repr__(self):
        chains = "{} chains".format(len(self.chains))
        if len(self.chains) == 1: chains = chains[:-1]
        ligands = "{} ligands".format(len(self.ligands))
        if len(self.ligands) == 1: ligands = ligands[:-1]
        return "<Model ({}, {})>".format(chains, ligands)



class Chain:

    def __init__(self, *residues, id=""):
        self.residues = residues
        self.id = id
    

    def __repr__(self):
        return f"<Chain {self.id} ({len(self.residues)} residues)>"



class Residue:

    def __init__(self, *atoms, id="", name=""):
        self.atoms = atoms
        self.id, self.name = id, name
    

    def __repr__(self):
        return f"<Residue {self.name} ({self.id})>"



class Ligand:

    def __init__(self, *atoms, id="", name=""):
        self.atoms = atoms
        self.id, self.name = id, name
    

    def __repr__(self):
        return f"<Ligand {self.name} ({self.id})>"



class Atom:

    def __init__(self, element, x, y, z, id, name, charge=0, bvalue=0, anisotropy=None):
        self.location = np.array([x, y, z])
        self.element, self.id, self.name = id, name, element
        self.charge, self.bvalue, self.anisotropy = charge, bvalue, anisotropy
    

    def __repr__(self):
        return f"<Atom {self.id} ({self.name})>"