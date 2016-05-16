import math
from collections import Counter
import warnings
from .exceptions import *

class PdbAtom:

    def __init__(self, x, y, z, element, atom_id, atom_name):
        for coordinate in (x, y, z):
            if not isinstance(coordinate, float):
                raise TypeError("'%s' is not a valid coordinate" % str(coordinate))
        self.x = x
        self.y = y
        self.z = z

        if not isinstance(element, str):
            raise TypeError("'%s' is not a valid element" % str(element))
        if len(element) == 0:
            raise InvalidElementError("Atom's element can't be an empty string")
        elif len(element) > 2:
            raise InvalidElementError("'%s' is not a valid element" % element)
        self.element = element

        if not isinstance(atom_id, int):
            raise TypeError("'%s' is not a valid atom_id" % str(atom_id))
        self.atom_id = atom_id

        if not isinstance(atom_name, str):
            raise TypeError("'%s' is not a valid atom_name" % str(atom_name))
        self.atom_name = atom_name

        self.covalent_bonds = set()
        self.molecule = None


    def __repr__(self):
        return "<Atom %i (%s)>" % (self.atom_id, self.element)


    def get_mass(self):
        return PERIODIC_TABLE.get(self.element.upper(), 0)


    def distance_to(self, other_atom):
        x_sum = math.pow((other_atom.x - self.x), 2)
        y_sum = math.pow((other_atom.y - self.y), 2)
        z_sum = math.pow((other_atom.z - self.z), 2)
        distance = math.sqrt(x_sum + y_sum + z_sum)
        return distance


    def covalent_bond_to(self, other_atom):
        CovalentBond(self, other_atom)


    def get_covalent_bonded_atoms(self):
        covalent_bonded_atoms = set()
        for bond in self.covalent_bonds:
            for atom in bond.atoms:
                if atom is not self: covalent_bonded_atoms.add(atom)
        return covalent_bonded_atoms


    def get_covalent_accessible_atoms(self, already_checked=None):
        already_checked = already_checked if already_checked else set()
        already_checked.add(self)
        while len(self.get_covalent_bonded_atoms().difference(already_checked)) > 0:
            picked = list(self.get_covalent_bonded_atoms().difference(already_checked))[0]
            picked.get_covalent_accessible_atoms(already_checked)
        already_checked_copy = already_checked.copy()
        already_checked_copy.discard(self)
        return already_checked_copy



class CovalentBond:

    def __init__(self, atom1, atom2):
        if not isinstance(atom1, PdbAtom) or not isinstance(atom2, PdbAtom):
            raise TypeError(
             "Can only bond atoms, not %s to '%s'" % (str(atom1), str(atom2))
            )
        if atom1.distance_to(atom2) > 5:
            warning = "The bond between atom %s and atom %s is %.2f Angstroms" % (
             str(atom1.atom_id) if atom1.atom_id else atom1.element,
             str(atom2.atom_id) if atom2.atom_id else atom2.element,
             atom1.distance_to(atom2)
            )
            warnings.warn(warning, LongBondWarning)

        self.atoms = set((atom1, atom2))
        atom1.covalent_bonds.add(self)
        atom2.covalent_bonds.add(self)


    def __repr__(self):
        atom1_element, atom2_element = [a.element for a in self.atoms]
        return "<CovalentBond (%s-%s)>" % (atom1_element, atom2_element)


    def get_bond_length(self):
        atom1, atom2 = self.atoms
        return atom1.distance_to(atom2)



class AtomicStructure:

    def __init__(self, *atoms):
        if not all(isinstance(atom, PdbAtom) for atom in atoms):
            non_atoms = [atom for atom in atoms if not isinstance(atom, PdbAtom)]
            raise TypeError("AtomicStructure needs atoms, not '%s'" % non_atoms[0])
        if not atoms:
            raise NoAtomsError("Cannot make AtomicStructure with zero atoms")
        atom_ids = [atom.atom_id for atom in atoms]
        if len(set(atom_ids)) < len(atom_ids):
            raise DuplicateAtomIdError("Cannot make AtomicStructure with duplicate atom_ids")
        self.atoms = set(atoms)


    def __repr__(self):
        return "<AtomicStructure (%i atoms)>" % len(self.atoms)


    def __contains__(self, atom):
        return atom in self.atoms


    def get_mass(self):
        return sum([atom.get_mass() for atom in self.atoms])


    def get_atom_by_name(self, atom_name):
        if not isinstance(atom_name, str):
            raise TypeError("Atom name search must be by str, not '%s'" % str(atom_name))
        for atom in self.atoms:
            if atom.atom_name == atom_name:
                return atom


    def get_atoms_by_name(self, atom_name):
        if not isinstance(atom_name, str):
            raise TypeError("Atom name search must be by str, not '%s'" % str(atom_name))
        return set([atom for atom in self.atoms if atom.atom_name == atom_name])


    def get_atoms_by_element(self, element):
        if not isinstance(element, str):
            raise TypeError("Atom element search must be by str, not '%s'" % str(element))
        return set([atom for atom in self.atoms if atom.element == element])


    def get_atom_by_id(self, atom_id):
        if not isinstance(atom_id, int):
            raise TypeError("Atom ID search must be by int, not '%s'" % str(atom_id))
        for atom in self.atoms:
            if atom.atom_id == atom_id:
                return atom


    def get_empirical_formula(self):
        return Counter([
         atom.element for atom in self.atoms if atom.element != "H"
        ])



class PdbSmallMolecule(AtomicStructure):

    def __init__(self, molecule_id, molecule_name, *atoms):
        if not isinstance(molecule_id, str):
            raise TypeError("'%s' is not a valid molecule_id" % str(molecule_id))
        self.molecule_id = molecule_id

        if not isinstance(molecule_name, str):
            raise TypeError("'%s' is not a valid molecule_name" % str(molecule_name))
        self.molecule_name = molecule_name

        AtomicStructure.__init__(self, *atoms)
        for atom in self.atoms:
            atom.molecule = self


    def __repr__(self):
        return "<SmallMolecule (%s)>" % self.molecule_name



class PdbResidue(AtomicStructure):

    def __init__(self, residue_id, residue_name, *atoms):
        if not isinstance(residue_id, int):
            raise TypeError("'%s' is not a valid residue_id" % str(residue_id))
        self.residue_id = residue_id

        if not isinstance(residue_name, str):
            raise TypeError("'%s' is not a valid molecule_name" % str(residue_name))
        self.residue_name = residue_name

        AtomicStructure.__init__(self, *atoms)
        for atom in self.atoms:
            atom.molecule = self

        self.chain = None


    def __repr__(self):
        return "<Residue (%s)>" % self.residue_name



class ResiduicStructure(AtomicStructure):

    def __init__(self, *residues):
        if not all(isinstance(residue, PdbResidue) for residue in residues):
            non_residues = [
             residue for residue in residues if not isinstance(residue, PdbResidue)
            ]
            raise TypeError("ResiduicStructure needs residues, not '%s'" % non_residues[0])
        if not residues:
            raise NoResiduesError("Cannot make ResiduicStructure with zero residues")
        residue_ids = [residue.residue_id for residue in residues]
        if len(set(residue_ids)) < len(residue_ids):
            raise DuplicateResidueIdError("Cannot make ResiduicStructure with duplicate residue_ids")
        self.residues = set(residues)


    def __repr__(self):
        return "<ResiduicStructure (%i residues)>" % len(self.residues)


    def __getattr__(self, attribute):
        if attribute == "atoms":
            atoms = set()
            for residue in self.residues:
                atoms.update(residue.atoms)
            return atoms
        else:
            return self.__getattribute__(attribute)


    def __contains__(self, obj):
        return obj in self.residues or obj in self.atoms


    def get_residue_by_name(self, residue_name):
        if not isinstance(residue_name, str):
            raise TypeError("Residue name search must be by str, not '%s'" % str(residue_name))
        for residue in self.residues:
            if residue.residue_name == residue_name:
                return residue


    def get_residues_by_name(self, residue_name):
        if not isinstance(residue_name, str):
            raise TypeError("Residue name search must be by str, not '%s'" % str(residue_name))
        return set([residue for residue in self.residues if residue.residue_name == residue_name])


    def get_residue_by_id(self, residue_id):
        if not isinstance(residue_id, int):
            raise TypeError("Residue ID search must be by int, not '%s'" % str(residue_id))
        for residue in self.residues:
            if residue.residue_id == residue_id:
                return residue



class ResiduicSequence(ResiduicStructure):

    def __init__(self, *residues):
        ResiduicStructure.__init__(self, *residues)
        self.residues = sorted(list(self.residues), key=lambda k: k.residue_id)


    def __repr__(self):
        return "<ResiduicSequence (%i residues)>" % len(self.residues)



class PdbChain(ResiduicSequence):

    def __init__(self, chain_id, *residues):
        if not isinstance(chain_id, str):
            raise TypeError("'%s' is not a valid chain_id" % str(chain_id))
        self.chain_id = chain_id

        ResiduicSequence.__init__(self, *residues)
        for residue in self.residues:
            residue.chain = self


    def __repr__(self):
        return "<Chain %s (%i residues)>" % (self.chain_id, len(self.residues))



class PdbModel(AtomicStructure):

    def __init__(self):
        self.small_molecules = set()
        self.chains = set()


    def __repr__(self):
        return "<Model (%i atoms)>" % len(self.atoms)


    def __getattr__(self, attribute):
        if attribute == "atoms":
            atoms = set()
            for small_molecule in self.small_molecules:
                atoms.update(small_molecule.atoms)
            for chain in self.chains:
                atoms.update(chain.atoms)
            return atoms
        else:
            return self.__getattribute__(attribute)


    def add_small_molecule(self, small_molecule):
        if not(isinstance(small_molecule, PdbSmallMolecule)):
            raise TypeError("Can only add SmallMolecules with add_small_molecule()")
        existing_small_molecule_ids = [mol.molecule_id for mol in self.small_molecules]
        if small_molecule.molecule_id in existing_small_molecule_ids:
            raise DuplicateSmallMoleculeError(
             "Cannot have two small molecules with ID %s" % small_molecule.molecule_id
            )
        self.small_molecules.add(small_molecule)


    def get_small_molecule_by_id(self, molecule_id):
        if not isinstance(molecule_id, str):
            raise TypeError("Can only search small molecule IDs by str")
        for small_molecule in self.small_molecules:
            if small_molecule.molecule_id == molecule_id:
                return small_molecule


    def get_small_molecule_by_name(self, molecule_name):
        if not isinstance(molecule_name, str):
            raise TypeError("Can only search small molecule names by string")
        for small_molecule in self.small_molecules:
            if small_molecule.molecule_name == molecule_name:
                return small_molecule
        pass


    def get_small_molecules_by_name(self, molecule_name):
        if not isinstance(molecule_name, str):
            raise TypeError("Can only search small molecule names by string")
        return set([
         mol for mol in self.small_molecules if mol.molecule_name == molecule_name
        ])


    def add_chain(self, chain):
        if not(isinstance(chain, PdbChain)):
            raise TypeError("Can only add Chain with add_chain()")
        existing_chain_ids = [c.chain_id for c in self.chains]
        if chain.chain_id in existing_chain_ids:
            raise DuplicateChainError(
             "Cannot have two chains with ID %s" % chain.chain_id
            )
        self.chains.add(chain)


    def get_chain_by_id(self, chain_id):
        if not isinstance(chain_id, str):
            raise TypeError("Can only search chain IDs by str")
        for chain in self.chains:
            if chain.chain_id == chain_id:
                return chain



PERIODIC_TABLE = {
 "H": 1.0079, "HE": 4.0026, "LI": 6.941, "BE": 9.0122, "B": 10.811,
  "C": 12.0107, "N": 14.0067, "O": 15.9994, "F": 18.9984, "NE": 20.1797,
   "NA": 22.9897, "MG": 24.305, "AL": 26.9815, "SI": 28.0855, "P": 30.9738,
    "S": 32.065, "CL": 35.453, "K": 39.0983, "AR": 39.948, "CA": 40.078,
     "SC": 44.9559, "TI": 47.867, "V": 50.9415, "CR": 51.9961, "MN": 54.938,
      "FE": 55.845, "NI": 58.6934, "CO": 58.9332, "CU": 63.546, "ZN": 65.39,
       "GA": 69.723, "GE": 72.64, "AS": 74.9216, "SE": 78.96, "BR": 79.904,
        "KR": 83.8, "RB": 85.4678, "SR": 87.62, "Y": 88.9059, "ZR": 91.224,
         "NB": 92.9064, "MO": 95.94, "TC": 98, "RU": 101.07, "RH": 102.9055,
          "PD": 106.42, "AG": 107.8682, "CD": 112.411, "IN": 114.818,
           "SN": 118.71, "SB": 121.76, "I": 126.9045, "TE": 127.6,
            "XE": 131.293, "CS": 132.9055, "BA": 137.327, "LA": 138.9055,
             "CE": 140.116, "PR": 140.9077, "ND": 144.24, "PM": 145,
              "SM": 150.36, "EU": 151.964, "GD": 157.25, "TB": 158.9253,
               "DY": 162.5, "HO": 164.9303, "ER": 167.259, "TM": 168.9342,
                "YB": 173.04, "LU": 174.967, "HF": 178.49, "TA": 180.9479,
                 "W": 183.84, "RE": 186.207, "OS": 190.23, "IR": 192.217,
                  "PT": 195.078, "AU": 196.9665, "HG": 200.59, "TL": 204.3833,
                   "PB": 207.2, "BI": 208.9804, "PO": 209, "AT": 210, "RN": 222,
                    "FR": 223, "RA": 226, "AC": 227, "PA": 231.0359,
                     "TH": 232.0381, "NP": 237, "U": 238.0289, "AM": 243,
                      "PU": 244, "CM": 247, "BK": 247, "CF": 251, "ES": 252,
                       "FM": 257, "MD": 258, "NO": 259, "RF": 261, "LR": 262,
                        "DB": 262, "BH": 264, "SG": 266, "MT": 268, "RG": 272,
                         "HS": 277
}
