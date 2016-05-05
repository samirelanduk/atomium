from unittest import TestCase
from molecupy import exceptions
from molecupy.atomic import Molecule, AtomicStructure, Atom

class MoleculeTest(TestCase):

    def setUp(self):
        self.atom1 = Atom(1.0, 1.0, 1.0, "H", atom_id=1, atom_name="H1")
        self.atom2 = Atom(1.0, 1.0, 2.0, "C", atom_id=2, atom_name="CA")
        self.atom3 = Atom(1.0, 1.0, 3.0, "O", atom_id=3, atom_name="OX1")
        self.atom1.covalent_bond_to(self.atom2)
        self.atom3.covalent_bond_to(self.atom2)


    def check_valid_molecule(self, molecule):
        self.assertIsInstance(molecule, Molecule)
        self.assertIsInstance(molecule, Molecule)
        self.assertRegex(str(molecule), r"<Molecule \((\d+) atoms\)>")



class MoleculeCreationTest(MoleculeTest):

    def test_can_create_molecule(self):
        molecule = Molecule(self.atom1, self.atom2, self.atom3)
        self.check_valid_molecule(molecule)
