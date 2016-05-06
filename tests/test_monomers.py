from unittest import TestCase
from molecupy import exceptions
from molecupy.atomic import Molecule, Atom
from molecupy.polymers import Monomer

class MonomerTest(TestCase):

    def setUp(self):
        self.atom1 = Atom(1.0, 1.0, 1.0, "H", atom_id=1, atom_name="H1")
        self.atom2 = Atom(1.0, 1.0, 2.0, "C", atom_id=2, atom_name="CA")
        self.atom3 = Atom(1.0, 1.0, 3.0, "O", atom_id=3, atom_name="OX1")
        self.atom1.covalent_bond_to(self.atom2)
        self.atom3.covalent_bond_to(self.atom2)


    def check_valid_monomer(self, monomer):
        self.assertIsInstance(monomer, Monomer)
        self.assertIsInstance(monomer, Molecule)
        for atom in monomer.atoms:
            self.assertEqual(atom.molecule, monomer)
        self.assertIsInstance(monomer.monomer_id, int)
        self.assertIsInstance(monomer.monomer_name, str)
        with self.assertRaises(AttributeError):
            monomer.molecule_id
        with self.assertRaises(AttributeError):
            monomer.molecule_name
        self.assertRegex(str(monomer), r"<Monomer \((.+)\)>")



class MonomerCreationTest(MonomerTest):

    def test_can_create_monomer(self):
        monomer = Monomer(1, "MON", self.atom1, self.atom2, self.atom3)
        self.check_valid_monomer(monomer)
