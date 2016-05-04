import warnings
from unittest import TestCase
from molecupy import exceptions
from molecupy.atomic import CovalentBond, Atom

class BondTest(TestCase):

    def check_valid_covalent_bond(self, covalent_bond):
        self.assertIsInstance(covalent_bond, CovalentBond)
        self.assertIsInstance(covalent_bond.atoms, set)
        self.assertEqual(len(covalent_bond.atoms), 2)
        for atom in covalent_bond.atoms:
            self.assertIsInstance(atom, Atom)
        self.assertRegex(
         str(covalent_bond),
         r"<CovalentBond \([a-zA-Z]{1,2}-[a-zA-Z]{1,2}\)>"
        )



class BondCreationTests(BondTest):

    def test_can_create_covalent_bond(self):
        atom1 = Atom(1.0, 1.0, 1.0, "H")
        atom2 = Atom(1.0, 1.0, 2.0, "H")
        covalent_bond = CovalentBond(atom1, atom2)
        self.check_valid_covalent_bond(covalent_bond)


    def test_can_only_covalently_bond_atoms(self):
        atom1 = Atom(1.0, 1.0, 1.0, "H")
        atom2 = "Hydrogen atom"
        with self.assertRaises(TypeError):
            covalent_bond = CovalentBond(atom1, atom2)


    def test_warning_on_ludicrous_bond_length(self):
        atom1 = Atom(1.0, 1.0, 1.0, "H")
        atom2 = Atom(1.0, 1.0, 10.0, "H")
        with warnings.catch_warnings(record=True) as warnings_given:
            covalent_bond = CovalentBond(atom1, atom2)
            self.assertEqual(len(warnings_given), 1)
            self.assertEqual(
             warnings_given[0].category,
             exceptions.LongBondWarning
            )



class BondBehaviorTests(BondTest):

    def test_can_get_bond_length(self):
        atom1 = Atom(-0.791, 64.789, 30.59, "O")
        atom2 = Atom(5.132, 63.307, 56.785, "C")
        covalent_bond = CovalentBond(atom1, atom2)
        self.assertEqual(
         covalent_bond.get_bond_length(),
         atom1.distance_to(atom2)
        )
