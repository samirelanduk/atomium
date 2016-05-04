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
