from unittest import TestCase
import unittest.mock
from molecupy.structures import CovalentBond

class CovalentBondCreationTests(TestCase):

    def test_can_create_covalent_bonds(self):
        atom1 = unittest.mock.Mock()
        atom2 = unittest.mock.Mock()
        bond = CovalentBond(atom1, atom2)
        self.assertEqual(bond._atoms, set((atom1, atom2)))
