from unittest import TestCase
import unittest.mock
from molecupy.structures import CovalentBond, PdbAtom

class CovalentBondCreationTests(TestCase):

    def setUp(self):
        self.atom1 = unittest.mock.Mock(spec=PdbAtom)
        self.atom2 = unittest.mock.Mock(spec=PdbAtom)


    def test_can_create_covalent_bonds(self):
        bond = CovalentBond(self.atom1, self.atom2)
        self.assertEqual(bond._atoms, set((self.atom1, self.atom2)))


    def test_covalent_bond_requires_pdb_atoms(self):
        with self.assertRaises(TypeError):
            bond = CovalentBond(self.atom1, "atom2")


    def test_atoms_must_be_different(self):
        with self.assertRaises(ValueError):
            bond = CovalentBond(self.atom1, self.atom1)
        with self.assertRaises(ValueError):
            bond = CovalentBond(self.atom2, self.atom2)



class CovalentBondPropertyTests(TestCase):

    def setUp(self):
        self.atom1 = unittest.mock.Mock(spec=PdbAtom)
        self.atom2 = unittest.mock.Mock(spec=PdbAtom)


    def test_can_get_atoms(self):
        bond = CovalentBond(self.atom1, self.atom2)
        self.assertEqual(bond.atoms(), set((self.atom1, self.atom2)))
