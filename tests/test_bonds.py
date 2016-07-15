from unittest import TestCase
import unittest.mock
from molecupy.structures import Bond, PdbAtom

class BondCreationTests(TestCase):

    def setUp(self):
        self.atom1 = unittest.mock.Mock(spec=PdbAtom)
        self.atom1.atom_id.return_value = 100
        self.atom2 = unittest.mock.Mock(spec=PdbAtom)
        self.atom2.atom_id.return_value = 101


    def test_can_create_covalent_bonds(self):
        bond = Bond(self.atom1, self.atom2)
        self.assertEqual(bond._atoms, set((self.atom1, self.atom2)))


    def test_bond_requires_pdb_atoms(self):
        with self.assertRaises(TypeError):
            bond = Bond(self.atom1, "atom2")


    def test_atoms_must_be_different(self):
        with self.assertRaises(ValueError):
            bond = Bond(self.atom1, self.atom1)
        with self.assertRaises(ValueError):
            bond = Bond(self.atom2, self.atom2)


    def test_bond_repr(self):
        bond = Bond(self.atom1, self.atom2)
        self.assertEqual(str(bond), "<Bond between Atom 100 and Atom 101>")



class BondPropertyTests(TestCase):

    def setUp(self):
        self.atom1 = unittest.mock.Mock(spec=PdbAtom)
        self.atom2 = unittest.mock.Mock(spec=PdbAtom)


    def test_can_get_atoms(self):
        bond = Bond(self.atom1, self.atom2)
        self.assertEqual(bond.atoms(), set((self.atom1, self.atom2)))


    def test_atoms_not_modifiable(self):
        bond = Bond(self.atom1, self.atom2)
        self.assertEqual(bond.atoms(), set((self.atom1, self.atom2)))
        bond.atoms().remove(self.atom2)
        self.assertEqual(bond.atoms(), set((self.atom1, self.atom2)))
