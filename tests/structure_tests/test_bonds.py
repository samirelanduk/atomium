from unittest import TestCase
import unittest.mock
from molecupy.structures import Bond, Atom
from molecupy.exceptions import LongBondWarning

class BondTest(TestCase):

    def setUp(self):
        self.atom1 = unittest.mock.Mock(Atom)
        self.atom1.atom_id.return_value = 100
        self.atom1.bonds.return_value = set()
        self.atom2 = unittest.mock.Mock(Atom)
        self.atom2.atom_id.return_value = 101
        self.atom2.bonds.return_value = set()
        self.atom1.distance_to.return_value = 2.1
        self.atom2.distance_to.return_value = 2.1



class BondCreationTests(BondTest):

    def test_can_create_covalent_bonds(self):
        bond = Bond(self.atom1, self.atom2)
        self.assertEqual(bond._atoms, set((self.atom1, self.atom2)))


    def test_bond_requires_atoms(self):
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


    def test_warning_issued_on_ludicrous_bond_length(self):
        self.atom1.distance_to.return_value = 22
        self.atom2.distance_to.return_value = 22
        with self.assertWarns(LongBondWarning):
            bond = Bond(self.atom1, self.atom2)



class BondPropertyTests(BondTest):

    def test_can_get_atoms(self):
        bond = Bond(self.atom1, self.atom2)
        self.assertEqual(bond.atoms(), set((self.atom1, self.atom2)))


    def test_atoms_not_directly_modifiable(self):
        bond = Bond(self.atom1, self.atom2)
        self.assertEqual(bond.atoms(), set((self.atom1, self.atom2)))
        bond.atoms().remove(self.atom2)
        self.assertEqual(bond.atoms(), set((self.atom1, self.atom2)))



class BondLengthTests(BondTest):

    def test_can_get_bond_length(self):
        bond = Bond(self.atom1, self.atom2)
        self.assertEqual(bond.bond_length(), 2.1)



class BondAtomTests(BondTest):

    def test_atoms_know_about_bond(self):
        bond = Bond(self.atom1, self.atom2)
        self.assertEqual(self.atom1.bonds(), set((bond,)))
        self.assertEqual(self.atom2.bonds(), set((bond,)))



class BondDeletionTests(BondTest):

    def test_can_delete_bond(self):
        bond = Bond(self.atom1, self.atom2)
        self.assertEqual(self.atom1.bonds(), set((bond,)))
        self.assertEqual(self.atom2.bonds(), set((bond,)))
        bond.delete()
        self.assertEqual(self.atom1.bonds(), set())
        self.assertEqual(self.atom2.bonds(), set())
