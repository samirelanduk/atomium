from unittest import TestCase
from molecupy.structures import Atom, PdbAtom

class PdbAtomCreationTests(TestCase):

    def test_can_create_pdb_atom(self):
        atom = PdbAtom(10.0, 20.0, 15.0, "C", 100, "CA")
        self.assertIsInstance(atom, Atom)
        self.assertEqual(atom._x, 10.0)
        self.assertEqual(atom._y, 20.0)
        self.assertEqual(atom._z, 15.0)
        self.assertEqual(atom._element, "C")
        self.assertEqual(atom._atom_id, 100)
        self.assertEqual(atom._atom_name, "CA")


    def test_repr(self):
        atom = PdbAtom(10.0, 20.0, 15.0, "C", 100, "CA")
        self.assertEqual(str(atom), "<PdbAtom 100 (CA)>")


    def test_coordinates_must_be_float(self):
        with self.assertRaises(TypeError):
            atom = PdbAtom("10.0", 20.0, 15.0, "C", 100, "CA")
        with self.assertRaises(TypeError):
            atom = PdbAtom(10.0, "20.0", 15.0, "C", 100, "CA")
        with self.assertRaises(TypeError):
            atom = PdbAtom(10.0, 20.0, "15.0", "C", 100, "CA")



class PdbAtomPropertyTests(TestCase):

    def setUp(self):
        self.atom = PdbAtom(10.0, 20.0, 15.0, "C", 100, "CA")


    def test_basic_properties(self):
        self.assertEqual(self.atom.x(), 10.0)
        self.assertEqual(self.atom.y(), 20.0)
        self.assertEqual(self.atom.z(), 15.0)


    def test_can_set_coordinates(self):
        self.atom.x(1000.0)
        self.assertEqual(self.atom.x(), 1000.0)
        self.atom.y(1000.0)
        self.assertEqual(self.atom.y(), 1000.0)
        self.atom.z(1000.0)
        self.assertEqual(self.atom.z(), 1000.0)


    def test_coordinates_must_be_float(self):
        with self.assertRaises(TypeError):
            self.atom.x("10.0")
        with self.assertRaises(TypeError):
            self.atom.y("10.0")
        with self.assertRaises(TypeError):
            self.atom.z("10.0")
