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



class PdbAtomDistanceTests(TestCase):

    def test_can_get_inter_atomic_distance(self):
        atom1 = PdbAtom(-0.791, 64.789, 30.59, "O", 2621, "OD1") # Atom 2621 in 1LOL
        atom2 = PdbAtom(5.132, 63.307, 56.785, "C", 1011, "CD") # Atom 1011 in 1LOL
        pymol_calculated_distance = 26.9
        self.assertAlmostEqual(
         atom1.distance_to(atom2),
         pymol_calculated_distance,
         delta=0.01
        )
