from unittest import TestCase
from atomium.structures.atoms import Atom

class AtomCreationTests(TestCase):

    def test_can_create_atom(self):
        atom = Atom("C", 2, 3, 5)
        self.assertEqual(atom._element, "C")
        self.assertEqual(atom._x, 2)
        self.assertEqual(atom._y, 3)
        self.assertEqual(atom._z, 5)


    def test_atom_element_must_be_str(self):
        with self.assertRaises(TypeError):
            Atom(1, 2, 3, 5)


    def test_atom_element_must_be_1_or_2_chars(self):
        with self.assertRaises(ValueError):
            Atom("", 2, 3, 5)
        with self.assertRaises(ValueError):
            Atom("XXX", 2, 3, 5)
        Atom("XX", 2, 3, 5)
        Atom("X", 2, 3, 5)


    def test_atom_x_coord_must_be_number(self):
        with self.assertRaises(TypeError):
            Atom("C", "2", 3, 5)
        Atom("C", 2.5, 3, 5)


    def test_atom_y_coord_must_be_number(self):
        with self.assertRaises(TypeError):
            Atom("C", 2, "3", 5)
        Atom("C", 2, 3.5, 5)


    def test_atom_z_coord_must_be_number(self):
        with self.assertRaises(TypeError):
            Atom("C", 2, 3, "5")
        Atom("C", 2, 3, 5.5)



class AtomReprTests(TestCase):

    def test_atom_repr(self):
        atom = Atom("C", 2, 3, 5)
        self.assertEqual(str(atom), "<C Atom at (2, 3, 5)>")



class AtomElementTests(TestCase):

    def test_element_property(self):
        atom = Atom("C", 2, 3, 5)
        self.assertIs(atom._element, atom.element())


    def test_can_update_element(self):
        atom = Atom("C", 2, 3, 5)
        atom.element("N")
        self.assertEqual(atom._element, "N")


    def test_atom_element_must_be_str(self):
        atom = Atom("C", 2, 3, 5)
        with self.assertRaises(TypeError):
            atom.element(1)


    def test_atom_element_must_be_1_or_2_chars(self):
        atom = Atom("C", 2, 3, 5)
        with self.assertRaises(ValueError):
            atom.element("")
        with self.assertRaises(ValueError):
            atom.element("XXX")
        atom.element("XX")



class AtomXTests(TestCase):

    def test_x_property(self):
        atom = Atom("C", 2, 3, 5)
        self.assertIs(atom._x, atom.x())


    def test_can_update_x(self):
        atom = Atom("C", 2, 3, 5)
        atom.x(6)
        self.assertEqual(atom._x, 6)


    def test_atom_x_must_be_numeric(self):
        atom = Atom("C", 2, 3, 5)
        with self.assertRaises(TypeError):
            atom.x("4")
        atom.x(4.5)



class AtomYTests(TestCase):

    def test_y_property(self):
        atom = Atom("C", 2, 3, 5)
        self.assertIs(atom._y, atom.y())


    def test_can_update_y(self):
        atom = Atom("C", 2, 3, 5)
        atom.y(6)
        self.assertEqual(atom._y, 6)


    def test_atom_y_must_be_numeric(self):
        atom = Atom("C", 2, 3, 5)
        with self.assertRaises(TypeError):
            atom.y("4")
        atom.y(4.5)



class AtomZTests(TestCase):

    def test_z_property(self):
        atom = Atom("C", 2, 3, 5)
        self.assertIs(atom._z, atom.z())


    def test_can_update_z(self):
        atom = Atom("C", 2, 3, 5)
        atom.z(6)
        self.assertEqual(atom._z, 6)


    def test_atom_z_must_be_numeric(self):
        atom = Atom("C", 2, 3, 5)
        with self.assertRaises(TypeError):
            atom.z("4")
        atom.z(4.5)



class AtomMassTests(TestCase):

    def test_known_element_mass(self):
        atom = Atom("C", 2, 3, 5)
        self.assertAlmostEqual(atom.mass(), 12, delta=0.1)
        atom._element = "H"
        self.assertAlmostEqual(atom.mass(), 1, delta=0.1)


    def test_atom_mass_case_insensitive(self):
        atom = Atom("he", 2, 3, 5)
        self.assertAlmostEqual(atom.mass(), 4, delta=0.1)
        atom = Atom("He", 2, 3, 5)
        self.assertAlmostEqual(atom.mass(), 4, delta=0.1)
        atom = Atom("hE", 2, 3, 5)
        self.assertAlmostEqual(atom.mass(), 4, delta=0.1)
        atom = Atom("HE", 2, 3, 5)
        self.assertAlmostEqual(atom.mass(), 4, delta=0.1)


    def test_unknown_atom_mass(self):
        atom = Atom("XX", 2, 3, 5)
        self.assertEqual(atom.mass(), 0)



class AtomDistanceToTests(TestCase):

    def test_can_get_distance_between_atoms(self):
        atom1 = Atom("C", 4, 8, 3)
        atom2 = Atom("H", 2, 3, 5)
        self.assertAlmostEqual(atom1.distance_to(atom2), 5.744, delta=0.001)


    def test_atom_distance_can_be_zero(self):
        atom1 = Atom("C", 4, 8, 3)
        atom2 = Atom("H", 4, 8, 3)
        self.assertEqual(atom1.distance_to(atom2), 0)


    def test_other_atom_must_be_atom(self):
        atom1 = Atom("C", 4, 8, 3)
        atom2 = "atom"
        with self.assertRaises(TypeError):
            atom1.distance_to(atom2)


    def test_can_get_distance_to_xyz_tuple(self):
        atom1 = Atom("C", 4, 8, 3)
        atom2 = (2, 3, 5)
        self.assertAlmostEqual(atom1.distance_to(atom2), 5.744, delta=0.001)
