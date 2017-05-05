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
