from unittest import TestCase
from molecupy.structures import Atom

class AtomCreationTests(TestCase):

    def test_can_create_atom(self):
        atom = Atom("C", 100, "CA")
        self.assertEqual(atom._element, "C")
        self.assertEqual(atom._atom_id, 100)
        self.assertEqual(atom._atom_name, "CA")


    def test_repr(self):
        atom = Atom("C", 100, "CA")
        self.assertEqual(str(atom), "<Atom 100 (CA)>")


    def test_element_must_be_str(self):
        with self.assertRaises(TypeError):
            atom = Atom(9, 100, "CA")


    def test_element_must_be_correct_length(self):
        with self.assertRaises(ValueError):
            atom = Atom("", 100, "CA")
        with self.assertRaises(ValueError):
            atom = Atom("CAC", 100, "CA")
        atom = Atom("MG", 100, "CA")
        atom = Atom("M", 100, "CA")


    def test_atom_id_must_be_int(self):
        with self.assertRaises(TypeError):
            atom = Atom("C", "100", "CA")


    def test_atom_name_must_be_str(self):
        with self.assertRaises(TypeError):
            atom = Atom("C", 100, 1.5)



class AtomPropertyTests(TestCase):

    def test_basic_properties(self):
        atom = Atom("C", 100, "CA")
        self.assertEqual(atom.element(), "C")
        self.assertEqual(atom.atom_id(), 100)
        self.assertEqual(atom.atom_name(), "CA")
