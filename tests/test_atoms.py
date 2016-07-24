from unittest import TestCase
from molecupy.structures import Atom

class AtomCreationTests(TestCase):

    def test_can_create_atom(self):
        atom = Atom("C", 100, "CA")
        self.assertEqual(atom._element, "C")
        self.assertEqual(atom._atom_id, 100)
        self.assertEqual(atom._atom_name, "CA")
        self.assertEqual(atom._molecule, None)


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

    def setUp(self):
        self.atom = Atom("C", 100, "CA")


    def test_basic_properties(self):
        self.assertEqual(self.atom.element(), "C")
        self.assertEqual(self.atom.atom_id(), 100)
        self.assertEqual(self.atom.atom_name(), "CA")
        self.assertEqual(self.atom.molecule(), None)


    def test_basic_property_update(self):
        self.atom.element("H")
        self.assertEqual(self.atom.element(), "H")
        self.atom.atom_name("C")
        self.assertEqual(self.atom.atom_name(), "C")


    def test_element_must_be_str(self):
        with self.assertRaises(TypeError):
            self.atom.element(100)


    def test_element_must_be_correct_length(self):
        with self.assertRaises(ValueError):
            self.atom.element("")
        with self.assertRaises(ValueError):
            self.atom.element("CAC")
        self.atom.element("MG")
        self.atom.element("M")


    def test_atom_name_must_be_str(self):
        with self.assertRaises(TypeError):
            self.atom.atom_name(100)




class AtomMassTests(TestCase):

    def test_can_get_atom_mass(self):
        lithium = Atom("Li", 1, "Li")
        sodium = Atom("Na", 1, "Na")
        iron = Atom("Fe", 1, "Fe")
        uranium = Atom("U", 1, "U")
        self.assertAlmostEqual(lithium.mass(), 7, delta=0.5)
        self.assertAlmostEqual(sodium.mass(), 23, delta=0.5)
        self.assertAlmostEqual(iron.mass(), 56, delta=0.5)
        self.assertAlmostEqual(uranium.mass(), 238, delta=0.5)


    def test_strange_elements_have_zero_mass(self):
        mysterium = Atom("My", 1, "My")
        self.assertEqual(mysterium.mass(), 0)
