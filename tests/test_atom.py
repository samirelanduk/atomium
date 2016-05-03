from unittest import TestCase
from molecupy import exceptions
from molecupy.atomic import Atom

class AtomTest(TestCase):

    def check_valid_atom(self, atom, check_atom_id=False):
        self.assertIsInstance(atom, Atom)
        self.assertIsInstance(atom.x, float)
        self.assertIsInstance(atom.y, float)
        self.assertIsInstance(atom.z, float)
        self.assertIsInstance(atom.element, str)
        if check_atom_id:
            self.assertIsInstance(atom.atom_id, int)
        self.assertRegex(str(atom), r"<Atom \([a-zA-Z]{1,2}\)>")



class AtomCreationTests(AtomTest):

    def test_can_create_atom(self):
        atom = Atom(1.0, 2.0, 3.0, "C")
        self.check_valid_atom(atom)


    def test_coordinates_must_be_floats(self):
        with self.assertRaises(TypeError):
            atom = Atom("1", 2.0, 3.0, "C")
        with self.assertRaises(TypeError):
            atom = Atom(1.0, "2", 3.0, "C")
        with self.assertRaises(TypeError):
            atom = Atom(1.0, 2.0, "3", "C")


    def test_element_must_be_element(self):
        with self.assertRaises(TypeError):
            atom = Atom(1.0, 2.0, 3.0, None)
        with self.assertRaises(exceptions.InvalidElementError):
            atom = Atom(1.0, 2.0, 3.0, "")
        with self.assertRaises(exceptions.InvalidElementError):
            atom = Atom(1.0, 2.0, 3.0, "XXX")


    def test_can_create_atom_with_id(self):
        atom = Atom(1.0, 2.0, 3.0, "C", atom_id=1001)
        self.check_valid_atom(atom, check_atom_id=True)


    def test_atom_id_must_be_int(self):
        with self.assertRaises(TypeError):
            atom = Atom(1.0, 2.0, 3.0, "C", atom_id=1.1)
        with self.assertRaises(TypeError):
            atom = Atom(1.0, 2.0, 3.0, "C", atom_id="1001")
