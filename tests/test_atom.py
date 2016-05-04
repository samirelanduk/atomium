from unittest import TestCase
from molecupy import exceptions
from molecupy.atomic import Atom

class AtomTest(TestCase):

    def check_valid_atom(self, atom, check_atom_id=False, check_atom_name=False):
        self.assertIsInstance(atom, Atom)
        self.assertIsInstance(atom.x, float)
        self.assertIsInstance(atom.y, float)
        self.assertIsInstance(atom.z, float)
        self.assertIsInstance(atom.element, str)
        if check_atom_id:
            self.assertIsInstance(atom.atom_id, int)
        if check_atom_name:
            self.assertIsInstance(atom.atom_name, str)
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


    def test_can_create_atom_with_name(self):
        atom = Atom(1.0, 2.0, 3.0, "C", atom_name="CA")
        self.check_valid_atom(atom, check_atom_name=True)


    def test_atom_name_must_be_str(self):
        with self.assertRaises(TypeError):
            atom = Atom(1.0, 2.0, 3.0, "C", atom_name=1001)



class AtomBehaviorTests(AtomTest):

    def test_can_get_atom_mass(self):
        lithium = Atom(1.0, 1.0, 1.0, "Li")
        sodium = Atom(1.0, 1.0, 1.0, "Na")
        iron = Atom(1.0, 1.0, 1.0, "Fe")
        uranium = Atom(1.0, 1.0, 1.0, "U")
        self.assertAlmostEqual(lithium.get_mass(), 7, delta=0.5)
        self.assertAlmostEqual(sodium.get_mass(), 23, delta=0.5)
        self.assertAlmostEqual(iron.get_mass(), 56, delta=0.5)
        self.assertAlmostEqual(uranium.get_mass(), 238, delta=0.5)


    def test_strange_elements_have_zero_mass(self):
        mysterium = Atom(1.0, 1.0, 1.0, "My")
        self.assertEqual(mysterium.get_mass(), 0)



class AtomInteractionTests(AtomTest):

    def test_can_determine_distance_between_atoms(self):
        atom1 = Atom(-0.791, 64.789, 30.59, "O") # Atom 2621 in 1LOL
        atom2 = Atom(5.132, 63.307, 56.785, "C") # Atom 1011 in 1LOL
        pymol_calculated_distance = 26.9
        self.assertAlmostEqual(
         atom1.distance_to(atom2),
         pymol_calculated_distance,
         delta=0.01
        )
