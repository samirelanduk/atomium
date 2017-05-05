from unittest import TestCase
from atomium.structures.atoms import Atom

class AtomCreationTests(TestCase):

    def test_can_create_atom(self):
        atom = Atom("C", 2, 3, 5)
        self.assertEqual(atom._element, "C")
        self.assertEqual(atom._x, 2)
        self.assertEqual(atom._y, 3)
        self.assertEqual(atom._z, 5)
