from re import S
from unittest import TestCase
import atomium

class DistanceTests(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.pdb = atomium.open("tests/integration/files/1lol.cif")


    def test_atom_distances(self):
        # Atom to atom
        atom1 = self.pdb.model.atom(1)
        atom2 = self.pdb.model.atom(2)
        self.assertEqual(round(atom1.distance_to(atom2), 3), 1.496)

        # Atom to point
        self.assertEqual(round(atom1.distance_to([1, 2, 3]), 1), 68.2)