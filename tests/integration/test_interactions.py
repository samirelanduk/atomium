import math
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
    

    def test_atom_angles(self):
        # Angles between atoms
        atom1 = self.pdb.model.atom(1)
        atom2 = self.pdb.model.atom(2)
        atom3 = self.pdb.model.atom(3)
        self.assertEqual(round(atom1.angle(atom2, atom3), 2), 0.63)
        self.assertEqual(round(atom2.angle(atom1, atom3), 2), 1.91)
        self.assertEqual(round(atom3.angle(atom1, atom2), 2), 0.6)
        self.assertAlmostEqual(sum([
            atom1.angle(atom2, atom3),
            atom2.angle(atom1, atom3),
            atom3.angle(atom1, atom2)
        ]), math.pi, delta=0.001)

        # Angles between points
        self.assertEqual(round(atom1.angle([1, 2, 3], [4, 5, 6]), 2), 0.05)