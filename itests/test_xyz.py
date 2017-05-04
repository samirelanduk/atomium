import atomium
from unittest import TestCase

class XyzReadingTests(TestCase):

    def test_can_read_xyz_file(self):
        # Read file
        xyz = atomium.get_xyz_from_file("itests/files/example.xyz")
        self.assertEqual(xyz.comment(), "glucose from 2gbp")

        # XYZ has model
        model = xyz.model()

        # The atoms are all there
        self.assertEqual(len(model.atoms()), 12)
        self.assertEqual(len(model.atoms(element="C")), 6)
        self.assertEqual(len(model.atoms(element="O")), 6)

        # It has the correct mass
        self.assertAlmostEqual(model.mass(), 168, delta=0.5)

        # The atoms themselves all look right
        for atom in model:
            self.assertAlmostEqual(
             atom.mass(), {"C": 12, "N": 14}[atom.emelent()], delta=0.2
            )
