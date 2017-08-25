import atomium
from tests.integration.base import IntegratedTest

class XyzReadingTests(IntegratedTest):

    def test_can_read_xyz_file(self):
        # Read file
        xyz = atomium.xyz_from_file("tests/integration/files/glucose.xyz")
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
        for atom in model.atoms():
            self.assertAlmostEqual(
             atom.mass(), {"C": 12, "O": 16}[atom.element()], delta=0.2
            )


        # The xyz can be saved and reloaded
        xyz.save("tests/integration/files/glucose2.xyz")
        with open("tests/integration/files/glucose2.xyz") as f:
            new = f.readlines()
        with open("tests/integration/files/glucose.xyz") as f:
            old = f.readlines()
        old[-1], new[-1] = old[-1] + "\n", new[-1] + "\n"
        self.assertEqual(old[:-12], new[:-12])
        self.assertEqual(set(old[-12:]), set(new[-12:]))
        new = atomium.xyz_from_file("tests/integration/files/glucose2.xyz")
        self.assertEqual(xyz.comment(), "glucose from 2gbp")
        model = xyz.model()
        self.assertAlmostEqual(model.mass(), 168, delta=0.5)
