import atomium
from tests.integration.base import IntegratedTest

class XyzReadingTests(IntegratedTest):

    def test_can_read_xyz_file(self):
        # Read file
        xyz = atomium.xyz_from_file("tests/integration/files/glucose.xyz")
        self.assertEqual(xyz.title, "glucose from 2gbp")

        # The model is correct
        model = xyz.model
        self.assertEqual(len(model.atoms()), 12)
        self.assertEqual(len(model.atoms(element="C")), 6)
        self.assertEqual(len(model.atoms(element="O")), 6)
        self.assertEqual(len(model.atoms(element_regex="[CO]")), 12)

        # It has the correct mass
        self.assertAlmostEqual(model.mass, 168, delta=0.5)

        # The atoms themselves all look right
        for atom in model.atoms():
            self.assertAlmostEqual(
             atom.mass, {"C": 12, "O": 16}[atom.element], delta=0.2
            )
            self.assertIs(atom.model, model)


    def test_can_read_xyz_data(self):
        xyz = atomium.xyz_data_from_file("tests/integration/files/glucose.xyz")
        self.assertEqual(xyz["title"], "glucose from 2gbp")



class XyzSavingTests(IntegratedTest):

    def test_can_save_xyz_file(self):
        # Open and save file
        xyz = atomium.xyz_from_file("tests/integration/files/glucose.xyz")
        xyz.save("tests/integration/files/glucose2.xyz")

        # The saved xyz is correct
        with open("tests/integration/files/glucose2.xyz") as f:
            new = [l.strip() for l in f.readlines()]
        with open("tests/integration/files/glucose.xyz") as f:
            old = [l.strip() for l in f.readlines()]
        self.assertEqual(old[:-12], new[:-12])
        self.assertEqual(set(old[-12:]), set(new[-12:]))
        new = atomium.xyz_from_file("tests/integration/files/glucose2.xyz")
        model = xyz.model
        self.assertAlmostEqual(model.mass, 168, delta=0.5)
