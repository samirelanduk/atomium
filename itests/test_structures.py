from unittest import TestCase
from atomium.structures import Model, Atom

class StructureTests(TestCase):

    def test_can_manipualte_model(self):
        # Create a model
        model = Model()
        self.assertEqual(model.atoms(), set())
        self.assertEqual(model.mass(), 0)

        # Create some atoms
        atom1 = Atom("N", 12.0, 11.5, 1.5)
        atom2 = Atom("C", 12.5, 10, 2)

        # The atoms can work out their distance to each other
        self.assertAlmostEqual(atom1.distance_to(atom2), 1.658312, delta=0.0005)

        # Give the model some atoms
        model.add_atom(atom1)
        model.add_atom(atom2)

        # The model can give a full account of itself
        self.assertEqual(model.atoms(), set([atom1, atom2]))
        self.assertAlmostEqual(model.mass(), 26, delta=0.1)
        self.assertIn(atom1, model)
        self.assertIn(atom2, model)

        # Atom retrieval can pick by element
        self.assertIs(model.atom(element="N"), atom1)
        self.assertIs(model.atom(element="C"), atom2)

        # The structure has a center of mass and radius of gyration
        self.assertAlmostEqual(model.center_of_mass()[0], 12.23, delta=0.005)
        self.assertAlmostEqual(model.center_of_mass()[1], 10.81, delta=0.005)
        self.assertAlmostEqual(model.center_of_mass()[2], 1.73, delta=0.005)
        self.assertAlmostEqual(model.radius_of_gyration(), 0.83, delta=0.005)

        # The model is fully transformable
        model.translate(-12, -11.5, -1.5)
        self.assertEqual((atom1.x(), atom1.y(), atom1.z()), (0, 0, 0))
        self.assertEqual((atom2.x(), atom2.y(), atom2.z()), (0.5, -1.5, 0.5))
        model.rotate("x", 180)
        self.assertEqual((atom1.x(), atom1.y(), atom1.z()), (0, 0, 0))

        # Atoms can be removed
        model.remove_atom(atom2)
        self.assertIn(atom1, model)
        self.assertNotIn(atom2, model)
