from base import IntegratedTest
from atomium.structures import Model, Atom
import atomium

class StructureTests(IntegratedTest):

    def test_can_manipualte_model(self):
        # Create a model
        model = Model()
        self.assertEqual(model.atoms(), set())
        self.assertEqual(model.mass(), 0)

        # Create some atoms
        atom1 = Atom("N", 12.0, 11.5, 1.5, atom_id=1)
        atom2 = Atom("C", 12.5, 10, 2, atom_id=2)

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

        # Atom retrieval can pick by element and ID
        self.assertIs(model.atom(element="N"), atom1)
        self.assertIs(model.atom(element="C"), atom2)
        self.assertIs(model.atom(atom_id=1), atom1)
        self.assertIs(model.atom(atom_id=2), atom2)

        # Atoms can be bonded to each other
        self.assertEqual(atom1.bonded_atoms(), set())
        self.assertEqual(atom2.bonded_atoms(), set())
        atom1.bond(atom2)
        self.assertEqual(atom1.bonded_atoms(), set([atom2]))
        self.assertEqual(atom2.bonded_atoms(), set([atom1]))
        atom2.unbond(atom1)
        self.assertEqual(atom1.bonded_atoms(), set())
        self.assertEqual(atom2.bonded_atoms(), set())

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
        model.remove_atom(atom1)
        self.assertIn(atom2, model)
        self.assertNotIn(atom1, model)

        # Model can be saved and reloaded
        model.save("itests/files/model.xyz", description="Some atoms")
        new = atomium.xyz_from_file("itests/files/model.xyz")
        self.assertEqual(new.comment(), "Some atoms")
        self.assertEqual(len(new.model().atoms()), 1)
        self.assertEqual(new.model().atom().x(), 0.5)
        self.assertEqual(new.model().atom().y(), 1.5)
        self.assertEqual(new.model().atom().z(), -0.5)
