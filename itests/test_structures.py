from unittest import TestCase
from atomium.structures import Model, Atom

class ModelBuildingTests(TestCase):

    def test_can_build_model(self):
        model = Model()
        self.assertEqual(model.atoms(), set())
        self.assertEqual(model.mass(), 0)

        atom1 = Atom("N", 12.0, 11.5, 1.5)
        atom2 = Atom("C", 12.5, 10, 2)

        model.add_atom(atom1)
        model.add_atom(atom2)

        self.assertEqual(model.atoms(), set([atom1, atom2]))
        self.assertAlmostEqual(model.mass(), 26, delta=0.1)



class AtomInteractionTests(TestCase):

    def test_can_get_atom_relations(self):
        atom1 = Atom("N", 12.0, 11.5, 1.5)
        atom2 = Atom("C", 12.5, 10, 2)

        self.assertAlmostEqual(atom1.distance_to(atom2), 1.658312, delta=0.0005)
