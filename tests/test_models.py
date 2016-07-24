from unittest import TestCase
import unittest.mock
from molecupy.structures import Model, AtomicStructure, SmallMolecule

class ModelTest(TestCase):

    def setUp(self):
        self.small_molecule1 = unittest.mock.Mock(spec=SmallMolecule)
        self.small_molecule1._model = None
        self.small_molecule2 = unittest.mock.Mock(spec=SmallMolecule)
        self.small_molecule2._model = None



class ModelCreationTest(ModelTest):

    def test_can_create_chain(self):
        model = Model()
        self.assertIsInstance(model, AtomicStructure)
        self.assertEqual(model._atoms, set())


    def test_model_repr(self):
        model = Model()
        self.assertEqual(str(model), "<Model (0 atoms)>")



class ModelSmallMoleculeTests(ModelTest):

    def test_can_add_small_molecules(self):
        model = Model()
        self.assertEqual(model.small_molecules(), set())
        model.add_small_molecule(self.small_molecule1)
        self.assertEqual(model.small_molecules(), set([self.small_molecule1]))
        model.add_small_molecule(self.small_molecule2)
        self.assertEqual(
         model.small_molecules(),
         set([self.small_molecule1, self.small_molecule2])
        )


    def test_must_use_method_to_add_small_molecule(self):
        model = Model()
        self.assertEqual(model.small_molecules(), set())
        model.small_molecules().add(self.small_molecule1)
        self.assertEqual(model.small_molecules(), set())


    def test_can_remove_small_molecules(self):
        model = Model()
        model.add_small_molecule(self.small_molecule1)
        self.assertEqual(model.small_molecules(), set([self.small_molecule1]))
        model.remove_small_molecule(self.small_molecule1)
        self.assertEqual(model.small_molecules(), set())


    def test_small_molecule_knows_about_model(self):
        model = Model()
        self.assertIs(self.small_molecule1._model, None)
        model.add_small_molecule(self.small_molecule1)
        self.assertIs(self.small_molecule1._model, model)
        model.remove_small_molecule(self.small_molecule1)
        self.assertIs(self.small_molecule1._model, None)
