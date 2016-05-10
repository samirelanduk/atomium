from unittest import TestCase
from molecupy import exceptions
from molecupy.molecules import Atom, Molecule, AtomicStructure, Model

class ModelTest(TestCase):

    def check_valid_model(self, model):
        self.assertIsInstance(model, Model)
        self.assertIsInstance(model, AtomicStructure)
        self.assertIsInstance(model.molecules, set)
        self.assertIsInstance(model.atoms, set)
        self.assertRegex(str(model), r"<Model \((\d+) atoms\)>")



class ModelCreationTests(ModelTest):

    def test_can_make_model(self):
        model = Model()
        self.check_valid_model(model)
