from unittest import TestCase
from molecupy import exceptions
from molecupy.structures import PdbAtom, AtomicStructure, PdbModel

class ModelTest(TestCase):

    def check_valid_model(self, model):
        self.assertIsInstance(model, PdbModel)
        self.assertIsInstance(model, AtomicStructure)
        self.assertIsInstance(model.atoms, set)
        self.assertRegex(str(model), r"<Model \((\d+) atoms\)>")



class ModelCreationTests(ModelTest):

    def test_can_make_model(self):
        model = PdbModel()
        self.check_valid_model(model)
