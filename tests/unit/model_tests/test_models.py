from unittest import TestCase
from unittest.mock import patch, Mock, PropertyMock
from atomium.models.molecules import Model
from atomium.models.structures import AtomStructure

class ModelCreationTests(TestCase):

    def test_can_create_model(self):
        model = Model()
        self.assertIsInstance(model, AtomStructure)
