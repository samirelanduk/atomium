from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.converters.model2xyzstring import model_to_xyz_string
from atomium.structures.models import Model
from atomium.structures.atoms import Atom

class ModelToXyzTest(TestCase):

    def setUp(self):
        self.model = Mock(Model)
        self.atoms = [Mock(Atom) for _ in range(5)]
        self.model.atoms.return_value = set(self.atoms)


    def test_can_get_string_from_model(self):
        self.assertIsInstance(model_to_xyz_string(self.model), str)
