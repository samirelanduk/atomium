from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.structures import Molecule

class MoleculeModelTests(TestCase):

    def test_can_refer_to_model(self):
        class Structure(Molecule): pass
        structure = Structure()
        structure._model = "MODEL"
        self.assertIs(structure._model, structure.model)
