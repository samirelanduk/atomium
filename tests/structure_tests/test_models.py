from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.structures.models import Model
from atomium.structures.molecules import AtomicStructure
from atomium.structures.atoms import Atom

class ModelTest(TestCase):

    def setUp(self):
        self.atom1 = Mock(Atom)
        self.atom2 = Mock(Atom)
        self.atom3 = Mock(Atom)
        self.atoms = [self.atom1, self.atom2, self.atom3]



class ModelCreationTests(ModelTest):

    @patch("atomium.structures.molecules.AtomicStructure.__init__")
    def test_model_is_atomic_structure(self, mock_init):
        model = Model(*self.atoms)
        self.assertIsInstance(model, AtomicStructure)
        self.assertTrue(mock_init.called)
