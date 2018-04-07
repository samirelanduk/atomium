from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.structures.models import Model
from atomium.structures.molecules import AtomicStructure
from atomium.structures.atoms import Atom

class ModelTest(TestCase):

    def setUp(self):
        self.atom1, self.atom2, self.atom3 = Mock(Atom), Mock(Atom), Mock(Atom)
        self.atoms = [self.atom1, self.atom2, self.atom3]
        def mock_init(obj, *args, **kwargs):
            obj._atoms = set(args)
        self.patch1 = patch("atomium.structures.molecules.AtomicStructure.__init__")
        self.mock_init = self.patch1.start()
        self.mock_init.side_effect = mock_init



class ModelCreationTests(ModelTest):

    @patch("atomium.structures.molecules.AtomicStructure.__init__")
    def test_model_is_atomic_structure(self, mock_init):
        mock_init.side_effect = self.mock_init
        model = Model(*self.atoms)
        self.assertIsInstance(model, AtomicStructure)
        self.assertTrue(mock_init.called)


    def test_atoms_are_linked_to_model(self):
        model = Model(*self.atoms)
        self.assertIs(self.atom1._model, model)
        self.assertIs(self.atom2._model, model)
        self.assertIs(self.atom3._model, model)



class ModelReprTests(ModelTest):

    def test_model_repr(self):
        model = Model(self.atom1, self.atom2, self.atom3)
        self.assertEqual(str(model), "<Model (3 atoms)>")
