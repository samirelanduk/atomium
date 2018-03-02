from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.structures.models import Model, ChainStructure, MoleculeStructure
from atomium.structures.molecules import AtomicStructure
from atomium.structures.chains import ResidueStructure
from atomium.structures.atoms import Atom

class ModelTest(TestCase):

    def setUp(self):
        self.atom1 = Mock(Atom)
        self.atom2 = Mock(Atom)
        self.atom3 = Mock(Atom)
        self.atoms = [self.atom1, self.atom2, self.atom3]
        def mock_init(obj, *args, **kwargs):
            obj._atoms = set()
            obj._id_atoms = {}
        self.mock_init = mock_init



class ModelCreationTests(ModelTest):

    @patch("atomium.structures.molecules.AtomicStructure.__init__")
    def test_model_is_atomic_structure(self, mock_init):
        mock_init.side_effect = self.mock_init
        model = Model(*self.atoms)
        self.assertIsInstance(model, AtomicStructure)
        self.assertIsInstance(model, ResidueStructure)
        self.assertIsInstance(model, ChainStructure)
        self.assertIsInstance(model, MoleculeStructure)
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



class ModelAtomAdditionTests(ModelTest):

    def test_adding_atom_updates_atom(self):
        model = Model(self.atom1, self.atom2)
        model.add_atom(self.atom3)
        self.assertIs(self.atom3._model, model)



class ModelAtomRemovalTests(ModelTest):

    def test_removing_atom_updates_atom(self):
        model = Model(self.atom1, self.atom2, self.atom3)
        model.remove_atom(self.atom3)
        self.assertIs(self.atom3._model, None)
