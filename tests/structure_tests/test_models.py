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



class ModelReprTests(ModelTest):

    def test_model_repr(self):
        structure = Model(self.atom1, self.atom2, self.atom3)
        self.assertEqual(str(structure), "<Model (3 atoms)>")



class ModelToStringTests(ModelTest):

    @patch("atomium.converters.model2xyzstring.model_to_xyz_string")
    def test_can_save_as_xyz_string(self, mock_convert):
        mock_convert.return_value = "filestring"
        model = Model(*self.atoms)
        s = model.to_file_string("xyz")
        self.assertEqual(s, "filestring")
        mock_convert.assert_called_with(model, "")


    @patch("atomium.converters.model2xyzstring.model_to_xyz_string")
    def test_can_save_as_xyz_string_with_comment(self, mock_convert):
        mock_convert.return_value = "filestring"
        model = Model(*self.atoms)
        s = model.to_file_string("xyz", description="A description")
        self.assertEqual(s, "filestring")
        mock_convert.assert_called_with(model, "A description")


    def test_invalid_file_format_is_error(self):
        model = Model(*self.atoms)
        with self.assertRaises(ValueError):
            model.to_file_string("nosuchfile")



class ModelSavingTests(ModelTest):

    @patch("atomium.structures.models.Model.to_file_string")
    @patch("atomium.converters.strings.string_to_file")
    def test_saving_uses_correct_functions(self, mock_save, mock_string):
        mock_string.return_value = "filestring"
        model = Model(*self.atoms)
        model.save("file.xyz", "a description")
        mock_string.assert_called_with("xyz", "a description")
        mock_save.assert_called_with("filestring", "file.xyz")


    @patch("atomium.structures.models.Model.to_file_string")
    @patch("atomium.converters.strings.string_to_file")
    def test_file_format_extracted(self, mock_save, mock_string):
        mock_string.return_value = "filestring"
        model = Model(*self.atoms)
        model.save("path/to/file.dfghdfg", "a description")
        mock_string.assert_called_with("dfghdfg", "a description")
        mock_save.assert_called_with("filestring", "path/to/file.dfghdfg")
