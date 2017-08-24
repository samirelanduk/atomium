from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.files.xyz import Xyz, xyz_from_file
from atomium.structures.models import Model

class XyzCreationTests(TestCase):

    def test_can_create_xyz(self):
        xyz = Xyz()
        self.assertEqual(xyz._model, None)
        self.assertEqual(xyz._comment, "")


    def test_can_create_xyz_with_comment(self):
        xyz = Xyz("Glucose molecule")
        self.assertEqual(xyz._model, None)
        self.assertEqual(xyz._comment, "Glucose molecule")


    def test_xyz_comment_must_be_str(self):
        with self.assertRaises(TypeError):
            Xyz(100)



class XyzReprTests(TestCase):

    def test_xyz_repr(self):
        xyz = Xyz("Glucose molecule")
        self.assertEqual(str(xyz), "<Xyz (Glucose molecule)>")



class XyzCommentTests(TestCase):

    def test_comment_property(self):
        xyz = Xyz("Glucose molecule")
        self.assertIs(xyz._comment, xyz.comment())


    def test_can_change_comment(self):
        xyz = Xyz("Glucose molecule")
        xyz.comment("Fructose molecule")
        self.assertEqual(xyz._comment, "Fructose molecule")


    def test_xyz_comment_must_be_str(self):
        xyz = Xyz("Glucose molecule")
        with self.assertRaises(TypeError):
            xyz.comment(100)



class XyzModelTests(TestCase):

    def test_model_property(self):
        xyz = Xyz("Glucose molecule")
        xyz._model = "totally a model"
        self.assertIs(xyz._model, xyz.model())


    def test_can_change_model(self):
        model = Mock(Model)
        xyz = Xyz("Glucose molecule")
        xyz.model(model)
        self.assertIs(xyz._model, model)


    def test_xyz_model_must_be_model(self):
        xyz = Xyz("Glucose molecule")
        with self.assertRaises(TypeError):
            xyz.model(100)



class XyzFromFileTests(TestCase):

    @patch("atomium.converters.strings.string_from_file")
    @patch("atomium.converters.string2xyz.string_to_xyz")
    def test_can_get_xyz_from_file(self, mock_xyz, mock_string):
        mock_string.return_value = "filestring"
        mock_xyz.return_value = "xyz"
        xyz = xyz_from_file("path")
        mock_string.assert_called_with("path")
        mock_xyz.assert_called_with("filestring")
        self.assertEqual(xyz, "xyz")



class XyzToStringTests(TestCase):

    @patch("atomium.converters.structure2xyzstring.structure_to_xyz_string")
    def test_can_get_string_from_xyz(self, mock_string):
        xyz = Xyz("Glucose molecule")
        xyz._model = Model()
        mock_string.return_value = "filecontents"
        s = xyz.to_file_string()
        mock_string.assert_called_with(xyz._model, xyz._comment)
        self.assertEqual(s, "filecontents")



class XyzToFileTests(TestCase):

    @patch("atomium.converters.strings.string_to_file")
    @patch("atomium.files.xyz.Xyz.to_file_string")
    def test_can_save_xyz_to_file(self, mock_string, mock_save):
        xyz = Xyz("Glucose molecule")
        xyz._model = Model()
        mock_string.return_value = "filestring"
        xyz.save("glucose.xyz")
        mock_save.assert_called_with("filestring", "glucose.xyz")
