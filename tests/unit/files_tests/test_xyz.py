from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.files.xyz import Xyz
from atomium.models.molecules import Model

class XyzCreationTests(TestCase):

    def test_can_create_xyz(self):
        xyz = Xyz()
        self.assertEqual(xyz._model, None)
        self.assertEqual(xyz._title, "")


    def test_can_create_xyz_with_title(self):
        xyz = Xyz("Glucose molecule")
        self.assertEqual(xyz._model, None)
        self.assertEqual(xyz._title, "Glucose molecule")



class XyzReprTests(TestCase):

    def test_xyz_repr(self):
        xyz = Xyz("Glucose molecule")
        self.assertEqual(str(xyz), "<Xyz (Glucose molecule)>")


    def test_xyz_repr_no_title(self):
        xyz = Xyz()
        self.assertEqual(str(xyz), "<Xyz>")



class XyzTitleTests(TestCase):

    def test_title_property(self):
        xyz = Xyz("Glucose molecule")
        self.assertIs(xyz._title, xyz.title)


    def test_can_change_title(self):
        xyz = Xyz("Glucose molecule")
        xyz.title = "Fructose molecule"
        self.assertEqual(xyz._title, "Fructose molecule")



class XyzModelTests(TestCase):

    def test_model_property(self):
        xyz = Xyz("Glucose molecule")
        xyz._model = "totally a model"
        self.assertIs(xyz._model, xyz.model)



class XyzToStringTests(TestCase):

    @patch("atomium.files.xyz2xyzdict.xyz_to_xyz_dict")
    @patch("atomium.files.xyzdict2xyzstring.xyz_dict_to_xyz_string")
    def test_can_get_string_from_xyz(self, mock_string, mock_dict):
        xyz = Xyz()
        xyz_dict = Mock()
        mock_string.return_value = "filecontents"
        mock_dict.return_value = xyz_dict
        s = xyz.to_file_string()
        mock_dict.assert_called_with(xyz)
        mock_string.assert_called_with(xyz_dict)
        self.assertEqual(s, "filecontents")



class XyzToFileTests(TestCase):

    @patch("atomium.files.utilities.string_to_file")
    @patch("atomium.files.xyz.Xyz.to_file_string")
    def test_can_save_xyz_to_file(self, mock_string, mock_save):
        xyz = Xyz()
        mock_string.return_value = "filestring"
        xyz.save("test.xyz")
        mock_save.assert_called_with("filestring", "test.xyz")
