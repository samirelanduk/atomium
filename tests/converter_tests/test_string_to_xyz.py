from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.converters.string2xyz import string_to_xyz
from atomium.xyz.xyz import Xyz

class StringToXyzTests(TestCase):

    @patch("atomium.converters.string2xyz.string2lines")
    def test_empty_string_is_blank_xyz(self, mock_lines):
        mock_lines.return_value = []
        xyz = string_to_xyz("")
        mock_lines.assert_called_with("")
        self.assertIsInstance(xyz, Xyz)
        self.assertEqual(xyz._comment, "")
        self.assertEqual(xyz._model._atoms, set())


    @patch("atomium.converters.string2xyz.string2lines")
    def test_basic_xyz_parsing(self, mock_lines):
        string = ("12\n"
        "glucose from 2gbp\n"
        "C  35.884  30.895  49.120\n"
        "C  36.177  29.853  50.124")
        mock_lines.return_value = string.split("\n")
        xyz = string_to_xyz(string)
        mock_lines.assert_called_with(string)
        self.assertEqual(xyz._comment, "glucose from 2gbp")
        atoms = xyz._model._atoms
        self.assertEqual([atom._element for atom in atoms], ["C", "C"])
        self.assertEqual(
         set(((35.884, 30.895, 49.120), (36.177, 29.853, 50.124))),
         set([(atom._x, atom._y, atom._z) for atom in atoms])
        )


    @patch("atomium.converters.string2xyz.string2lines")
    def test_can_skip_numatom_parsing(self, mock_lines):
        string = ("glucose from 2gbp\n"
        "C  35.884  30.895  49.120\n"
        "C  36.177  29.853  50.124")
        mock_lines.return_value = string.split("\n")
        xyz = string_to_xyz(string)
        mock_lines.assert_called_with(string)
        self.assertEqual(xyz._comment, "glucose from 2gbp")
        atoms = xyz._model._atoms
        self.assertEqual([atom._element for atom in atoms], ["C", "C"])
        self.assertEqual(
         set(((35.884, 30.895, 49.120), (36.177, 29.853, 50.124))),
         set([(atom._x, atom._y, atom._z) for atom in atoms])
        )


    @patch("atomium.converters.string2xyz.string2lines")
    def test_can_skip_comment(self, mock_lines):
        string = ("C  35.884  30.895  49.120\n"
        "C  36.177  29.853  50.124")
        mock_lines.return_value = string.split("\n")
        xyz = string_to_xyz(string)
        mock_lines.assert_called_with(string)
        self.assertEqual(xyz._comment, "")
        atoms = xyz._model._atoms
        self.assertEqual([atom._element for atom in atoms], ["C", "C"])
        self.assertEqual(
         set(((35.884, 30.895, 49.120), (36.177, 29.853, 50.124))),
         set([(atom._x, atom._y, atom._z) for atom in atoms])
        )
