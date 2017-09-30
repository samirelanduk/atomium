from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.files.xyzstring2xyzdict import *

class XyzStringToPdbDictTests(TestCase):

    @patch("atomium.files.utilities.string_to_lines")
    @patch("atomium.files.xyzstring2xyzdict.extract_header")
    @patch("atomium.files.xyzstring2xyzdict.extract_structure")
    def test_can_convert_pdb_string_to_dict(self, mock_struc, mock_head, mock_lines):
        mock_lines.return_value = ["line1", "line2"]
        xyz_dict = xyz_string_to_xyz_dict("filestring")
        mock_lines.assert_called_with("filestring")
        mock_head.assert_called_with({}, ["line1", "line2"])
        mock_struc.assert_called_with({}, ["line1", "line2"])
        self.assertEqual(xyz_dict, {})



class HeaderExtractionTests(TestCase):

    def test_can_extract_header(self):
        lines = ["2", " Title line  ",  "C 1.5", "D 4.5"]
        xyz_dict = {}
        extract_header(xyz_dict, lines)
        self.assertEqual(xyz_dict, {"title": "Title line"})
        self.assertEqual(lines, ["C 1.5", "D 4.5"])


    def test_can_miss_title(self):
        lines = ["C 1.5", "C 4.5"]
        xyz_dict = {}
        extract_header(xyz_dict, lines)
        self.assertEqual(xyz_dict, {"title": None})
        self.assertEqual(lines, ["C 1.5", "C 4.5"])



class StructureExtractionTests(TestCase):

    def test_can_extract_structure(self):
        lines = [
         "C     38.553     30.400     50.259",
         "P     33.553     30.500     0.259"
        ]
        xyz_dict = {}
        extract_structure(xyz_dict, lines)
        self.assertEqual(xyz_dict, {"atoms": [{
         "element": "C", "x": 38.553, "y": 30.4, "z": 50.259
        }, {
         "element": "P", "x": 33.553, "y": 30.5, "z": 0.259
        }]})
