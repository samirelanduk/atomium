from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.files.xyzdict2xyzstring import *

class XyzDictToXyzStringTests(TestCase):

    @patch("atomium.files.xyzdict2xyzstring.pack_header")
    @patch("atomium.files.xyzdict2xyzstring.pack_structure")
    @patch("atomium.files.utilities.lines_to_string")
    def test_can_convert_xyz_dict_to_string(self, mock_str, mock_struct, mock_head):
        xyz_dict = {"xyz": "dict"}
        mock_str.return_value = "filestring"
        filestring = xyz_dict_to_xyz_string(xyz_dict)
        mock_head.assert_called_with([], xyz_dict)
        mock_struct.assert_called_with([], xyz_dict)
        mock_str.assert_called_with([])
        self.assertEqual(filestring, "filestring")



class HeaderPackingTests(TestCase):

    def test_can_pack_header(self):
        xyz_dict = {"title": "Title!!!", "atoms": list(range(100))}
        lines = []
        pack_header(lines, xyz_dict)
        self.assertEqual(lines, ["100", "Title!!!"])


    def test_can_pack_header_no_title(self):
        xyz_dict = {"title": None, "atoms": list(range(100))}
        lines = []
        pack_header(lines, xyz_dict)
        self.assertEqual(lines, ["100"])



class StructurePackingTests(TestCase):

    def test_can_pack_structure(self):
        xyz_dict = {"atoms": [
         {"element": "C", "x": 3.45, "y": -12.42, "z": 197.25},
         {"element": "N", "x": 13.45, "y": -123.42, "z": 1.1}
        ]}
        lines = []
        pack_structure(lines, xyz_dict)
        self.assertEqual(lines, [
         "C      3.450    -12.420    197.250",
         "N     13.450   -123.420      1.100"
        ])
