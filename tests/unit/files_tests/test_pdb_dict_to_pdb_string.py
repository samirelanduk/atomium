from datetime import datetime
from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.files.pdbdict2pdbstring import *

class PdbDictToPdbStringTests(TestCase):

    @patch("atomium.files.pdbdict2pdbstring.pack_header")
    @patch("atomium.files.pdbdict2pdbstring.pack_structure")
    @patch("atomium.files.pdbdict2pdbstring.pdb_lines_to_string")
    def test_can_convert_pdb_dict_to_string(self, mock_str, mock_struct, mock_head):
        pdb_dict = {"pdb": "dict"}
        mock_str.return_value = "filestring"
        filestring = pdb_dict_to_pdb_string(pdb_dict)
        mock_head.assert_called_with([], pdb_dict)
        mock_struct.assert_called_with([], pdb_dict)
        mock_str.assert_called_with([])
        self.assertEqual(filestring, "filestring")



class HeaderPackingTests(TestCase):

    def setUp(self):
        self.pdb_dict = {
         "deposition_date": datetime(1990, 9, 1).date(),
         "code": "1XYZ", "title": "ABC" * 40
        }
        self.lines = []


    def test_can_pack_full_header(self):
        pack_header(self.lines, self.pdb_dict)
        self.assertEqual(self.lines, [
         "HEADER" + " " * 44 + "01-SEP-90   1XYZ" + " " * 14,
         "TITLE     " + "ABC" * 23 + "A",
         "TITLE    2 " + "BC" + "ABC" * 16 + " " * 19
        ])


    def test_can_pack_deposition_date(self):
        self.pdb_dict["code"], self.pdb_dict["title"] = None, None
        pack_header(self.lines, self.pdb_dict)
        self.assertEqual(self.lines, [
         "HEADER" + " " * 44 + "01-SEP-90       " + " " * 14,
        ])


    def test_can_pack_code(self):
        self.pdb_dict["deposition_date"], self.pdb_dict["title"] = None, None
        pack_header(self.lines, self.pdb_dict)
        self.assertEqual(self.lines, [
         "HEADER" + " " * 44 + "            1XYZ" + " " * 14,
        ])


    def test_can_pack_title(self):
        self.pdb_dict["deposition_date"], self.pdb_dict["code"] = None, None
        pack_header(self.lines, self.pdb_dict)
        self.assertEqual(self.lines, [
         "TITLE     " + "ABC" * 23 + "A",
         "TITLE    2 " + "BC" + "ABC" * 16 + " " * 19
        ])
