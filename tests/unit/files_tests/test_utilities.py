from unittest import TestCase
from unittest.mock import Mock, MagicMock, patch
from atomium.files.utilities  import *

class StringFromFileTests(TestCase):

    @patch("builtins.open")
    def test_gets_string_from_file(self, mock_open):
        open_return = MagicMock()
        mock_file = Mock()
        open_return.__enter__.return_value = mock_file
        mock_file.read.return_value = "returnstring"
        mock_open.return_value = open_return
        string = string_from_file("path/to/file")
        mock_open.assert_called_with("path/to/file")
        self.assertEqual(string, "returnstring")



class StringToLinesTests(TestCase):

    def test_can_convert_filestring_to_lines(self):
        filestring = "line1\nline2"
        lines = string_to_lines(filestring)
        self.assertEqual(lines, ["line1", "line2"])


    def test_can_handle_windows_line_endings(self):
        filestring = "line1\r\nline2"
        lines = string_to_lines(filestring)
        self.assertEqual(lines, ["line1", "line2"])


    def test_can_remove_empty_lines(self):
        filestring = "line1\n\nline2\n"
        lines = string_to_lines(filestring)
        self.assertEqual(lines, ["line1", "line2"])


    def test_can_make_lines_fixed_width(self):
        filestring = "line1\nline2"
        lines = string_to_lines(filestring, width=80)
        self.assertEqual(lines, ["line1".ljust(80), "line2".ljust(80)])



class StringFromWebServicesTests(TestCase):

    @patch("atomium.files.utilities.get")
    def test_can_fetch_string(self, mock_get):
        response = Mock()
        response.status_code = 200
        response.text = "filestring"
        mock_get.return_value = response
        returned_string = fetch_string("1XXX")
        mock_get.assert_called_with("https://files.rcsb.org/view/1xxx.pdb")
        self.assertEqual(returned_string, "filestring")


    @patch("atomium.files.utilities.get")
    def test_can_fetch_pdbe_string(self, mock_get):
        response = Mock()
        response.status_code = 200
        response.text = "filestring"
        mock_get.return_value = response
        returned_string = fetch_string("1XXX", pdbe=True)
        mock_get.assert_called_with(
         "https://www.ebi.ac.uk/pdbe/entry-files/pdb1xxx.ent"
        )
        self.assertEqual(returned_string, "filestring")


    @patch("atomium.files.utilities.get")
    def test_can_get_none_if_no_file_found(self, mock_get):
        response = Mock()
        response.status_code = 404
        self.assertIsNone(fetch_string("1XXX"))


    def test_pdb_code_must_be_string(self):
        with self.assertRaises(TypeError):
            fetch_string(4000)


    def test_pdb_code_must_be_len_4(self):
        with self.assertRaises(ValueError):
            fetch_string("1xx")



class PdbDictFromFileTests(TestCase):

    @patch("atomium.files.utilities.string_from_file")
    @patch("atomium.files.utilities.pdb_string_to_pdb_dict")
    def test_can_get_data_file_from_file(self, mock_dict, mock_str):
        mock_str.return_value = "filestring"
        mock_dict.return_value = {"pdb": "dict"}
        pdb_dict = pdb_data_from_file("path")
        mock_str.assert_called_with("path")
        mock_dict.assert_called_with("filestring")
        self.assertEqual(pdb_dict, {"pdb": "dict"})



class PdbDictFetchingTests(TestCase):

    @patch("atomium.files.utilities.fetch_string")
    @patch("atomium.files.utilities.pdb_string_to_pdb_dict")
    def test_can_get_data_file_from_file(self, mock_dict, mock_str):
        mock_str.return_value = "filestring"
        mock_dict.return_value = {"pdb": "dict"}
        pdb_dict = fetch_data("1xxx", a="blorg")
        mock_str.assert_called_with("1xxx", a="blorg")
        mock_dict.assert_called_with("filestring")
        self.assertEqual(pdb_dict, {"pdb": "dict"})



class PdbFromFileTests(TestCase):

    @patch("atomium.files.utilities.pdb_data_from_file")
    @patch("atomium.files.utilities.pdb_dict_to_pdb")
    def test_can_get_data_file_from_file(self, mock_pdb, mock_dict):
        mock_dict.return_value = {"pdb": "dict"}
        mock_pdb.return_value = "PDB"
        pdb = pdb_from_file("path")
        mock_dict.assert_called_with("path")
        mock_pdb.assert_called_with({"pdb": "dict"})
        self.assertEqual(pdb, "PDB")



class PdbFetchingTests(TestCase):

    @patch("atomium.files.utilities.fetch_data")
    @patch("atomium.files.utilities.pdb_dict_to_pdb")
    def test_can_get_data_file_from_file(self, mock_pdb, mock_dict):
        mock_dict.return_value = {"pdb": "dict"}
        mock_pdb.return_value = "PDB"
        pdb = fetch("1xxx", a="blorg")
        mock_dict.assert_called_with("1xxx", a="blorg")
        mock_pdb.assert_called_with({"pdb": "dict"})
        self.assertEqual(pdb, "PDB")



class XyzDictFromFileTests(TestCase):

    @patch("atomium.files.utilities.string_from_file")
    @patch("atomium.files.utilities.xyz_string_to_xyz_dict")
    def test_can_get_data_file_from_file(self, mock_dict, mock_str):
        mock_str.return_value = "filestring"
        mock_dict.return_value = {"xyz": "dict"}
        xyz_dict = xyz_data_from_file("path")
        mock_str.assert_called_with("path")
        mock_dict.assert_called_with("filestring")
        self.assertEqual(xyz_dict, {"xyz": "dict"})



class PdbFromFileTests(TestCase):

    @patch("atomium.files.utilities.xyz_data_from_file")
    @patch("atomium.files.utilities.xyz_dict_to_xyz")
    def test_can_get_data_file_from_file(self, mock_xyz, mock_dict):
        mock_dict.return_value = {"xyz": "dict"}
        mock_xyz.return_value = "XYZ"
        xyz = xyz_from_file("path")
        mock_dict.assert_called_with("path")
        mock_xyz.assert_called_with({"xyz": "dict"})
        self.assertEqual(xyz, "XYZ")



class LinesToStringTests(TestCase):

    def test_can_convert_lines_to_string(self):
        lines = ["line1", "line2", "line3"]
        self.assertEqual(lines_to_string(lines), "line1\nline2\nline3")



class StringToFileTests(TestCase):

    @patch("builtins.open")
    def test_saves_string_to_file(self, mock_open):
        open_return = MagicMock()
        mock_file = Mock()
        mock_write = MagicMock()
        mock_file.write = mock_write
        open_return.__enter__.return_value = mock_file
        mock_open.return_value = open_return
        string_to_file("filestring", "filename")
        mock_open.assert_called_once_with("filename", "w")
        mock_write.assert_called_once_with("filestring")
