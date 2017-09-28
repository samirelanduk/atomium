from unittest import TestCase
from unittest.mock import Mock, MagicMock, patch
from atomium.files.loaders import string_from_file, fetch_string
from atomium.files.loaders import pdb_data_from_file, fetch_data
from atomium.files.loaders import pdb_from_file, fetch

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



class StringFromWebServicesTests(TestCase):

    @patch("atomium.files.loaders.get")
    def test_can_fetch_string(self, mock_get):
        response = Mock()
        response.status_code = 200
        response.text = "filestring"
        mock_get.return_value = response
        returned_string = fetch_string("1XXX")
        mock_get.assert_called_with("https://files.rcsb.org/view/1xxx.pdb")
        self.assertEqual(returned_string, "filestring")


    @patch("atomium.files.loaders.get")
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


    @patch("atomium.files.loaders.get")
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

    @patch("atomium.files.loaders.string_from_file")
    @patch("atomium.files.loaders.pdb_string_to_pdb_dict")
    def test_can_get_data_file_from_file(self, mock_dict, mock_str):
        mock_str.return_value = "filestring"
        mock_dict.return_value = {"pdb": "dict"}
        pdb_dict = pdb_data_from_file("path")
        mock_str.assert_called_with("path")
        mock_dict.assert_called_with("filestring")
        self.assertEqual(pdb_dict, {"pdb": "dict"})



class PdbDictFetchingTests(TestCase):

    @patch("atomium.files.loaders.fetch_string")
    @patch("atomium.files.loaders.pdb_string_to_pdb_dict")
    def test_can_get_data_file_from_file(self, mock_dict, mock_str):
        mock_str.return_value = "filestring"
        mock_dict.return_value = {"pdb": "dict"}
        pdb_dict = fetch_data("1xxx", a="blorg")
        mock_str.assert_called_with("1xxx", a="blorg")
        mock_dict.assert_called_with("filestring")
        self.assertEqual(pdb_dict, {"pdb": "dict"})



class PdbFromFileTests(TestCase):

    @patch("atomium.files.loaders.pdb_data_from_file")
    @patch("atomium.files.loaders.pdb_dict_to_pdb")
    def test_can_get_data_file_from_file(self, mock_pdb, mock_dict):
        mock_dict.return_value = {"pdb": "dict"}
        mock_pdb.return_value = "PDB"
        pdb = pdb_from_file("path")
        mock_dict.assert_called_with("path")
        mock_pdb.assert_called_with({"pdb": "dict"})
        self.assertEqual(pdb, "PDB")



class PdbFetchingTests(TestCase):

    @patch("atomium.files.loaders.fetch_data")
    @patch("atomium.files.loaders.pdb_dict_to_pdb")
    def test_can_get_data_file_from_file(self, mock_pdb, mock_dict):
        mock_dict.return_value = {"pdb": "dict"}
        mock_pdb.return_value = "PDB"
        pdb = fetch("1xxx", a="blorg")
        mock_dict.assert_called_with("1xxx", a="blorg")
        mock_pdb.assert_called_with({"pdb": "dict"})
        self.assertEqual(pdb, "PDB")


'''
class PdbDataFileFetchingTests(TestCase):

    @patch("atomium.files.pdbfile.fetch_file")
    @patch("atomium.converters.pdbfile2pdbdatafile.pdb_file_to_pdb_data_file")
    def test_can_fetch_data_file(self, mock_data, mock_file):
        pdb_file, data_file = Mock(), Mock()
        mock_file.return_value = pdb_file
        mock_data.return_value = data_file
        returned_data_file = fetch_data_file("1xxx", a="blorg")
        mock_file.assert_called_with("1xxx", a="blorg")
        mock_data.assert_called_with(pdb_file)
        self.assertIs(data_file, returned_data_file)'''
