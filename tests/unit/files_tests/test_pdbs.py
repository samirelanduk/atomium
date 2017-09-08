from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.files.pdb import Pdb, pdb_from_file, pdb_from_string, fetch
from atomium.structures.models import Model

class PdbCreationTests(TestCase):

    def test_can_create_pdb(self):
        pdb = Pdb()
        self.assertIsNone(pdb._model)



class PdbReprTests(TestCase):

    def test_pdb_repr(self):
        pdb = Pdb()
        self.assertEqual(str(pdb), "<Pdb>")



class PdbModelTests(TestCase):

    def test_model_property(self):
        pdb = Pdb()
        pdb._model = "snarglefargle"
        self.assertIs(pdb._model, pdb.model())


    def test_can_change_model(self):
        model = Mock(Model)
        pdb = Pdb()
        pdb._model = "snarglefargle"
        pdb.model(model)
        self.assertIs(pdb._model, model)


    def test_pdb_model_must_be_model(self):
        pdb = Pdb()
        pdb._model = "snarglefargle"
        with self.assertRaises(TypeError):
            pdb.model(100)



class PdbToStringTests(TestCase):

    @patch("atomium.converters.pdb2pdbdatafile.pdb_to_pdb_data_file")
    @patch("atomium.converters.pdbdatafile2pdbfile.pdb_data_file_to_pdb_file")
    @patch("atomium.converters.pdbfile2pdbstring.pdb_file_to_pdb_string")
    def test_can_get_string_from_pdb(self, mock_string, mock_file, mock_data):
        pdb = Pdb()
        data_file, pdb_file = Mock(), Mock()
        mock_string.return_value = "filecontents"
        mock_file.return_value = pdb_file
        mock_data.return_value = data_file
        s = pdb.to_file_string()
        mock_data.assert_called_with(pdb)
        mock_file.assert_called_with(data_file)
        mock_string.assert_called_with(pdb_file)
        self.assertEqual(s, "filecontents")



class PdbToFileTests(TestCase):

    @patch("atomium.converters.strings.string_to_file")
    @patch("atomium.files.pdb.Pdb.to_file_string")
    def test_can_save_xyz_to_file(self, mock_string, mock_save):
        pdb = Pdb()
        mock_string.return_value = "filestring"
        pdb.save("test.pdb")
        mock_save.assert_called_with("filestring", "test.pdb")



class PdbFromStringTests(TestCase):

    @patch("atomium.converters.string2pdbfile.string_to_pdb_file")
    @patch("atomium.converters.pdbfile2pdbdatafile.pdb_file_to_pdb_data_file")
    @patch("atomium.converters.pdbdatafile2pdb.pdb_data_file_to_pdb")
    def test_can_get_pdb_from_file(self, mock_pdb, mock_data, mock_file):
        pdb_file, data_file, pdb = Mock(), Mock(), Mock()
        mock_file.return_value = pdb_file
        mock_data.return_value = data_file
        mock_pdb.return_value = pdb
        returned_pdb = pdb_from_string("filestring")
        mock_file.assert_called_with("filestring")
        mock_data.assert_called_with(pdb_file)
        mock_pdb.assert_called_with(data_file)
        self.assertIs(pdb, returned_pdb)



class PdbFromFileTests(TestCase):

    @patch("atomium.converters.strings.string_from_file")
    @patch("atomium.files.pdb.pdb_from_string")
    def test_can_get_pdb_from_file(self, mock_pdb, mock_string):
        mock_string.return_value = "filestring"
        pdb = Mock()
        mock_pdb.return_value = pdb
        returned_pdb = pdb_from_file("path")
        mock_string.assert_called_with("path")
        mock_pdb.assert_called_with("filestring")
        self.assertIs(pdb, returned_pdb)



class PdbfetchingTests(TestCase):

    @patch("requests.get")
    @patch("atomium.files.pdb.pdb_from_string")
    def test_can_get_pdb_from_file(self, mock_pdb, mock_get):
        response = Mock()
        response.status_code = 200
        response.text = "filestring"
        mock_get.return_value = response
        pdb = Mock()
        mock_pdb.return_value = pdb
        returned_pdb = fetch("1XXX")
        mock_get.assert_called_with("https://files.rcsb.org/view/1XXX.pdb")
        mock_pdb.assert_called_with("filestring")
        self.assertIs(pdb, returned_pdb)


    @patch("requests.get")
    def test_can_get_none_if_no_file_found(self, mock_get):
        response = Mock()
        response.status_code = 404
        self.assertIsNone(fetch("1XXX"))


    def test_pdb_code_must_be_string(self):
        with self.assertRaises(TypeError):
            fetch(4000)


    def test_pdb_code_must_be_len_4(self):
        with self.assertRaises(ValueError):
            fetch("1xx")
