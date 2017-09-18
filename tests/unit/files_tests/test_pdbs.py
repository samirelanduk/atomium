from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.files.pdb import Pdb, pdb_from_file, fetch
from atomium.structures.models import Model

class PdbCreationTests(TestCase):

    def test_can_create_pdb(self):
        pdb = Pdb()
        self.assertEqual(pdb._models, [])



class PdbReprTests(TestCase):

    def test_pdb_repr_no_models(self):
        pdb = Pdb()
        self.assertEqual(str(pdb), "<Pdb (0 models)>")


    def test_pdb_repr_one_model(self):
        pdb = Pdb()
        pdb._models = ["1"]
        self.assertEqual(str(pdb), "<Pdb (1 model)>")


    def test_pdb_repr_multiple_models(self):
        pdb = Pdb()
        pdb._models = ["1", "2", "3"]
        self.assertEqual(str(pdb), "<Pdb (3 models)>")



class PdbModelsTests(TestCase):

    def test_can_get_pdb_models(self):
        pdb = Pdb()
        pdb._models = ["1", "2", "3"]
        self.assertEqual(pdb.models(), ("1", "2", "3"))



class PdbModelTests(TestCase):

    def test_model_gets_first_model(self):
        pdb = Pdb()
        pdb._models = ["1", "2", "3"]
        self.assertEqual(pdb.model(), "1")


    def test_can_get_no_model(self):
        pdb = Pdb()
        self.assertIsNone(pdb.model())



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



class PdbFromFileTests(TestCase):

    @patch("atomium.files.pdbdatafile.pdb_data_file_from_file")
    @patch("atomium.converters.pdbdatafile2pdb.pdb_data_file_to_pdb")
    def test_can_get_pdb_from_file(self, mock_pdb, mock_data):
        pdb, data_file = Mock(), Mock()
        mock_pdb.return_value = pdb
        mock_data.return_value = data_file
        returned_pdb = pdb_from_file("path")
        mock_data.assert_called_with("path")
        mock_pdb.assert_called_with(data_file)
        self.assertIs(pdb, returned_pdb)



class PdbFetchingTests(TestCase):

    @patch("atomium.files.pdbdatafile.fetch_data_file")
    @patch("atomium.converters.pdbdatafile2pdb.pdb_data_file_to_pdb")
    def test_can_fetch_pdb(self, mock_pdb, mock_data):
        pdb, data_file = Mock(), Mock()
        mock_pdb.return_value = pdb
        mock_data.return_value = data_file
        returned_pdb = fetch("1xxx", a="blorg")
        mock_data.assert_called_with("1xxx", a="blorg")
        mock_pdb.assert_called_with(data_file)
        self.assertIs(pdb, returned_pdb)
