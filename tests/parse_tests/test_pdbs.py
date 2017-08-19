from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.parse.pdb import Pdb, pdb_from_file
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



class PdbFromFileTests(TestCase):

    @patch("atomium.converters.strings.string_from_file")
    @patch("atomium.converters.string2pdbfile.string_to_pdb_file")
    @patch("atomium.converters.pdbfile2pdbdatafile.pdb_file_to_pdb_data_file")
    @patch("atomium.converters.pdbdatafile2pdb.pdb_data_file_to_pdb")
    def test_can_get_pdb_from_file(self, mock_pdb, mock_data, mock_file, mock_string):
        mock_string.return_value = "filestring"
        pdb_file, data_file, pdb = Mock(), Mock(), Mock()
        mock_file.return_value = pdb_file
        mock_data.return_value = data_file
        mock_pdb.return_value = pdb
        returned_pdb = pdb_from_file("path")
        mock_string.assert_called_with("path")
        mock_file.assert_called_with("filestring")
        mock_data.assert_called_with(pdb_file)
        mock_pdb.assert_called_with(data_file)
        self.assertIs(pdb, returned_pdb)
