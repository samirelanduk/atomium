from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.converters.pdbdatafile2pdb import *
from atomium.files.pdbdatafile import PdbDataFile
from atomium.files.pdb import Pdb

class PdbDataFileToPdbTests(TestCase):

    @patch("atomium.converters.pdbdatafile2pdb.pdb_data_file_to_model")
    def test_can_get_pdb_from_data_file(self, mock_model):
        model = Mock()
        mock_model.return_value = model
        data_file = Mock(PdbDataFile)
        pdb = pdb_data_file_to_pdb(data_file)
        self.assertIsInstance(pdb, Pdb)
        self.assertIs(pdb._model, model)
        mock_model.assert_called_with(data_file)
