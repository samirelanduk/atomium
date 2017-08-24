from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.converters.pdb2pdbdatafile import *
from atomium.files.pdbdatafile import PdbDataFile
from atomium.files.pdb import Pdb

class PdbToPdbDataFileTests(TestCase):

    @patch("atomium.converters.pdb2pdbdatafile.structure_to_pdb_data_file")
    def test_can_get_data_file_from_pdb(self, mock_data):
        data_file = Mock()
        mock_data.return_value = data_file
        pdb = Mock(Pdb)
        pdb.model.return_value = "model"
        returned_data_file = pdb_to_pdb_data_file(pdb)
        self.assertIs(returned_data_file, data_file)
        mock_data.assert_called_with("model")
