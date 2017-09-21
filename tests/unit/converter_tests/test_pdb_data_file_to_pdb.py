from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.converters.pdbdatafile2pdb import *
from atomium.files.pdbdatafile import PdbDataFile
from atomium.files.pdb import Pdb

class PdbDataFileToPdbTests(TestCase):

    @patch("atomium.converters.pdbdatafile2pdb.pdb_data_file_to_models")
    def test_can_get_pdb_from_data_file(self, mock_models):
        models = [Mock(), Mock()]
        mock_models.return_value = models
        data_file = Mock(PdbDataFile)
        data_file.code, data_file.deposition_date, data_file.title = "C", "D", "T"
        pdb = pdb_data_file_to_pdb(data_file)
        self.assertIsInstance(pdb, Pdb)
        self.assertEqual(pdb._models, models)
        self.assertEqual(pdb._code, "C")
        self.assertEqual(pdb._deposition_date, "D")
        self.assertEqual(pdb._title, "T")
        mock_models.assert_called_with(data_file)
