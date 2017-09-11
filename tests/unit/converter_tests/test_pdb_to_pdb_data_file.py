from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.converters.pdb2pdbdatafile import *
from atomium.files.pdbdatafile import PdbDataFile
from atomium.files.pdb import Pdb

class PdbToPdbDataFileTests(TestCase):

    @patch("atomium.converters.pdb2pdbdatafile.structure_to_pdb_data_file")
    def test_can_get_data_file_from_pdb(self, mock_data):
        data_file1, data_file2, data_file3 = [Mock(), Mock(), Mock()]
        data_file1.atoms = ["a1", "a2"]
        data_file2.atoms = ["a3", "a4"]
        data_file3.atoms = ["a5", "a6"]
        data_file1.heteroatoms = ["h1", "h2"]
        data_file2.heteroatoms = ["h3", "h4"]
        data_file3.heteroatoms = ["h5", "h6"]
        mock_data.side_effect = [data_file1, data_file2, data_file3]
        pdb = Mock(Pdb)
        pdb.models.return_value = ("model1", "model2", "model3")
        returned_data_file = pdb_to_pdb_data_file(pdb)
        self.assertIs(returned_data_file, data_file1)
        mock_data.assert_any_call("model1", model=1)
        mock_data.assert_any_call("model2", model=2)
        mock_data.assert_any_call("model3", model=3)
        self.assertEqual(
         returned_data_file.atoms, ["a1", "a2", "a3", "a4", "a5", "a6"]
        )
        self.assertEqual(
         returned_data_file.heteroatoms, ["h1", "h2", "h3", "h4", "h5", "h6"]
        )
