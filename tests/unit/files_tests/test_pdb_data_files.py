from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.files.pdbdatafile import PdbDataFile, pdb_data_file_from_file
from atomium.files.pdbdatafile import fetch_data_file

class PdbDataFileSlotsTests(TestCase):

    def test_data_file_slots(self):
        data_file = PdbDataFile()
        data_file.atoms = "aaa"
        data_file.heteroatoms = "aaa"
        data_file.connections = "aaa"
        with self.assertRaises(AttributeError):
            data_file.snorgleborgle = "aaa"



class PdbDataReprTests(TestCase):

    def test_data_file_repr(self):
        data_file = PdbDataFile()
        self.assertEqual(str(data_file), "<PdbDataFile>")



class PdbDataFileFromFileTests(TestCase):

    @patch("atomium.files.pdbfile.pdb_file_from_file")
    @patch("atomium.converters.pdbfile2pdbdatafile.pdb_file_to_pdb_data_file")
    def test_can_get_data_file_from_file(self, mock_data, mock_file):
        pdb_file, data_file = Mock(), Mock()
        mock_file.return_value = pdb_file
        mock_data.return_value = data_file
        returned_data_file = pdb_data_file_from_file("path")
        mock_file.assert_called_with("path")
        mock_data.assert_called_with(pdb_file)
        self.assertIs(data_file, returned_data_file)



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
        self.assertIs(data_file, returned_data_file)
