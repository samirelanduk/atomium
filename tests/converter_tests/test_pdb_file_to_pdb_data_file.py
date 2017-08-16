from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.converters.pdbfile2pdbdatafile import pdb_file_to_pdb_data_file
from atomium.parse.pdbfile import PdbFile, PdbRecord
from atomium.parse.pdbdatafile import PdbDataFile

class PdbFileToPdbDataFileTests(TestCase):

    def test_can_create_data_file_from_pdb_file(self):
        pdb_file = Mock(PdbFile)
        data_file = pdb_file_to_pdb_data_file(pdb_file)
        self.assertIsInstance(data_file, PdbDataFile)
