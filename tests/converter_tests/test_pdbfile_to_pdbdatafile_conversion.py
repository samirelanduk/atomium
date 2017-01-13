from unittest import TestCase
from unittest.mock import Mock
from molecupy.pdb.pdbfile import PdbFile
from molecupy.pdb.pdbdatafile import PdbDataFile
from molecupy.converters.pdbfile2pdbdatafile import pdb_data_file_from_pdb_file

class PdbFile2PdbDataFileTest(TestCase):

    def setUp(self):
        self.pdb_file = Mock(PdbFile)



class BasicPdbDataFileCreationTests(PdbFile2PdbDataFileTest):

    def test_can_create_pdb_data_file(self):
        data_file = pdb_data_file_from_pdb_file(self.pdb_file)
        self.assertIsInstance(data_file, PdbDataFile)


    def test_can_only_convert_pdb_files(self):
        with self.assertRaises(TypeError):
            pdb_data_file_from_pdb_file("PDB file")
