from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.pdb.pdbfile import PdbRecord, PdbFile

class PdbFileTest(TestCase):

    def setUp(self):
        self.records = [Mock(PdbRecord) for i in (1, 2, 3)]



class PdbFileCreationTests(PdbFileTest):

    def test_can_create_pdb_file(self):
        pdb_file = PdbFile(self.records[0], self.records[1], self.records[2])
        self.assertEqual(pdb_file._records, self.records)


    def test_pdb_file_requires_pdb_records(self):
        with self.assertRaises(TypeError):
            PdbFile(self.records[0], "self.records[1]", self.records[2])



class PdbFileReprTests(PdbFileTest):

    @patch("atomium.pdb.pdbfile.PdbFile.length")
    def test_pdb_file_repr(self, mock_length):
        mock_length.return_value = 3
        pdb_file = PdbFile(self.records[0], self.records[1], self.records[2])
        self.assertEqual(str(pdb_file), "<PdbFile (3 records)>")



class PdbFileLengthTests(PdbFileTest):

    def test_can_get_length_of_pdb_file(self):
        pdb_file = PdbFile(self.records[1], self.records[2])
        self.assertEqual(pdb_file.length(), 2)
        pdb_file = PdbFile(self.records[0], self.records[1], self.records[2])
        self.assertEqual(pdb_file.length(), 3)
