from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.pdb.pdbfile import PdbRecord, PdbFile

class PdbFileTest(TestCase):

    def setUp(self):
        self.records = [Mock(PdbRecord) for i in (1, 2, 3, 4)]



class PdbFileCreationTests(PdbFileTest):

    def test_can_create_pdb_file(self):
        pdb_file = PdbFile(self.records[0], self.records[1], self.records[2])
        self.assertEqual(pdb_file._records, self.records[:-1])


    def test_pdb_file_requires_pdb_records(self):
        with self.assertRaises(TypeError):
            PdbFile(self.records[0], "self.records[1]", self.records[2])



class PdbFileReprTests(PdbFileTest):

    @patch("atomium.pdb.pdbfile.PdbFile.length")
    def test_pdb_file_repr(self, mock_length):
        mock_length.return_value = 3
        pdb_file = PdbFile(self.records[0], self.records[1], self.records[2])
        self.assertEqual(str(pdb_file), "<PdbFile (3 records)>")



class PdbFileEqualityTests(PdbFileTest):

    @patch("atomium.pdb.pdbfile.PdbFile.length")
    def test_pdb_files_of_equal_length_and_equal_records_are_equal(self, mock_len):
        mock_len.return_value = 2
        pdb_file1 = PdbFile(self.records[0], self.records[2])
        pdb_file2 = PdbFile(self.records[1], self.records[3])
        for record in self.records:
            record.__eq__ = lambda s, o: True
        self.assertEqual(pdb_file1, pdb_file2)


    @patch("atomium.pdb.pdbfile.PdbFile.length")
    def test_pdb_files_of_equal_length_and_unequal_records_are_unequal(self, mock_len):
        mock_len.return_value = 2
        pdb_file1 = PdbFile(self.records[0], self.records[2])
        pdb_file2 = PdbFile(self.records[1], self.records[3])
        self.assertNotEqual(pdb_file1, pdb_file2)


    @patch("atomium.pdb.pdbfile.PdbFile.length")
    def test_pdb_files_of_unequal_length_and_equal_records_are_equal(self, mock_len):
        mock_len.side_effect = (2, 1)
        pdb_file1 = PdbFile(self.records[0], self.records[2])
        pdb_file2 = PdbFile(self.records[1])
        for record in self.records:
            record.__eq__ = lambda s, o: True
        self.assertNotEqual(pdb_file1, pdb_file2)


    def test_pdb_file_not_equal_other_objects(self):
        self.assertNotEqual(PdbFile(self.records[0]), "PdbFile")



class PdbFileLengthTests(PdbFileTest):

    def test_can_get_length_of_pdb_file(self):
        pdb_file = PdbFile(self.records[1], self.records[2])
        self.assertEqual(pdb_file.length(), 2)
        pdb_file = PdbFile(self.records[0], self.records[1], self.records[2])
        self.assertEqual(pdb_file.length(), 3)


    @patch("atomium.pdb.pdbfile.PdbFile.length")
    def test_len_uses_length(self, mock_length):
        mock_length.return_value = 3
        pdb_file = PdbFile(self.records[0], self.records[1], self.records[2])
        self.assertEqual(len(pdb_file), 3)
