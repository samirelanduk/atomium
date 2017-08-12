from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.parse.pdbfile import PdbRecord, PdbFile

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

    @patch("atomium.parse.pdbfile.PdbFile.length")
    def test_pdb_file_repr(self, mock_length):
        mock_length.return_value = 3
        pdb_file = PdbFile(self.records[0], self.records[1], self.records[2])
        self.assertEqual(str(pdb_file), "<PdbFile (3 records)>")



class PdbFileEqualityTests(PdbFileTest):

    @patch("atomium.parse.pdbfile.PdbFile.length")
    def test_pdb_files_of_equal_length_and_equal_records_are_equal(self, mock_len):
        mock_len.return_value = 2
        pdb_file1 = PdbFile(self.records[0], self.records[2])
        pdb_file2 = PdbFile(self.records[1], self.records[3])
        for record in self.records:
            record.__eq__ = lambda s, o: True
        self.assertEqual(pdb_file1, pdb_file2)


    @patch("atomium.parse.pdbfile.PdbFile.length")
    def test_pdb_files_of_equal_length_and_unequal_records_are_unequal(self, mock_len):
        mock_len.return_value = 2
        pdb_file1 = PdbFile(self.records[0], self.records[2])
        pdb_file2 = PdbFile(self.records[1], self.records[3])
        self.assertNotEqual(pdb_file1, pdb_file2)


    @patch("atomium.parse.pdbfile.PdbFile.length")
    def test_pdb_files_of_unequal_length_and_equal_records_are_equal(self, mock_len):
        mock_len.side_effect = (2, 1)
        pdb_file1 = PdbFile(self.records[0], self.records[2])
        pdb_file2 = PdbFile(self.records[1])
        for record in self.records:
            record.__eq__ = lambda s, o: True
        self.assertNotEqual(pdb_file1, pdb_file2)


    def test_pdb_file_not_equal_other_objects(self):
        self.assertNotEqual(PdbFile(self.records[0]), "PdbFile")



class PdbFileContainerTests(PdbFileTest):

    def test_pdb_file_container(self):
        pdb_file = PdbFile(self.records[0], self.records[2])
        self.assertIn(self.records[0], pdb_file)
        self.assertNotIn(self.records[1], pdb_file)
        self.assertIn(self.records[2], pdb_file)
        self.assertNotIn(self.records[3], pdb_file)



class PdbFileIterableTests(PdbFileTest):

    def test_can_iterate_through_pdb_file(self):
        pdb_file = PdbFile(*self.records)
        records = []
        for record in pdb_file:
            records.append(record)
        self.assertEqual(self.records, records)



class PdbFileIndexingTests(PdbFileTest):

    def test_can_get_records_by_index(self):
        pdb_file = PdbFile(*self.records)
        self.assertEqual(pdb_file[0], self.records[0])
        self.assertEqual(pdb_file[2], self.records[2])
        self.assertEqual(pdb_file[-1], self.records[-1])


    def test_can_get_records_by_slice(self):
        pdb_file = PdbFile(*self.records)
        self.assertEqual(pdb_file[0:2], self.records[0:2])
        self.assertEqual(pdb_file[2:4], self.records[2:4])
        self.assertEqual(pdb_file[-3:-1:-1], self.records[-3:-1:-1])



class PdbFileLengthTests(PdbFileTest):

    def test_can_get_length_of_pdb_file(self):
        pdb_file = PdbFile(self.records[1], self.records[2])
        self.assertEqual(pdb_file.length(), 2)
        pdb_file = PdbFile(self.records[0], self.records[1], self.records[2])
        self.assertEqual(pdb_file.length(), 3)


    @patch("atomium.parse.pdbfile.PdbFile.length")
    def test_len_uses_length(self, mock_length):
        mock_length.return_value = 3
        pdb_file = PdbFile(self.records[0], self.records[1], self.records[2])
        self.assertEqual(len(pdb_file), 3)



class PdbFileRecordTests(PdbFileTest):

    def test_can_get_records(self):
        pdb_file = PdbFile(*self.records)
        self.assertEqual(pdb_file._records, pdb_file.records())
        self.assertIsNot(pdb_file._records, pdb_file.records())



class PdbFileRecordAdditionTests(PdbFileTest):

    def test_can_add_pdb_record(self):
        pdb_file = PdbFile(*self.records[0:3])
        pdb_file.add_record(self.records[3])
        self.assertEqual(pdb_file._records, self.records)


    def test_can_only_add_pdb_records(self):
        pdb_file = PdbFile(*self.records[0:2])
        with self.assertRaises(TypeError):
            pdb_file.add_record("self.records[2]")



class PdbFileRecordRemovalTests(PdbFileTest):

    def test_can_remove_pdb_record(self):
        pdb_file = PdbFile(*self.records)
        pdb_file.remove_record(self.records[1])
        self.assertEqual(pdb_file._records, [self.records[0]] + self.records[2:])
        pdb_file.remove_record(self.records[0])
        self.assertEqual(pdb_file._records, self.records[2:])



class RecordRetrievalByNameTests(PdbFileTest):

    def test_can_get_record_by_name(self):
        self.records[0].name.return_value = "AAA"
        self.records[1].name.return_value = "BBB"
        self.records[2].name.return_value = "CCC"
        self.records[3].name.return_value = "DDD"
        pdb_file = PdbFile(*self.records)
        self.assertIs(pdb_file.get_record_by_name("AAA"), self.records[0])
        self.assertIs(pdb_file.get_record_by_name("BBB"), self.records[1])
        self.assertIs(pdb_file.get_record_by_name("CCC"), self.records[2])
        self.assertIs(pdb_file.get_record_by_name("DDD"), self.records[3])


    def test_no_match_returns_none(self):
        self.records[0].name.return_value = "AAA"
        self.records[1].name.return_value = "BBB"
        self.records[2].name.return_value = "CCC"
        self.records[3].name.return_value = "DDD"
        pdb_file = PdbFile(*self.records)
        self.assertIs(pdb_file.get_record_by_name("EEE"), None)


    def test_can_only_search_records_by_text(self):
        pdb_file = PdbFile(*self.records)
        with self.assertRaises(TypeError):
            pdb_file.get_record_by_name(100)
