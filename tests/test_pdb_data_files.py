import datetime
from unittest import TestCase
from molecupy.pdbfile import PdbFile, PdbRecord
from molecupy.pdbdatafile import PdbDataFile, date_from_string, merge_records

class PdbDataFileTest(TestCase):

    def setUp(self):
        self.empty = PdbDataFile(PdbFile(""))


    def test_has_pdb_file(self):
        self.assertIsInstance(self.empty.pdb_file, PdbFile)


    def test_repr(self):
        self.assertRegex(
         str(self.empty),
         r"<PdbDataFile \(([^\s]{4})\)>"
        )



class DateFromStringTests(TestCase):

    def test_can_get_date_from_string(self):
        self.assertEqual(
         date_from_string("01-JAN-00"),
         datetime.datetime(2000, 1, 1).date()
        )
        self.assertEqual(
         date_from_string("28-SEP-99"),
         datetime.datetime(1999, 9, 28).date()
        )



class RecordMergingTests(TestCase):

    def setUp(self):
        self.records = [PdbRecord(l, 1) for l in [
         "0123456789",
         "abcdefghij",
         "0123456789"
        ]]

        self.punc_records = [PdbRecord(l, 1) for l in [
         "0123, 456789",
         "abcd  efghij",
         "0123; 456789"
        ]]


    def test_can_merge_records(self):
        self.assertEqual(
         merge_records(self.records, 5),
         "56789 fghij 56789"
        )
        self.assertEqual(
         merge_records(self.records, 8),
         "89 ij 89"
        )


    def test_can_vary_join(self):
        self.assertEqual(
         merge_records(self.records, 5, join=""),
         "56789fghij56789"
        )
        self.assertEqual(
         merge_records(self.records, 8, join="."),
         "89.ij.89"
        )


    def test_can_condense(self):
        self.assertEqual(
         merge_records(self.punc_records, 2),
         "23,456789 cd efghij 23;456789"
        )


    def test_can_ignore_consensors(self):
        self.assertEqual(
         merge_records(self.punc_records, 2, dont_condense=","),
         "23, 456789 cd efghij 23;456789"
        )
        self.assertEqual(
         merge_records(self.punc_records, 2, dont_condense=";"),
         "23,456789 cd efghij 23; 456789"
        )
        self.assertEqual(
         merge_records(self.punc_records, 2, dont_condense=";,"),
         "23, 456789 cd efghij 23; 456789"
        )
        self.assertEqual(
         merge_records(self.punc_records, 2, dont_condense=";, "),
         "23, 456789 cd  efghij 23; 456789"
        )


class HeaderRecordTests(PdbDataFileTest):

    def test_header_processing(self):
        data_file = PdbDataFile(PdbFile(
         "HEADER    LYASE                                   06-MAY-02   1LOL"
        ))
        self.assertEqual(data_file.classification, "LYASE")
        self.assertEqual(
         data_file.deposition_date,
         datetime.datetime(2002, 5, 6).date()
        )
        self.assertEqual(data_file.pdb_code, "1LOL")


    def test_missing_header_processing(self):
        self.assertEqual(self.empty.classification, None)
        self.assertEqual(self.empty.deposition_date, None)
        self.assertEqual(self.empty.pdb_code, None)



class ObslteRecordTests(PdbDataFileTest):

    def test_obslte_processing(self):
        data_file = PdbDataFile(PdbFile(
         "OBSLTE     30-SEP-93 1LOL      1SAM"
        ))
        self.assertTrue(data_file.is_obsolete)
        self.assertEqual(
         data_file.obsolete_date,
         datetime.datetime(1993, 9, 30).date()
        )
        self.assertEqual(
         data_file.replacement_code,
         "1SAM"
        )


    def test_missing_obslte_processing(self):
        self.assertFalse(self.empty.is_obsolete)
        self.assertEqual(self.empty.obsolete_date, None)
        self.assertEqual(self.empty.replacement_code, None)



class TitleRecordTests(PdbDataFileTest):

    def test_title_processing(self):
        data_file = PdbDataFile(PdbFile(
         "TITLE     CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE\n"
         "TITLE    2 COMPLEX WITH XMP"
        ))
        self.assertEqual(
         data_file.title,
         "CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE COMPLEX WITH XMP"
        )


    def test_missing_title_processing(self):
        self.assertEqual(self.empty.title, None)
