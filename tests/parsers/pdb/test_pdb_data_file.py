import datetime
from unittest import TestCase
from molecupy.parsers.pdb.pdb_file import PdbFile
from molecupy.parsers.pdb.pdb_data_file import PdbDataFile, date_from_string

class PdbDataFileTest(TestCase):

    def setUp(self):
         with open("tests/parsers/pdb/test.pdb") as f:
            pdb_file = PdbFile(f.read())
            self.data_file = PdbDataFile(pdb_file)
            self.empty_data_file = PdbDataFile(PdbFile(""))


    def check_pdb_data_file(self, data_file):
        self.assertIsInstance(data_file, PdbDataFile)
        self.assertIsInstance(data_file.pdb_file, PdbFile)
        self.assertRegex(
         str(data_file),
         r"<PdbDataFile \(([^\s]{4})\)>"
        )


    def test_can_create_pdb_data_file(self):
        self.check_pdb_data_file(self.data_file)
        self.check_pdb_data_file(self.empty_data_file)


    def test_can_get_date_from_string(self):
        self.assertEqual(
         date_from_string("01-JAN-00"),
         datetime.datetime(2000, 1, 1).date()
        )
        self.assertEqual(
         date_from_string("28-SEP-99"),
         datetime.datetime(1999, 9, 28).date()
        )



class TitleSectionTest(PdbDataFileTest):

    def test_header(self):
        self.assertEqual(self.data_file.classification, "LYASE")
        self.assertEqual(self.empty_data_file.classification, None)
        self.assertEqual(
         self.data_file.deposition_date, datetime.datetime(1990, 9, 28).date()
        )
        self.assertEqual(self.empty_data_file.deposition_date, None)
        self.assertEqual(self.data_file.pdb_code, "1SAM")
        self.assertEqual(self.empty_data_file.pdb_code, None)
