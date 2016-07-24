import datetime
from unittest import TestCase
from molecupy.pdbfile import PdbFile, PdbRecord
from molecupy.pdbdatafile import PdbDataFile, date_from_string

class PdbDataFileTest(TestCase):

    def setUp(self):
        self.empty = PdbDataFile(PdbFile(""))



class PdbdataFilePropertiesTests(PdbDataFileTest):

    def test_has_pdb_file(self):
        self.assertIsInstance(self.empty.pdb_file(), PdbFile)


    def test_repr(self):
        self.assertEqual(str(self.empty), "<PdbDataFile (????)>")



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
