import datetime
from unittest import TestCase
from molecupy.parsers.pdb.pdb_file import PdbFile
from molecupy.parsers.pdb.pdb_data_file import PdbDataFile
from molecupy.parsers.pdb.pdb import Pdb

class PdbTest(TestCase):

    def setUp(self):
        with open("tests/parsers/pdb/test.pdb") as f:
           pdb_file = PdbFile(f.read())
           data_file = PdbDataFile(pdb_file)
           self.pdb = Pdb(data_file)
           self.empty_pdb = Pdb(PdbDataFile(PdbFile("")))


    def check_valid_pdb(self, pdb):
        self.assertIsInstance(pdb, Pdb)
        self.assertIsInstance(pdb.data_file, PdbDataFile)



class PdbCreationTests(PdbTest):

    def test_can_make_pdb(self):
        self.check_valid_pdb(self.pdb)
        self.check_valid_pdb(self.empty_pdb)
