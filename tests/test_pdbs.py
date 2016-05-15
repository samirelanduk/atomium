import datetime
from unittest import TestCase
from molecupy.pdbfile import PdbFile
from molecupy.pdbdatafile import PdbDataFile
from molecupy.pdb import Pdb

class PdbTest(TestCase):

    def setUp(self):
        self.empty = Pdb(PdbDataFile(PdbFile("")))



class PdbPropertiesTests(PdbTest):

    def test_has_pdb_file(self):
        self.assertIsInstance(self.empty.data_file, PdbDataFile)


    def test_repr(self):
        self.assertRegex(
         str(self.empty),
         r"<Pdb \(([^\s]{4})\)>"
        )
