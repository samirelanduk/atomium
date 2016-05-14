from unittest import TestCase
from molecupy.parsers.pdb import pdb_from_string, get_pdb_from_file
from molecupy.parsers.pdb.pdb import Pdb

class PdbFromStringTests(TestCase):

    def test_can_get_pdb_from_string(self):
        pdb = pdb_from_string("TITLE     CRYSTAL")
        self.assertIsInstance(pdb, Pdb)
        self.assertEqual(pdb.title, "CRYSTAL")



class PdbFromFileTests(TestCase):

    def test_can_get_pdb_from_file(self):
        pdb = get_pdb_from_file("tests/parsers/pdb/test.pdb")
        self.assertIsInstance(pdb, Pdb)
        self.assertEqual(pdb.classification, "LYASE")
