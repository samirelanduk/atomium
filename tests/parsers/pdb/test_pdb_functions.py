from unittest import TestCase
from molecupy.parsers.pdb import pdb_from_string, get_pdb_from_file, get_pdb_remotely
from molecupy.parsers.pdb.pdb import Pdb
from molecupy import exceptions

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



class PdbFromFileTests(TestCase):

    def test_can_get_pdb_remotely(self):
        pdb = get_pdb_remotely("1NVQ")
        self.assertIsInstance(pdb, Pdb)
        self.assertEqual(pdb.pdb_code, "1NVQ")


    def test_invalid_code_raises_error(self):
        with self.assertRaises(exceptions.InvalidPdbCodeError):
            pdb = get_pdb_remotely("XXXX")
