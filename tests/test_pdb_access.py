from unittest import TestCase
from molecupy import pdb_from_string, get_pdb_remotely, get_pdb_from_file
from molecupy.pdb import Pdb
from molecupy.exceptions import *

class PdbFromStringTests(TestCase):

    def test_can_get_pdb_from_string(self):
        pdb = pdb_from_string("TITLE     CRYSTAL")
        self.assertIsInstance(pdb, Pdb)
        self.assertEqual(pdb.title, "CRYSTAL")



class PdbFromFileTests(TestCase):

    def test_can_get_pdb_from_file(self):
        pdb = get_pdb_from_file("tests/1AOI.pdb")
        self.assertIsInstance(pdb, Pdb)
        self.assertEqual(pdb.classification, "DNA BINDING PROTEIN/DNA")



class PdbFromRemotelyTests(TestCase):

    def test_can_get_pdb_remotely(self):
        pdb = get_pdb_remotely("1NVQ")
        self.assertIsInstance(pdb, Pdb)
        self.assertEqual(pdb.pdb_code, "1NVQ")


    def test_invalid_code_raises_error(self):
        with self.assertRaises(InvalidPdbCodeError):
            pdb = get_pdb_remotely("XXXX")
