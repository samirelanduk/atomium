import os
from unittest import TestCase
from molecupy.access import pdb_from_string, get_pdb_from_file
from molecupy.pdb import Pdb

class PdbFromStringTests(TestCase):

    def test_can_make_pdb_from_string(self):
        text = "KEYWDS    TIM BARREL, LYASE"
        pdb = pdb_from_string(text)
        self.assertIsInstance(pdb, Pdb)
        self.assertEqual(pdb.keywords(), ["TIM BARREL", "LYASE"])



class PdbFromFileTests(TestCase):

    def test_can_make_pdb_from_file(self):
        pdb = get_pdb_from_file(os.path.sep.join(["tests", "1LOL.pdb"]))
        self.assertIsInstance(pdb, Pdb)
        self.assertEqual(pdb.keywords(), ["TIM BARREL", "LYASE"])
