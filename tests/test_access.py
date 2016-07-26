from unittest import TestCase
from molecupy.access import pdb_from_string
from molecupy.pdb import Pdb

class PdbFromStringTests(TestCase):

    def test_can_make_pdb_from_string(self):
        text = "KEYWDS    TIM BARREL, LYASE"
        pdb = pdb_from_string(text)
        self.assertIsInstance(pdb, Pdb)
        self.assertEqual(pdb.keywords(), ["TIM BARREL", "LYASE"])
