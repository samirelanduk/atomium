import os
from unittest import TestCase
from unittest.mock import patch
import unittest
from molecupy.access import pdb_from_string, get_pdb_from_file, get_pdb_remotely
from molecupy.pdb import Pdb
from molecupy.exceptions import InvalidPdbCodeError

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



class PdbRemoteTests(TestCase):


    @patch("requests.get")
    def test_can_get_pdb_remotely(self, mock_get):
        response = unittest.mock.Mock()
        with open("tests/1LOL.pdb") as f:
            response.text = f.read()
        response.status_code = 200
        mock_get.return_value = response
        pdb = get_pdb_remotely("1LOL")
        self.assertIsInstance(pdb, Pdb)
        self.assertEqual(pdb.keywords(), ["TIM BARREL", "LYASE"])


    @patch("requests.get")
    def test_invalid_pdb_code(self, mock_get):
        response = unittest.mock.Mock()
        response.text = "<html></html>"
        response.status_code = 404
        mock_get.return_value = response
        with self.assertRaises(InvalidPdbCodeError):
            pdb = get_pdb_remotely("1LOL")
