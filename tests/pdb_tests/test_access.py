import os
from unittest import TestCase
from unittest.mock import patch
import unittest
from molecupy.pdb.access import *
from molecupy.pdb.pdb import Pdb
from molecupy.exceptions import InvalidPdbCodeError

class PdbFromStringTests(TestCase):

    def test_can_make_pdb_from_string(self):
        text = "KEYWDS    TIM BARREL, LYASE"
        pdb = pdb_from_string(text)
        self.assertIsInstance(pdb, Pdb)
        self.assertEqual(pdb.keywords(), ["TIM BARREL", "LYASE"])


    def test_can_make_pdb_data_file_from_string(self):
        text = "KEYWDS    TIM BARREL, LYASE"
        pdb = pdb_data_file_from_string(text)
        self.assertIsInstance(pdb, PdbDataFile)
        self.assertEqual(pdb.keywords(), ["TIM BARREL", "LYASE"])


    def test_can_make_pdb_file_from_string(self):
        text = "KEYWDS    TIM BARREL, LYASE"
        pdb = pdb_file_from_string(text)
        self.assertIsInstance(pdb, PdbFile)
        self.assertEqual(
         pdb.get_record_by_name("KEYWDS").text(),
         "KEYWDS    TIM BARREL, LYASE" + (" " * 53)
        )



class PdbFromFileTests(TestCase):

    def test_can_make_pdb_from_file(self):
        pdb = get_pdb_from_file(os.path.sep.join(["tests", "pdb_tests", "1LOL.pdb"]))
        self.assertIsInstance(pdb, Pdb)
        self.assertEqual(pdb.keywords(), ["TIM BARREL", "LYASE"])


    def test_can_make_pdb_data_file_from_file(self):
        pdb = get_pdb_from_file(
         os.path.sep.join(["tests", "pdb_tests", "1LOL.pdb"]),
         processing="datafile"
        )
        self.assertIsInstance(pdb, PdbDataFile)
        self.assertEqual(pdb.keywords(), ["TIM BARREL", "LYASE"])


    def test_can_make_pdb_file_from_file(self):
        pdb = get_pdb_from_file(
         os.path.sep.join(["tests", "pdb_tests", "1LOL.pdb"]),
         processing="pdbfile"
        )
        self.assertIsInstance(pdb, PdbFile)
        self.assertEqual(
         pdb.get_record_by_name("KEYWDS").text(),
         "KEYWDS    TIM BARREL, LYASE" + (" " * 53)
        )


    def test_processing_arg_must_be_valid(self):
        with self.assertRaises(ValueError):
            pdb = get_pdb_from_file(
             os.path.sep.join(["tests", "pdb_tests", "1LOL.pdb"]),
             processing="....."
            )



class PdbRemoteTests(TestCase):

    @patch("requests.get")
    def test_can_get_pdb_remotely(self, mock_get):
        response = unittest.mock.Mock()
        with open(os.path.sep.join(["tests", "pdb_tests", "1LOL.pdb"])) as f:
            response.text = f.read()
        response.status_code = 200
        mock_get.return_value = response
        pdb = get_pdb_remotely("1LOL")
        self.assertIsInstance(pdb, Pdb)
        self.assertEqual(pdb.keywords(), ["TIM BARREL", "LYASE"])


    @patch("requests.get")
    def test_can_get_pdb_data_file_remotely(self, mock_get):
        response = unittest.mock.Mock()
        with open(os.path.sep.join(["tests", "pdb_tests", "1LOL.pdb"])) as f:
            response.text = f.read()
        response.status_code = 200
        mock_get.return_value = response
        pdb = get_pdb_remotely("1LOL", processing="datafile")
        self.assertIsInstance(pdb, PdbDataFile)
        self.assertEqual(pdb.keywords(), ["TIM BARREL", "LYASE"])


    @patch("requests.get")
    def test_can_get_pdb_file_remotely(self, mock_get):
        response = unittest.mock.Mock()
        with open(os.path.sep.join(["tests", "pdb_tests", "1LOL.pdb"])) as f:
            response.text = f.read()
        response.status_code = 200
        mock_get.return_value = response
        pdb = get_pdb_remotely("1LOL", processing="pdbfile")
        self.assertIsInstance(pdb, PdbFile)
        self.assertEqual(
         pdb.get_record_by_name("KEYWDS").text(),
         "KEYWDS    TIM BARREL, LYASE" + (" " * 53)
        )


    def test_processing_arg_must_be_valid(self):
        with self.assertRaises(ValueError):
            pdb = get_pdb_remotely("1LOL", processing=".....")


    @patch("requests.get")
    def test_invalid_pdb_code(self, mock_get):
        response = unittest.mock.Mock()
        response.text = "<html></html>"
        response.status_code = 404
        mock_get.return_value = response
        with self.assertRaises(InvalidPdbCodeError):
            pdb = get_pdb_remotely("1LOL")
