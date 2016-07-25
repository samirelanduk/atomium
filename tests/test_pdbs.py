from unittest import TestCase
import unittest.mock
from molecupy.pdb import Pdb
from molecupy.pdbdatafile import PdbDataFile

class PdbCreationTests(TestCase):

    def test_can_create_pdb(self):
        data_file = unittest.mock.Mock(spec=PdbDataFile)
        pdb = Pdb(data_file)
        self.assertIs(pdb.data_file(), data_file)
