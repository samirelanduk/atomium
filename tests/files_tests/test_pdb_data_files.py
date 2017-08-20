from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.files.pdbdatafile import PdbDataFile

class PdbDataFileSlotsTests(TestCase):

    def test_data_file_slots(self):
        data_file = PdbDataFile()
        data_file.atoms = "aaa"
        data_file.heteroatoms = "aaa"
        data_file.connections = "aaa"
        with self.assertRaises(AttributeError):
            data_file.snorgleborgle = "aaa"



class PdbDataReprTests(TestCase):

    def test_data_file_repr(self):
        data_file = PdbDataFile()
        self.assertEqual(str(data_file), "<PdbDataFile>")
