from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.converters.pdbdatafile2model import *
from atomium.parse.pdbdatafile import PdbDataFile
from atomium.structures.models import Model

class PdbDataFileToModelTests(TestCase):

    def test_can_get_model_from_data_file(self):
        data_file = Mock(PdbDataFile)
        model = pdb_data_file_to_model(data_file)
        self.assertIsInstance(model, Model)
