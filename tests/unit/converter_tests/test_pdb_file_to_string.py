from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.converters.pdbfile2pdbstring import pdb_file_to_pdb_string
from atomium.files.pdbfile import PdbFile, PdbRecord

class PdbFileToPdbStringTests(TestCase):

    def test_can_make_string_from_pdb_file(self):
        pdb_file = Mock(PdbFile)
        records = (Mock(PdbRecord), Mock(PdbRecord), Mock(PdbRecord))
        records[0].text.return_value = "Line 1"
        records[1].text.return_value = "Line 2"
        records[2].text.return_value = "Line 3"
        pdb_file.records.return_value = records
        string = pdb_file_to_pdb_string(pdb_file)
        self.assertEqual(string, "\n".join([
         "Line 1" + " " * 74, "Line 2" + " " * 74, "Line 3" + " " * 74
        ]))
