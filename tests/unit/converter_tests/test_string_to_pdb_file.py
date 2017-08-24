from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.converters.string2pdbfile import *
from atomium.files.pdbfile import PdbFile, PdbRecord

class StringToPdbFileTests(TestCase):

    @patch("atomium.converters.string2pdbfile.string2lines")
    @patch("atomium.converters.string2pdbfile.PdbRecord")
    def test_can_convert_string_to_pdb_file(self, mock_record, mock_lines):
        mock_lines.return_value = ["Line1", "Line2", "Line3", "Line4"]
        records = [Mock(), Mock(), Mock(), Mock()]
        mock_record.side_effect = records
        string = "filestring"
        pdb_file = string_to_pdb_file(string)
        mock_lines.assert_called_with("filestring")
        mock_record.assert_any_call("Line1")
        mock_record.assert_any_call("Line2")
        mock_record.assert_any_call("Line3")
        mock_record.assert_any_call("Line4")
        self.assertIsInstance(pdb_file, PdbFile)
        self.assertEqual(pdb_file._records, records)
