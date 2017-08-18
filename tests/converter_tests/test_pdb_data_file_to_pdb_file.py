from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.converters.pdbdatafile2pdbfile import *

class AtomDictToAtomRecordTests(TestCase):

    def setUp(self):
        self.atom_dict = {
         "atom_id": 107, "atom_name": "N", "alt_loc": None,
         "residue_name": "GLY",
         "chain_id": "A", "residue_id": 13, "insert_code": "A",
         "x": 12.681, "y": 37.302, "z": -25.211,
         "occupancy": 1.0, "temperature_factor": 15.56,
         "element": "N", "charge": -2
        }


    @patch("atomium.converters.pdbdatafile2pdbfile.PdbRecord")
    def test_can_convert_atom_dict(self, mock_record):
        record = Mock()
        mock_record.return_value = record
        converted_record = atom_dict_to_record(self.atom_dict)
        mock_record.assert_called_with("ATOM    107  N   GLY A  13A     " +
         "12.681  37.302 -25.211 1.000 15.56           N-2")
        self.assertIs(record, converted_record)


    @patch("atomium.converters.pdbdatafile2pdbfile.PdbRecord")
    def test_can_handle_different_atom_name_sizes(self, mock_record):
        self.atom_dict["atom_name"] = "CA"
        atom_dict_to_record(self.atom_dict)
        mock_record.assert_called_with("ATOM    107  CA  GLY A  13A     " +
         "12.681  37.302 -25.211 1.000 15.56           N-2")
        self.atom_dict["atom_name"] = "CG2"
        atom_dict_to_record(self.atom_dict)
        mock_record.assert_called_with("ATOM    107  CG2 GLY A  13A     " +
         "12.681  37.302 -25.211 1.000 15.56           N-2")
        self.atom_dict["atom_name"] = "HCG2"
        atom_dict_to_record(self.atom_dict)
        mock_record.assert_called_with("ATOM    107 HCG2 GLY A  13A     " +
         "12.681  37.302 -25.211 1.000 15.56           N-2")


    @patch("atomium.converters.pdbdatafile2pdbfile.PdbRecord")
    def test_can_handle_empty_dict(self, mock_record):
        atom_dict_to_record({})
        mock_record.assert_called_with("ATOM" + " " * 51 + "1.000" + " " * 19 + "0")


    @patch("atomium.converters.pdbdatafile2pdbfile.PdbRecord")
    def test_can_convert_heteroatom_dict(self, mock_record):
        converted_record = atom_dict_to_record(self.atom_dict, hetero=True)
        mock_record.assert_called_with("HETATM  107  N   GLY A  13A     " +
         "12.681  37.302 -25.211 1.000 15.56           N-2")
