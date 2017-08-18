from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.converters.pdbdatafile2pdbfile import *
from atomium.parse.pdbfile import PdbFile
from atomium.parse.pdbdatafile import PdbDataFile

class StructurePackingTests(TestCase):

    @patch("atomium.converters.pdbdatafile2pdbfile.atom_dict_to_record")
    @patch("atomium.converters.pdbdatafile2pdbfile.conections_list_to_records")
    def test_can_pack_structure(self, mock_connections, mock_atom):
        mock_connections.return_value = ["c1", "c2"]
        mock_atom.side_effect = ["a1", "a2", "h1", "h2"]
        pdb_file = Mock(PdbFile)
        pdb_file._records = []
        data_file = Mock(PdbDataFile)
        data_file.connections = "ccc"
        data_file.atoms = ["1", "2"]
        data_file.heteroatoms = ["3", "4"]
        pack_structure(data_file, pdb_file)
        mock_connections.assert_called_with("ccc")
        mock_atom.assert_any_call("1")
        mock_atom.assert_any_call("2")
        mock_atom.assert_any_call("3", hetero=True)
        mock_atom.assert_any_call("4", hetero=True)
        self.assertEqual(pdb_file._records, ["a1", "a2", "h1", "h2", "c1", "c2"])



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



class ConectListToRecordsTests(TestCase):

    @patch("atomium.converters.pdbdatafile2pdbfile.PdbRecord")
    def test_can_convert_connections_list_to_records(self, mock_record):
        records = [Mock(), Mock(), Mock()]
        mock_record.side_effect = records
        connections = [{
         "atom": 1179, "bond_to": [746, 1184, 1195, 1203, 1211, 1222]
        }, {
         "atom": 1221, "bond_to": [544, 1017, 1020, 1022]
        }]
        converted_records = conections_list_to_records(connections)
        args1, kwargs = mock_record.call_args_list[0]
        args2, kwargs = mock_record.call_args_list[1]
        args3, kwargs = mock_record.call_args_list[2]
        self.assertEqual(args1[0], "CONECT 1179  746 1184 1195 1203")
        self.assertEqual(args2[0], "CONECT 1179 1211 1222")
        self.assertEqual(args3[0], "CONECT 1221  544 1017 1020 1022")
        self.assertEqual(converted_records, records)
