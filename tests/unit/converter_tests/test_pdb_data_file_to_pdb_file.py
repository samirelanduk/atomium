from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.converters.pdbdatafile2pdbfile import *
from atomium.files.pdbfile import PdbFile
from atomium.files.pdbdatafile import PdbDataFile

class PdbDataFileToPdbFileTests(TestCase):

    @patch("atomium.converters.pdbdatafile2pdbfile.pack_structure")
    def test_can_create_pdb_file_from_data_file(self, mock_pack):
        data_file = Mock(PdbDataFile)
        pdb_file = pdb_data_file_to_pdb_file(data_file)
        self.assertIsInstance(pdb_file, PdbFile)
        mock_pack.assert_called_with(data_file, pdb_file)



class StructurePackingTests(TestCase):

    def setUp(self):
        self.atoms = [Mock(), Mock(), Mock(), Mock()]
        self.pdb_file = Mock(PdbFile)
        self.pdb_file._records = []
        self.data_file = Mock(PdbFile)
        self.data_file.connections = "ccc"
        self.data_file.atoms = [{"model": 1, "a": 1}, {"model": 1, "a": 2}]
        self.data_file.heteroatoms = [{"model": 1, "a": 3}, {"model": 1, "a": 4}]


    @patch("atomium.converters.pdbdatafile2pdbfile.atom_dict_to_record")
    @patch("atomium.converters.pdbdatafile2pdbfile.conections_list_to_records")
    def test_can_pack_structure_one_model(self, mock_connections, mock_atom):
        mock_connections.return_value = ["c1", "c2"]
        mock_atom.side_effect = self.atoms
        pack_structure(self.data_file, self.pdb_file)
        mock_atom.assert_any_call({"model": 1, "a": 1})
        mock_atom.assert_any_call({"model": 1, "a": 2})
        mock_atom.assert_any_call({"model": 1, "a": 3}, hetero=True)
        mock_atom.assert_any_call({"model": 1, "a": 4}, hetero=True)
        mock_connections.assert_called_with("ccc")
        self.assertEqual(self.pdb_file._records, self.atoms + ["c1", "c2"])


    @patch("atomium.converters.pdbdatafile2pdbfile.atom_dict_to_record")
    @patch("atomium.converters.pdbdatafile2pdbfile.conections_list_to_records")
    @patch("atomium.converters.pdbdatafile2pdbfile.PdbRecord")
    def test_can_pack_structure_two_models(self, mock_record, mock_connections, mock_atom):
        self.data_file.atoms[1]["model"] = 2
        self.data_file.heteroatoms[1]["model"] = 2
        mock_record.side_effect = ["model", "end", "model", "end"]
        mock_connections.return_value = ["c1", "c2"]
        mock_atom.side_effect = self.atoms
        pack_structure(self.data_file, self.pdb_file)
        mock_atom.assert_any_call({"model": 1, "a": 1})
        mock_atom.assert_any_call({"model": 2, "a": 2})
        mock_atom.assert_any_call({"model": 1, "a": 3}, hetero=True)
        mock_atom.assert_any_call({"model": 2, "a": 4}, hetero=True)
        mock_connections.assert_called_with("ccc")
        mock_record.assert_any_call("MODEL        1")
        mock_record.assert_any_call("MODEL        2")
        mock_record.assert_any_call("ENDMDL")
        self.assertEqual(self.pdb_file._records, [
         "model", self.atoms[0], self.atoms[2], "end",
         "model", self.atoms[1], self.atoms[3], "end", "c1", "c2"
        ])



class AtomDictToAtomRecordTests(TestCase):

    def setUp(self):
        self.atom_dict = {
         "atom_id": 107, "atom_name": "N", "alt_loc": None,
         "residue_name": "GLY",
         "chain_id": "A", "residue_id": 13, "insert_code": "A",
         "x": 12.681, "y": 37.302, "z": -25.211,
         "occupancy": 1.0, "temperature_factor": 15.56,
         "element": "N", "charge": -2, "model": 5
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
        mock_record.assert_called_with("ATOM" + " " * 51 + "1.000" + " " * 20)


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
