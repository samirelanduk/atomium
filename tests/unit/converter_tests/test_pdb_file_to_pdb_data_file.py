from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.converters.pdbfile2pdbdatafile import *
from atomium.files.pdbfile import PdbFile, PdbRecord
from atomium.files.pdbdatafile import PdbDataFile

class PdbFileToPdbDataFileTests(TestCase):

    @patch("atomium.converters.pdbfile2pdbdatafile.extract_structure")
    def test_can_create_data_file_from_pdb_file(self, mock_extract):
        pdb_file = Mock(PdbFile)
        data_file = pdb_file_to_pdb_data_file(pdb_file)
        self.assertIsInstance(data_file, PdbDataFile)
        mock_extract.assert_called_with(pdb_file, data_file)



class StructureExtractionTests(TestCase):

    @patch("atomium.converters.pdbfile2pdbdatafile.atom_record_to_dict")
    @patch("atomium.converters.pdbfile2pdbdatafile.conect_records_to_list")
    def test_can_extract_structure_from_pdb_file(self, mock_conect, mock_atom):
        pdb_file = Mock(PdbFile)
        data_file = Mock(PdbDataFile)
        atom_recs = (Mock(), Mock())
        hetatom_recs = (Mock(), Mock())
        conect_recs = (Mock(), Mock())
        model_recs = (Mock(), Mock())
        mock_conect.return_value = "connections"
        mock_atom.side_effect = ["a1", "a2", "h1", "h2"]
        pdb_file.records.side_effect = [model_recs, atom_recs, hetatom_recs, conect_recs]
        extract_structure(pdb_file, data_file)
        pdb_file.records.assert_any_call(name="model")
        pdb_file.records.assert_any_call(name="atom")
        pdb_file.records.assert_any_call(name="hetatm")
        pdb_file.records.assert_any_call(name="conect")
        mock_conect.assert_called_with(conect_recs)
        mock_atom.assert_any_call(atom_recs[0], model_recs)
        mock_atom.assert_any_call(atom_recs[1], model_recs)
        mock_atom.assert_any_call(hetatom_recs[0], model_recs)
        mock_atom.assert_any_call(hetatom_recs[1], model_recs)
        self.assertEqual(data_file.atoms, ["a1", "a2"])
        self.assertEqual(data_file.heteroatoms, ["h1", "h2"])
        self.assertEqual(data_file.connections, "connections")



class AtomRecordToAtomDictTests(TestCase):

    def setUp(self):
        self.record = PdbRecord("ATOM    107  N   GLY A  13A     12.681  " +
         "37.302 -25.211 1.000 15.56           N-2")
        self.record.number = lambda: 121
        self.model_records = [Mock(), Mock(), Mock()]
        self.model_records[0].number.return_value = 101
        self.model_records[1].number.return_value = 119
        self.model_records[2].number.return_value = 213


    def test_can_convert_atom_record_no_models(self):
        self.assertEqual(atom_record_to_dict(self.record, []), {
         "atom_id": 107, "atom_name": "N", "alt_loc": None,
         "residue_name": "GLY",
         "chain_id": "A", "residue_id": 13, "insert_code": "A",
         "x": 12.681, "y": 37.302, "z": -25.211,
         "occupancy": 1.0, "temperature_factor": 15.56,
         "element": "N", "charge": -2, "model": 1
        })


    def test_can_convert_atom_record_with_models(self):
        self.assertEqual(atom_record_to_dict(self.record, self.model_records), {
         "atom_id": 107, "atom_name": "N", "alt_loc": None,
         "residue_name": "GLY",
         "chain_id": "A", "residue_id": 13, "insert_code": "A",
         "x": 12.681, "y": 37.302, "z": -25.211,
         "occupancy": 1.0, "temperature_factor": 15.56,
         "element": "N", "charge": -2, "model": 2
        })


    def test_can_convert_opposite_charge_order(self):
        self.record = PdbRecord("ATOM    107  N   GLY A  13A     12.681  " +
         "37.302 -25.211 1.000 15.56           N2-")
        d = atom_record_to_dict(self.record, [])
        self.assertEqual(d["charge"], -2)



class ConectRecordsToConectListTests(TestCase):

    @patch("atomium.converters.pdbfile2pdbdatafile.merge_records")
    def test_can_convert_conect_records_to_list(self, mock_merge):
        mock_merge.side_effect = [
         "746 1184 1195 1203 1211 1222", "544 1017 1020 1022"
        ]
        record1 = PdbRecord("CONECT 1179  746 1184 1195 1203")
        record2 = PdbRecord("CONECT 1179 1211 1222")
        record3 = PdbRecord("CONECT 1221  544 1017 1020 1022")
        connections = conect_records_to_list((record1, record2, record3))
        self.assertEqual(connections, [{
         "atom": 1179, "bond_to": [746, 1184, 1195, 1203, 1211, 1222]
        }, {
         "atom": 1221, "bond_to": [544, 1017, 1020, 1022]
        }])



class RecordMergingTests(TestCase):

    def setUp(self):
        self.records = [PdbRecord(l) for l in [
         "0123456789",
         "abcdefghij",
         "0123456789"
        ]]
        self.punc_records = [PdbRecord(l) for l in [
         "0123, 456789",
         "abcd  efghij",
         "0123; 456789"
        ]]


    def test_can_merge_records(self):
        self.assertEqual(
         merge_records(self.records, 5),
         "56789 fghij 56789"
        )
        self.assertEqual(
         merge_records(self.records, 8),
         "89 ij 89"
        )


    def test_can_vary_join(self):
        self.assertEqual(
         merge_records(self.records, 5, join=""),
         "56789fghij56789"
        )
        self.assertEqual(
         merge_records(self.records, 8, join="."),
         "89.ij.89"
        )


    def test_can_condense(self):
        self.assertEqual(
         merge_records(self.punc_records, 2),
         "23,456789 cd efghij 23;456789"
        )


    def test_can_ignore_consensors(self):
        self.assertEqual(
         merge_records(self.punc_records, 2, dont_condense=","),
         "23, 456789 cd efghij 23;456789"
        )
        self.assertEqual(
         merge_records(self.punc_records, 2, dont_condense=";"),
         "23,456789 cd efghij 23; 456789"
        )
        self.assertEqual(
         merge_records(self.punc_records, 2, dont_condense=";,"),
         "23, 456789 cd efghij 23; 456789"
        )
        self.assertEqual(
         merge_records(self.punc_records, 2, dont_condense=";, "),
         "23, 456789 cd  efghij 23; 456789"
        )
