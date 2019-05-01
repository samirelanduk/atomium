from collections import deque
import struct
import msgpack
from datetime import date
from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.mmtf import *

class MmtfBytesToMmtfDictTests(TestCase):

    @patch("msgpack.unpackb")
    @patch("atomium.mmtf.decode_dict")
    def test_can_turn_mmtf_bytes_to_mmtf_dict(self, mock_dec, mock_unpack):
        d = mmtf_bytes_to_mmtf_dict(b"1234")
        mock_unpack.assert_called_with(b"1234")
        mock_dec.assert_called_with(mock_unpack.return_value)
        self.assertIs(d, mock_dec.return_value)



class BytesDictDecodingTests(TestCase):

    def test_can_handle_normal_dictionary_with_byte_keys(self):
        d = decode_dict({
         b"key": "value", b"key2": ["a", "b"], b"key3": 19.5
        })
        self.assertEqual(d, {
         "key": "value", "key2": ["a", "b"], "key3": 19.5
        })


    def test_can_handle_dictionary_with_byte_lists(self):
        d = decode_dict({
         b"key": "value", b"key2": [b"a", b"b"], b"key3": 19.5
        })
        self.assertEqual(d, {
         "key": "value", "key2": ["a", "b"], "key3": 19.5
        })


    @patch("atomium.mmtf.parse_binary_field")
    def test_can_handle_dictionary_with_binary_fields(self, mock_parse):
        mock_parse.side_effect = [1, 1]
        d = decode_dict({
         b"key": b"\x00\x08", b"key2": [b"a", b"b"], b"key3": "\x00123"
        })
        mock_parse.assert_any_call(b"\x00\x08")
        mock_parse.assert_any_call(b"\x00123")
        self.assertEqual(d, {
         "key": 1, "key2": ["a", "b"], "key3": 1
        })


    def test_can_handle_dictionary_with_recursive_fields(self):
        d = decode_dict({
         b"key": "value", b"key2": [{b"keya": [b"a", b"b"], b"key3": 19.5}]
        })
        self.assertEqual(d, {
         "key": "value", "key2": [{"keya": ["a", "b"], "key3": 19.5}]
        })



class BinaryFieldParsingTests(TestCase):

    def test_can_raise_exception_on_invalid_codec(self):
        with self.assertRaises(ValueError):
            parse_binary_field(struct.pack(">iii", 20, 0, 0))


    def test_can_parse_with_codec_1(self):
        b = parse_binary_field(
         struct.pack(">iii", 1, 3, 0) + struct.pack("fff", 4.5, 5.5, 6.5)
        )
        self.assertEqual(b, (4.5, 5.5, 6.5))


    def test_can_parse_with_codec_2(self):
        b = parse_binary_field(
         struct.pack(">iii", 2, 3, 0) + b"\x04\x05\x06"
        )
        self.assertEqual(b, (4, 5, 6))


    def test_can_parse_with_codec_3(self):
        b = parse_binary_field(
         struct.pack(">iii", 3, 3, 0) + struct.pack(">hhh", 4, 5, 6)
        )
        self.assertEqual(b, (4, 5, 6))


    def test_can_parse_with_codec_4(self):
        b = parse_binary_field(
         struct.pack(">iii", 4, 3, 0) + struct.pack(">iii", 4, 5, 6)
        )
        self.assertEqual(b, (4, 5, 6))


    def test_can_parse_with_codec_5(self):
        b = parse_binary_field(
         struct.pack(">iii", 5, 3, 0) + bytes([65, 0, 0, 0, 66, 0, 0, 0, 67, 0, 0, 0])
        )
        self.assertEqual(b, ["A", "B", "C"])
        b = parse_binary_field(
         struct.pack(">iii", 5, 2, 0) + bytes([65, 0, 0, 0, 68, 65, 0, 0])
        )
        self.assertEqual(b, ["A", "DA"])


    @patch("atomium.mmtf.run_length_decode")
    def test_can_parse_with_codec_6(self, mock_rldec):
        mock_rldec.return_value = [100, 0, 0, 105, 0, 0]
        b = parse_binary_field(
         struct.pack(">iii", 6, 3, 0) + struct.pack(">iii", 0, 5, 65)
        )
        mock_rldec.assert_called_with((0, 5, 65))
        self.assertEqual(b, ["d", "", "", "i", "", ""])


    @patch("atomium.mmtf.run_length_decode")
    def test_can_parse_with_codec_7(self, mock_rldec):
        mock_rldec.return_value = [100, 0, 0, 105, 0, 0]
        b = parse_binary_field(
         struct.pack(">iii", 7, 3, 0) + struct.pack(">iii", 0, 5, 65)
        )
        mock_rldec.assert_called_with((0, 5, 65))
        self.assertEqual(b, [100, 0, 0, 105, 0, 0])


    @patch("atomium.mmtf.run_length_decode")
    @patch("atomium.mmtf.delta_decode")
    def test_can_parse_with_codec_8(self, mock_deldec, mock_rldec):
        mock_rldec.return_value = [100, 0, 0, 105, 0, 0]
        mock_deldec.return_value = [1, 2, 3]
        b = parse_binary_field(
         struct.pack(">iii", 8, 3, 0) + struct.pack(">iii", 0, 5, 65)
        )
        mock_rldec.assert_called_with((0, 5, 65))
        mock_deldec.assert_called_with([100, 0, 0, 105, 0, 0])
        self.assertEqual(b, [1, 2, 3])


    @patch("atomium.mmtf.run_length_decode")
    def test_can_parse_with_codec_9(self, mock_rldec):
        mock_rldec.return_value = [1000, 0, 0, 1050, 0, 0]
        b = parse_binary_field(
         struct.pack(">iii", 9, 3, 100) + struct.pack(">iii", 0, 5, 65)
        )
        mock_rldec.assert_called_with((0, 5, 65))
        self.assertEqual(b, [10, 0, 0, 10.5, 0, 0])


    @patch("atomium.mmtf.recursive_decode")
    @patch("atomium.mmtf.delta_decode")
    def test_can_parse_with_codec_10(self, mock_deldec, mock_recdec):
        mock_recdec.return_value = [100, 0, 0, 105, 0, 0]
        mock_deldec.return_value = [1000, 2000, 3000]
        b = parse_binary_field(
         struct.pack(">iii", 10, 3, 1000) + struct.pack(">hhh", 0, 5, 65)
        )
        mock_recdec.assert_called_with((0, 5, 65))
        mock_deldec.assert_called_with([100, 0, 0, 105, 0, 0])
        self.assertEqual(b, [1, 2, 3])



class RunLengthDecodingTests(TestCase):

    def test_can_run_length_decode(self):
        self.assertEqual(
         run_length_decode([1, 10, 2, 1, 1, 4]),
         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1]
        )



class DeltaDecodingTests(TestCase):

    def test_can_delta_decode(self):
        self.assertEqual(
         delta_decode([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1]),
         [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16]
        )



class RecursiveDecodingTests(TestCase):

    def test_can_recursive_decode(self):
        self.assertEqual(
         recursive_decode([32767, 32767, 32767, 6899, 0, 2, -1, 100, -3, 5]),
         [105200, 0, 2, -1, 100, -3, 5]
        )


    def test_can_recursive_decode_other_powers(self):
        self.assertEqual(recursive_decode([
         127, 41, 34, 1, 0, -50, -128, 0, 7, 127, 0, 127, 127, 14
        ], bits=8), [168, 34, 1, 0, -50, -128, 7, 127, 268])



class MmtfDictToDataDictTests(TestCase):

    @patch("atomium.mmtf.mmtf_to_data_transfer")
    @patch("atomium.mmtf.update_models_list")
    def test_can_convert_mmtf_dict_to_data_dict(self, mock_up, mock_trans):
        m = {"A": "B", "bioAssemblyList": [{
         "name": "1", "transformList": [{
          "chainIndexList": [1, 2], "matrix": "abcdefghijklmnop"
         }, {
          "chainIndexList": [0, 2], "matrix": "ABCDEFGHIJKLMNOP"
         }]
        }], "chainIdList": ["A", "B", "C"]}
        d = mmtf_dict_to_data_dict(m)
        mock_trans.assert_any_call(m, d, "description", "code", "structureId")
        mock_trans.assert_any_call(m, d, "description", "title", "title")
        mock_trans.assert_any_call(m, d, "description", "deposition_date", "depositionDate", date=True)
        mock_trans.assert_any_call(m, d, "experiment", "technique", "experimentalMethods", first=True)
        mock_trans.assert_any_call(m, d, "quality", "resolution", "resolution", trim=3)
        mock_trans.assert_any_call(m, d, "quality", "rvalue", "rWork", trim=3)
        mock_trans.assert_any_call(m, d, "quality", "rfree", "rFree", trim=3)
        mock_trans.assert_any_call(m, d["geometry"], "crystallography", "space_group", "spaceGroup")
        mock_trans.assert_any_call(m, d["geometry"], "crystallography", "unit_cell", "unitCell", trim=3)
        mock_up.assert_called_with(m, d)
        self.assertEqual(d, {
         "description": {
          "code": None, "title": None, "deposition_date": None,
          "classification": None, "keywords": [], "authors": []
         }, "experiment": {
          "technique": None, "source_organism": None, "expression_system": None,
          "missing_residues": []
         }, "quality": {"resolution": None, "rvalue": None, "rfree": None},
         "geometry": {"assemblies": [{
          "id": 1, "software": None, "delta_energy": None,
          "buried_surface_area": None, "surface_area": None, "transformations": [{
           "chains": ["B", "C"], "matrix": ["abc", "efg", "ijk"], "vector": "dhl"
          }, {
           "chains": ["A", "C"], "matrix": ["ABC", "EFG", "IJK"], "vector": "DHL"
          }]
         }], "crystallography": {}}, "models": []
        })



class ModelListUpdatingTests(TestCase):

    @patch("atomium.mmtf.get_atoms_list")
    @patch("atomium.mmtf.get_type")
    @patch("atomium.mmtf.make_chain")
    @patch("atomium.mmtf.make_group_id")
    @patch("atomium.mmtf.make_molecule")
    @patch("atomium.mmtf.add_atom")
    def test_can_update_one_model(self, mock_ad, mock_ml, mock_gi, mock_ch, mock_tp, mock_gt):
        m = {
         "numAtoms": 18, "numGroups": 10, "numChains": 8, "numModels": 1,
         "chainsPerModel": [8], "groupsPerChain": [2, 2, 1, 1, 1, 1, 1, 1],
         "chainNameList": ["A", "B", "A", "A", "B", "B", "A", "B"],
         "groupList": [
          {"groupName": "ASN", "atomNameList": ["N", "CA"]},
          {"groupName": "XMP", "atomNameList": ["P", "O1P"]},
          {"groupName": "VAL", "atomNameList": ["N", "CA"]},
          {"groupName": "BU2", "atomNameList": ["S", "DH"]},
          {"groupName": "HOH", "atomNameList": ["O"]}
         ],
         "groupTypeList": [0, 2, 0, 2, 1, 3, 1, 3, 4, 4],
         "atomIdList": range(20), "sequenceIndexList": range(10)
        }
        mock_gt.return_value = [{"id": n} for n in range(20)]
        mock_tp.side_effect = [
         "polymer", "polymer", "non-polymer", "non-polymer", "non-polymer",
         "non-polymer", "water", "water"
        ]
        d = {"models": []}
        update_models_list(m, d)
        mock_gt.assert_called_with(m)
        for n in range(8): mock_tp.assert_any_call(m, n)
        for n in range(2): mock_ch.assert_any_call(m, n)
        mock_gi.assert_any_call(m, 0, 0)
        mock_ml.assert_any_call("polymer", "ASN", m, 0)
        mock_ml.assert_any_call("polymer", "VAL", m, 0)
        mock_ml.assert_any_call("polymer", "ASN", m, 1)
        mock_ml.assert_any_call("polymer", "VAL", m, 1)
        mock_ml.assert_any_call("non-polymer", "XMP", m, 2)
        mock_ml.assert_any_call("non-polymer", "BU2", m, 3)
        mock_ml.assert_any_call("water", "HOH", m, 6)
        mock_ml.assert_any_call("water", "HOH", m, 7)
        self.assertEqual(d, {"models": [{
         "polymer": {"A": mock_ch.return_value, "B": mock_ch.return_value},
         "non-polymer": {mock_gi.return_value: mock_ml.return_value},
         "water": {mock_gi.return_value: mock_ml.return_value}
        }]})



class AtomsListGettingTests(TestCase):

    def test_can_get_atoms(self):
        self.assertEqual(get_atoms_list({
         "xCoordList": range(3), "yCoordList": range(3, 6), "zCoordList": range(6, 9),
         "bFactorList": range(9, 12), "occupancyList": range(12, 15),
         "altLocList": ["A", "B", ""], "atomIdList":  range(15, 18),
        }), [{
         "x": 0, "y": 3, "z": 6, "bvalue": 9, "occupancy": 12, "alt_loc": "A", "id": 15
        }, {
         "x": 1, "y": 4, "z": 7, "bvalue": 10, "occupancy": 13, "alt_loc": "B", "id": 16
        }, {
         "x": 2, "y": 5, "z": 8, "bvalue": 11, "occupancy": 14, "alt_loc": None, "id": 17
        }])



class TypeGettingTests(TestCase):

    def test_can_get_polymer_type(self):
        self.assertEqual(get_type({
         "entityList": [
         {"type": "polymer", "chainIndexList": [0, 1]},
         {"type": "non-polymer", "chainIndexList": [2, 3]}
        ]}, 0), "polymer")
        self.assertEqual(get_type({
         "entityList": [
         {"type": "polymer", "chainIndexList": [0, 1]},
         {"type": "non-polymer", "chainIndexList": [2, 3]}
        ]}, 1), "polymer")


    def test_can_get_non_polymer_type(self):
        self.assertEqual(get_type({
         "entityList": [
         {"type": "polymer", "chainIndexList": [0, 1]},
         {"type": "non-polymer", "chainIndexList": [2, 3]}
        ]}, 2), "non-polymer")
        self.assertEqual(get_type({
         "entityList": [
         {"type": "polymer", "chainIndexList": [0, 1]},
         {"type": "non-polymer", "chainIndexList": [2, 3]}
        ]}, 3), "non-polymer")



class ChainCreationTests(TestCase):

    @patch("atomium.mmtf.get_chain_sequence")
    def test_can_make_chain(self, mock_seq):
        m = {"chainIdList": ["A", "B", "C", "D", "E", "F", "G", "H"]}
        self.assertEqual(make_chain(m, 2), {
         "internal_id": "C",
         "residues": {}, "sequence": mock_seq.return_value
        })
        mock_seq.assert_called_with(m, 2)



class ChainSequenceGettingTests(TestCase):

    def test_can_get_chain_sequence(self):
        self.assertEqual(get_chain_sequence({
         "entityList": [
         {"type": "polymer", "chainIndexList": [0, 1], "sequence": "TYU"},
         {"type": "polymer", "chainIndexList": [5, 6], "sequence": "APB"},
         {"type": "non-polymer", "chainIndexList": [2, 3]}
        ]}, 0), "TYU")
        self.assertEqual(get_chain_sequence({
         "entityList": [
         {"type": "polymer", "chainIndexList": [0, 1], "sequence": "TYU"},
         {"type": "polymer", "chainIndexList": [5, 6], "sequence": "APB"},
         {"type": "non-polymer", "chainIndexList": [2, 3]}
        ]}, 6), "APB")


    def test_can_get_no_sequence(self):
        self.assertEqual(get_chain_sequence({
         "entityList": [
         {"type": "polymer", "chainIndexList": [0, 1], "sequence": "TYU"},
         {"type": "polymer", "chainIndexList": [5, 6], "sequence": "APB"},
         {"type": "non-polymer", "chainIndexList": [2, 3]}
        ]}, 10), "")



class GroupIdMakingTests(TestCase):

    def test_can_make_group_id(self):
        self.assertEqual(make_group_id({
         "chainNameList": ["A", "B", "A", "A", "B", "B", "A", "B"],
         "groupIdList": [22, 23, 24], "insCodeList": ["", "A", "B"]
        }, 4, 1), "B.23A")



class MoleculeMakingTests(TestCase):

    def test_can_return_residue(self):
        self.assertEqual(make_molecule("polymer", "VAL", {}, 1), {
         "name": "VAL", "atoms": {}
        })


    def test_can_return_ligand(self):
        self.assertEqual(make_molecule("non-polymer", "XMP", {
         "chainNameList": ["A", "B", "A", "A", "B", "B", "A", "B"],
         "chainIdList": ["A", "B", "C", "D", "E", "F", "G", "H"]
        }, 3), {
         "name": "XMP", "atoms": {}, "internal_id": "D", "polymer": "A"
        })



class AtomAddingTests(TestCase):

    def test_can_add_atoms(self):
        atoms = {"atoms": {1: 2}}
        add_atom(atoms, [{"id": 1}, {"id": 2}, {"id": 3}, {"id": 4}], {
         "atomNameList": "AB", "elementList": "CD", "formalChargeList": [0, 3]
        }, 1, 2)
        self.assertEqual(atoms, {"atoms": {1: 2, 4: {
         "name": "B", "element": "D", "charge": 3, "anisotropy": [0] * 6
        }}})



class MmtfDictTransferTests(TestCase):

    def setUp(self):
        self.d = {"A": {1: None, 2: None, 3: None}, "B": {4: None, 5: None}}
        self.m = {"M": 10.127, "C": 100, "D": "2018-09-17", "L": [8, 9]}


    def test_can_do_nothing(self):
        mmtf_to_data_transfer(self.m, self.d, "A", 1, "X", 10)
        self.assertEqual(self.d["A"][1], None)
        mmtf_to_data_transfer(self.m, self.d, "A", 1, "M", 100)
        self.assertEqual(self.d["A"][1], None)


    def test_can_transfer_from_mmtf_to_data_dict(self):
        mmtf_to_data_transfer(self.m, self.d, "A", 1, "M")
        self.assertEqual(self.d["A"][1], 10.127)
        mmtf_to_data_transfer(self.m, self.d, "B", 5, "C")
        self.assertEqual(self.d["B"][5], 100)


    def test_can_transfer_from_mmtf_to_data_dict_date(self):
        mmtf_to_data_transfer(self.m, self.d, "B", 5, "D", date=True)
        self.assertEqual(self.d["B"][5], date(2018, 9, 17))


    def test_can_transfer_from_mmtf_to_data_dict_first(self):
        mmtf_to_data_transfer(self.m, self.d, "B", 5, "L", first=True)
        self.assertEqual(self.d["B"][5], 8)


    def test_can_transfer_from_mmtf_to_data_dict_round(self):
        mmtf_to_data_transfer(self.m, self.d, "B", 5, "M", trim=1)
        self.assertEqual(self.d["B"][5], 10.1)
        mmtf_to_data_transfer(self.m, self.d, "B", 5, "M", trim=2)
        self.assertEqual(self.d["B"][5], 10.13)
        self.m["M"] = [10.13122334, 1.119973]
        mmtf_to_data_transfer(self.m, self.d, "B", 5, "M", trim=2)
        self.assertEqual(self.d["B"][5], [10.13, 1.12])



class StructureToMmtfStringTests(TestCase):

    @patch("atomium.mmtf.get_structures")
    @patch("atomium.mmtf.get_entity_list")
    @patch("atomium.mmtf.get_chain_ids_and_names")
    @patch("atomium.mmtf.get_groups_per_chain")
    @patch("atomium.mmtf.get_groups")
    def test_can_convert_structure_to_mmtf_string(self, mock_gr, mock_gpc, mock_ch, mock_en, mock_st):
        structure = Mock()
        mock_st.return_value = ("AB", "CD", [
         Mock(chain="A"), Mock(chain="A"), Mock(chain="B"), Mock(chain="B")
        ], [
         [1, 2, 3, 4, 5, 6, 7], [8, 9, 10, 11, 12, 13, 14]
        ], [])
        mock_en.return_value = "GH"
        mock_ch.return_value = ([10, 20], [30, 40])
        mock_gpc.return_value = [8, 4]
        mock_gr.return_value = (100, 200, 300, 400)
        s = structure_to_mmtf_string(structure)
        self.assertEqual(s, msgpack.packb({
         "numModels": 1, "numChains": 6, "chainsPerModel": [6],
         "xCoordList": [1, 8], "yCoordList": [2, 9], "zCoordList": [3, 10], "altLocList": [4, 11],
         "bFactorList": [5, 12], "atomIdList": [6, 13], "occupancyList": [7, 14],
         "entityList": "GH", "chainIdList": [10, 20], "insCodeList": 400,
         "chainNameList": [30, 40], "groupsPerChain": [8, 4],
         "groupList": 300, "groupIdList": 200, "groupTypeList": 100
        }))



class StructureGettingTests(TestCase):

    @patch("atomium.mmtf.get_structure_from_atom")
    @patch("atomium.mmtf.create_entities")
    def test_can_get_structures(self, mock_cr, mock_get):
        s = Mock()
        s.atoms.return_value = [
         Mock(id=2, x=10, y=11, z=13, bvalue=14), Mock(id=1, x=20, y=21, z=23, bvalue=24)
        ]
        chains, ligands, waters, atom_properties, entities  = get_structures(s)
        mock_get.assert_any_call(s.atoms.return_value[0], set(), set(), set())
        mock_get.assert_any_call(s.atoms.return_value[1], set(), set(), set())
        mock_cr.assert_called_with([], [], [])
        self.assertEqual(chains, [])
        self.assertEqual(ligands, [])
        self.assertEqual(waters, [])
        self.assertEqual(atom_properties, [[20, 21, 23, "", 24, 1, 1], [10, 11, 13, "", 14, 2, 1]])
        self.assertEqual(entities, mock_cr.return_value)



class EntityListGenerationTests(TestCase):

    def test_can_get_entity_list(self):
        chains = [Mock(sequence="E"), Mock(sequence="E")]
        ligands = [Mock(_name="A"), Mock(_name="B")]
        waters = [Mock(), Mock()]
        entities = [Mock(Chain, sequence="E"), Mock(Ligand, water=False, _name="A"),
         Mock(Ligand, water=False, _name="B"), Mock(Ligand)]
        e = get_entity_list(entities, chains, ligands, waters)
        self.assertEqual(e, [{
         "type": "polymer", "sequence": "E", "chainIndexList": [0, 1]
        }, {
         "type": "non-polymer", "chainIndexList": [2]
        }, {
         "type": "non-polymer", "chainIndexList": [3]
        }, {
         "type": "water", "chainIndexList": [4, 5]
        }])



class ChainIdAndNameGettingTests(TestCase):

    def test_can_get_chain_ids_and_names(self):
        chains = [Mock(_internal_id="A", id="A"), Mock(_internal_id="B", id="B")]
        ligands = [Mock(_internal_id="C", chain=Mock(id="A")), Mock(_internal_id="D", chain=Mock(id="B"))]
        waters = [Mock(_internal_id="E", chain=Mock(id="A")), Mock(_internal_id="E", id="B")]
        ids, names = get_chain_ids_and_names(chains, ligands, waters)
        self.assertEqual(ids, ["A", "B", "C", "D", "E"])
        self.assertEqual(names, ["A", "B", "A", "B", "A"])



class GroupPerChainCountTests(TestCase):

    def test_can_get_chain_counts(self):
        chains = [Mock(), Mock()]
        chains[0].residues.return_value = "AB"
        chains[1].residues.return_value = "CDE"
        ligands = [Mock(), Mock()]
        waters = [Mock(chain="A", _internal_id=1), Mock(chain="A", _internal_id=1),
         Mock(chain="B", _internal_id=2), Mock(chain="B", _internal_id=2), Mock(chain="B", _internal_id=2)]
        counts = get_groups_per_chain(chains, ligands, waters)
        self.assertEqual(counts, [2, 3, 1, 1, 2, 3])



class GroupCreationTests(TestCase):

    @patch("atomium.mmtf.add_het_to_groups")
    def test_can_get_groups(self, mock_add):
        chains = [Mock(), Mock()]
        chains[0].residues.return_value = "AB"
        chains[1].residues.return_value = "CDE"
        ligands = [Mock(), Mock()]
        waters = [Mock(), Mock()]
        a, b, c, d = get_groups(chains, ligands, waters)
        mock_add.assert_any_call("A", [], [], [], [])
        mock_add.assert_any_call("B", [], [], [], [])
        mock_add.assert_any_call("C", [], [], [], [])
        mock_add.assert_any_call("D", [], [], [], [])
        mock_add.assert_any_call("E", [], [], [], [])
        mock_add.assert_any_call(ligands[0], [], [], [], [])
        mock_add.assert_any_call(ligands[1], [], [], [], [])
        mock_add.assert_any_call(waters[0], [], [], [], [])
        mock_add.assert_any_call(waters[1], [], [], [], [])
        self.assertEqual((a, b, c, d), ([], [], [], []))



class HetToGroupAdditionTests(TestCase):

    @patch("atomium.mmtf.split_residue_id")
    def test_can_add_new_het_to_group(self, mock_split):
        a, b, c, d = [], [], [], []
        het = Mock(_name="HET")
        mock_split.return_value = "AB"
        het.atoms.return_value = [Mock(id=2, _name="AA", element="A", charge=-1), Mock(id=1, _name="BB", element="B", charge=-2)]
        add_het_to_groups(het, a, b, c, d)
        self.assertEqual(a, [0])
        self.assertEqual(b, ["A"])
        self.assertEqual(c, [{
         "groupName": "HET", "atomNameList": ["BB", "AA"],
         "elementList": ["B", "A"], "formalChargeList": [-2, -1]
        }])
        self.assertEqual(d, ["B"])


    @patch("atomium.mmtf.split_residue_id")
    def test_can_add_existing_het_to_group(self, mock_split):
        a, b, c, d = [0], ["A"], [{
         "groupName": "HET", "atomNameList": ["BB", "AA"],
         "elementList": ["B", "A"], "formalChargeList": [-2, -1]
        }], ["B"]
        het = Mock(_name="HET")
        mock_split.return_value = "A?"
        het.atoms.return_value = [Mock(id=2, _name="AA", element="A", charge=-1), Mock(id=1, _name="BB", element="B", charge=-2)]
        add_het_to_groups(het, a, b, c, d)
        self.assertEqual(a, [0, 0])
        self.assertEqual(b, ["A", "A"])
        self.assertEqual(c, [{
         "groupName": "HET", "atomNameList": ["BB", "AA"],
         "elementList": ["B", "A"], "formalChargeList": [-2, -1]
        }])
        self.assertEqual(d, ["B", ""])
