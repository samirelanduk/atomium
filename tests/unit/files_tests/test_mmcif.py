from copy import deepcopy
from datetime import date
from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.files.mmcif import *

class MmcifStringToMmcifDictTests(TestCase):

    @patch("atomium.files.mmcif.consolidate_strings")
    @patch("atomium.files.mmcif.mmcif_lines_to_mmcif_blocks")
    @patch("atomium.files.mmcif.loop_block_to_list")
    @patch("atomium.files.mmcif.category_block_to_dict")
    @patch("atomium.files.mmcif.strip_quotes")
    def test_can_turn_mmcif_string_to_mmcif_dict(self, mock_strp, mock_cat, mock_loop, mock_block, mock_con):
        mock_con.return_value = "CONSOLIDATED"
        mock_block.return_value = [
         {"category": "A", "lines": ["1", "2"]},
         {"category": "B", "lines": ["loop_", "3"]},
         {"category": "C", "lines": ["loop_", "4"]},
         {"category": "D", "lines": ["5", "6"]},
        ]
        mock_cat.side_effect = ["CAT1", "CAT2"]
        mock_loop.side_effect = ["LOOP1", "LOOP2"]
        d = mmcif_string_to_mmcif_dict("1\n\n2\n3\n")
        mock_con.assert_called_with(["1", "2", "3"])
        mock_block.assert_called_with("CONSOLIDATED")
        mock_loop.assert_any_call({"category": "B", "lines": ["loop_", "3"]})
        mock_loop.assert_any_call({"category": "C", "lines": ["loop_", "4"]})
        mock_cat.assert_any_call({"category": "A", "lines": ["1", "2"]})
        mock_cat.assert_any_call({"category": "D", "lines": ["5", "6"]})
        mock_strp.assert_called_with({
         "A": "CAT1", "B": "LOOP1", "C": "LOOP2", "D": "CAT2"
        })
        self.assertEqual(d, {
         "A": "CAT1", "B": "LOOP1", "C": "LOOP2", "D": "CAT2"
        })



class StringConsolidationTests(TestCase):

    def test_consolidation_can_do_nothing(self):
        lines = consolidate_strings(["line1", "line2", "line3"])
        self.assertEqual(lines, ["line1", "line2", "line3"])


    def test_can_consolidate_string(self):
        lines = consolidate_strings(["line1", "line2", ";STRING", ";", "line3"])
        self.assertEqual(lines, ["line1", "line2 \"STRING\"", "line3"])


    def test_can_consolidate_multi_strings(self):
        lines = consolidate_strings([
         "line1", "line2", "; STRING 1", "CONT CONT", ";", "line3", ";STRING2", ";", "line4"
        ])
        self.assertEqual(lines, [
         "line1", "line2 \"STRING 1 CONT CONT\"", "line3 \"STRING2\"", "line4"
        ])



class MmcifLinesToMmcifBlockTests(TestCase):

    def test_can_make_category_blocks(self):
        blocks = mmcif_lines_to_mmcif_blocks([
         "data_1", "_AAA.x 100", "_BBB.x 200", "_BBB.y 200", "_CCC.z 300"
        ])
        self.assertEqual(blocks, [
         {"category": "AAA", "lines": ["_AAA.x 100"]},
         {"category": "BBB", "lines": ["_BBB.x 200", "_BBB.y 200"]},
         {"category": "CCC", "lines": ["_CCC.z 300"]}
        ])


    def test_can_make_loop_blocks(self):
        blocks = mmcif_lines_to_mmcif_blocks([
         "data_1", "loop_", "_AAA.x", "100", "loop_", "_BBB.x", "_BBB.y", "11", "12", "13"
        ])
        self.assertEqual(blocks, [
         {"category": "AAA", "lines": ["loop_", "_AAA.x", "100"]},
         {"category": "BBB", "lines": ["loop_", "_BBB.x", "_BBB.y", "11", "12", "13"]},
        ])


    def test_can_make_mixed_blocks(self):
        blocks = mmcif_lines_to_mmcif_blocks([
         "data_1", "_AAA.x 100", "loop_", "_111.x", "100", "_BBB.x 200",
         "_BBB.y 200", "loop_", "_222.x", "_222.y", "11", "12", "13", "_CCC.z 300"
        ])
        self.assertEqual(blocks, [
         {"category": "AAA", "lines": ["_AAA.x 100"]},
          {"category": "111", "lines": ["loop_", "_111.x", "100"]},
         {"category": "BBB", "lines": ["_BBB.x 200", "_BBB.y 200"]},
          {"category": "222", "lines": ["loop_", "_222.x", "_222.y", "11", "12", "13"]},
         {"category": "CCC", "lines": ["_CCC.z 300"]}
        ])



class CategoryBlockToDictTests(TestCase):

    def test_can_convert_category_block_to_dict(self):
        d = category_block_to_dict({"category": "AA", "lines": [
         "_AA.one XXX", "_AA.two 'YYY ZZZ'", "_AA.three 1.2.3"
        ]})
        self.assertEqual(d, [{
         "one": "XXX", "two": "YYY ZZZ", "three": "1.2.3"
        }])


    def test_can_move_string_in_category_block(self):
        d = category_block_to_dict({"category": "AA", "lines": [
         "_AA.one XXX", "_AA.two", "'YYY ZZZ'", "_AA.three 1.2.3"
        ]})
        self.assertEqual(d, [{
         "one": "XXX", "two": "YYY ZZZ", "three": "1.2.3"
        }])


class LoopBlockToListTests(TestCase):

    @patch("atomium.files.mmcif.split_values")
    def test_can_convert_loop_block_to_list(self, mock_split):
        mock_split.side_effect = [["1", "2", "3"], ["Wu, N.", "5", "6"]]
        l = loop_block_to_list({"category": "AA", "lines": [
         "loop_", "_AA.one", "_AA.two", "_AA.three", "1 2 3", "'Wu, N.'    5  6"
        ]})
        self.assertEqual(l, [{
         "one": "1", "two": "2", "three": "3"
        }, {
         "one": "Wu, N.", "two": "5", "three": "6"
        }])
        mock_split.assert_any_call("1 2 3")
        mock_split.assert_any_call("'Wu, N.'    5  6")


    @patch("atomium.files.mmcif.split_values")
    def test_can_convert_loop_block_to_list_broken_lines(self, mock_split):
        mock_split.side_effect = [["1", "2", "3"], ["Wu, N."], ["5", "6"]]
        l = loop_block_to_list({"category": "AA", "lines": [
         "loop_", "_AA.one", "_AA.two", "_AA.three", "1 2 3", "'Wu, N.'", "    5  6"
        ]})
        self.assertEqual(l, [{
         "one": "1", "two": "2", "three": "3"
        }, {
         "one": "Wu, N.", "two": "5", "three": "6"
        }])
        mock_split.assert_any_call("1 2 3")
        mock_split.assert_any_call("'Wu, N.'")
        mock_split.assert_any_call("    5  6")



class ValueSplittingTests(TestCase):

    def test_can_split_basic_line(self):
        self.assertEqual(split_values("1 2 3"), ["1", "2", "3"])


    def test_can_split_string_line(self):
        self.assertEqual(split_values("1 '2 5' 3"), ["1", "2 5", "3"])
        self.assertEqual(split_values("1 \"2 5\" 3"), ["1", "2 5", "3"])


    def test_can_split_string_line_with_quote(self):
        self.assertEqual(split_values("1 '2 '5' 3"), ["1", "2 '5", "3"])



class QuoteStrippingTests(TestCase):

    def test_can_strip_quotes(self):
        d = {
         "A": {
          "1": "20", "2": '"34"'
         }, "B": [{
          "3": "30", "4": '"34"'
         }, {
          "3": "50", "4": '"34 dfe"'
         }]
        }
        strip_quotes(d)
        self.assertEqual(d, {
         "A": {
          "1": "20", "2": '34'
         }, "B": [{
          "3": "30", "4": '34'
         }, {
          "3": "50", "4": '34 dfe'
         }]
        })



class MmcifDictToDataDictTests(TestCase):

    @patch("atomium.files.mmcif.update_description_dict")
    @patch("atomium.files.mmcif.update_experiment_dict")
    @patch("atomium.files.mmcif.update_quality_dict")
    @patch("atomium.files.mmcif.update_geometry_dict")
    @patch("atomium.files.mmcif.update_models_list")
    def test_can_convert_mmcif_dict_to_data_dict(self, mock_mod, mock_geom, mock_qual, mock_exp, mock_desc):
        self.assertEqual(mmcif_dict_to_data_dict({"A": "B"}), DATA_DICT)
        mock_desc.assert_called_with({"A": "B"}, DATA_DICT)
        mock_exp.assert_called_with({"A": "B"}, DATA_DICT)
        mock_qual.assert_called_with({"A": "B"}, DATA_DICT)
        mock_geom.assert_called_with({"A": "B"}, DATA_DICT)
        mock_mod.assert_called_with({"A": "B"}, DATA_DICT)



class DescriptionDictUpdatingTests(TestCase):

    def setUp(self):
        self.mmcif_dict = {
         "entry": [{"id": "1LOL"}],
         "struct": [{"title": "TTT"}],
         "pdbx_database_status": [{"recvd_initial_deposition_date": "1990-09-28"}],
         "struct_keywords": [{"pdbx_keywords": "CLS", "text": "One,Two, Three"}],
         "audit_author": [{"name": "A"}, {"name": "B"}]
        }


    def test_can_create_description_dict(self):
        d = deepcopy(DATA_DICT)
        update_description_dict(self.mmcif_dict, d)
        self.assertEqual(d["description"], {
         "code": "1LOL",
         "title": "TTT",
         "deposition_date": date(1990, 9, 28),
         "classification": "CLS",
         "keywords": ["One", "Two", "Three"],
         "authors": ["A", "B"]
        })


    def test_can_handle_missing_keys(self):
        for key in self.mmcif_dict:
            self.mmcif_dict[key] = {}
        self.mmcif_dict["audit_author"] = [{}, {}]
        d = deepcopy(DATA_DICT)
        update_description_dict(self.mmcif_dict, d)
        self.assertEqual(d["description"], DATA_DICT["description"])
        update_description_dict({}, d)
        self.assertEqual(d["description"], DATA_DICT["description"])



class ExperimentDictUpdatingTests(TestCase):

    def setUp(self):
        self.mmcif_dict = {
         "exptl": [{
          "method": "NMR"
         }], "entity_src_gen": [{
          "pdbx_gene_src_scientific_name": "HS",
          "pdbx_host_org_scientific_name": "EC"
         }],
        }


    def test_can_create_quality_dict(self):
        d = deepcopy(DATA_DICT)
        update_experiment_dict(self.mmcif_dict, d)
        self.assertEqual(d["experiment"], {
         "technique": "NMR", "source_organism": "HS", "expression_system": "EC"
        })


    def test_can_handle_missing_keys(self):
        for key in self.mmcif_dict:
            self.mmcif_dict[key] = {}
        d = deepcopy(DATA_DICT)
        update_experiment_dict(self.mmcif_dict, d)
        self.assertEqual(d["experiment"], DATA_DICT["experiment"])
        update_experiment_dict({}, d)
        self.assertEqual(d["experiment"], DATA_DICT["experiment"])



class QualityDictUpdatingTests(TestCase):

    def setUp(self):
        self.mmcif_dict = {
         "reflns": [{
          "d_resolution_high": "4.5"
         }], "refine": [{
          "ls_R_factor_R_work": "0.3",
          "ls_R_factor_R_free": "0.2"
         }],
        }


    def test_can_create_quality_dict(self):
        d = deepcopy(DATA_DICT)
        update_quality_dict(self.mmcif_dict, d)
        self.assertEqual(d["quality"], {
         "resolution": 4.5, "rvalue": 0.3, "rfree": 0.2
        })


    def test_can_handle_missing_keys(self):
        for key in self.mmcif_dict:
            self.mmcif_dict[key] = {}
        d = deepcopy(DATA_DICT)
        update_quality_dict(self.mmcif_dict, d)
        self.assertEqual(d["quality"], DATA_DICT["quality"])
        update_quality_dict({}, d)
        self.assertEqual(d["quality"], DATA_DICT["quality"])



class GeometryDictUpdatingTests(TestCase):

    @patch("atomium.files.mmcif.assign_software_to_assembly")
    @patch("atomium.files.mmcif.assign_metrics_to_assembly")
    @patch("atomium.files.mmcif.assign_transformations_to_assembly")
    def test_can_update_geometry_dict(self, mock_tran, mock_met, mock_soft):
        mmcif_dict = {"pdbx_struct_assembly": [{"id": "1"}, {"id": "2"}]}
        d = deepcopy(DATA_DICT)
        update_geometry_dict(mmcif_dict, d)
        self.assertEqual(d["geometry"]["assemblies"], [{
         "id": 1, "software": None, "transformations": [],
         "delta_energy": None, "surface_area": None, "buried_surface_area": None,
        }, {
         "id": 2, "software": None, "transformations": [],
         "delta_energy": None, "surface_area": None, "buried_surface_area": None,
        }])
        mock_soft.assert_any_call(d["geometry"]["assemblies"][0], mmcif_dict)
        mock_soft.assert_any_call(d["geometry"]["assemblies"][1], mmcif_dict)
        mock_met.assert_any_call(d["geometry"]["assemblies"][0], mmcif_dict)
        mock_met.assert_any_call(d["geometry"]["assemblies"][1], mmcif_dict)
        mock_tran.assert_any_call(d["geometry"]["assemblies"][0], mmcif_dict)
        mock_tran.assert_any_call(d["geometry"]["assemblies"][1], mmcif_dict)



class ModelsListUpdatingTests(TestCase):

    @patch("atomium.files.mmcif.atom_line_to_atom_dict")
    @patch("atomium.files.mmcif.generate_higher_structures")
    def test_can_create_models(self, mock_gen, mock_atom):
        mmcif_dict = {
         "atom_site": [
          {"pdbx_PDB_model_num": 1, "id": 1},
          {"pdbx_PDB_model_num": 1, "id": 2},
          {"pdbx_PDB_model_num": 1, "id": 3},
          {"pdbx_PDB_model_num": 2, "id": 4},
          {"pdbx_PDB_model_num": 2, "id": 5},
          {"pdbx_PDB_model_num": 2, "id": 6}
         ],
         "pdbx_poly_seq_scheme": [
          {"pdb_strand_id": "A", "mon_id": "VAL"}, {"pdb_strand_id": "A", "mon_id": "TYR"},
          {"pdb_strand_id": "B", "mon_id": "TRP"}, {"pdb_strand_id": "B", "mon_id": "MET"},
         ],
         "atom_site_anisotrop": [
          {"id": 3, "U[1][1]": 123, "U[1][2]": 234, "U[1][3]": 345,
           "U[2][2]": 456, "U[2][3]": 567, "U[3][3]": 678}
         ]
        }
        mock_gen.side_effect = lambda models: [
         [m["chains"].append({"id": n, "full_sequence": []}) for n in "AB"] for m in models
        ]
        mock_atom.side_effect = [{"id": n, "anisotropy": []} for n in range(1, 7)]
        d = deepcopy(DATA_DICT)
        update_models_list(mmcif_dict, d)
        self.assertEqual(d["models"], [{
         "atoms": [{"id": n, "anisotropy": [0, 0, 0, 0, 0, 0]} for n in range(1, 3)] +
          [{"id": 3, "anisotropy": [0.0123, 0.0234, 0.0345, 0.0456, 0.0567, 0.0678]}],
         "residues": [], "ligands": [], "connections": [],
         "chains": [{"id": "A", "full_sequence": ["VAL", "TYR"]}, {"id": "B", "full_sequence": ["TRP", "MET"]}]
        }, {
         "atoms": [{"id": n, "anisotropy": [0, 0, 0, 0, 0, 0]} for n in range(4, 7)],
         "residues": [], "ligands": [], "connections": [],
         "chains": [{"id": "A", "full_sequence": ["VAL", "TYR"]}, {"id": "B", "full_sequence": ["TRP", "MET"]}]
        }])
        mock_gen.assert_called_with(d["models"])



class AssemblySoftwareAssigningTests(TestCase):

    def test_can_handle_no_software(self):
        assembly = {"id": 2, "software": None}
        mmcif_dict = {"pdbx_struct_assembly": [{
         "id": "1", "method_details": "?"
        }, {
         "id": "2", "method_details": "?"
        }]}
        assign_software_to_assembly(assembly, mmcif_dict)
        self.assertEqual(assembly["software"], None)
        assign_software_to_assembly(assembly, {})
        self.assertEqual(assembly["software"], None)


    def test_can_update_software(self):
        assembly = {"id": 2, "software": None}
        mmcif_dict = {"pdbx_struct_assembly": [{
         "id": "1", "method_details": "AAA"
        }, {
         "id": "2", "method_details": "BBB"
        }]}
        assign_software_to_assembly(assembly, mmcif_dict)
        self.assertEqual(assembly["software"], "BBB")



class AssemblyMetricsAssigningTests(TestCase):

    def test_can_handle_no_metrics(self):
        assembly = {"id": 2, "delta_energy": None,
         "surface_area": None, "buried_surface_area": None}
        mmcif_dict = {"pdbx_struct_assembly_prop": [{
         "biol_id": "1", "type": "MORE", "value": "0.5"
        }]}
        assign_metrics_to_assembly(assembly, mmcif_dict)
        self.assertEqual(assembly["delta_energy"], None)
        self.assertEqual(assembly["surface_area"], None)
        self.assertEqual(assembly["buried_surface_area"], None)
        assign_metrics_to_assembly(assembly, {})
        self.assertEqual(assembly["delta_energy"], None)
        self.assertEqual(assembly["surface_area"], None)
        self.assertEqual(assembly["buried_surface_area"], None)


    def test_can_update_metrics(self):
        assembly = {"id": 2, "delta_energy": None,
         "surface_area": None, "buried_surface_area": None}
        mmcif_dict = {"pdbx_struct_assembly_prop": [{
         "biol_id": "1", "type": "MORE", "value": "0.5"
        }, {
         "biol_id": "2", "type": "MORE", "value": "0.5"
        }, {
         "biol_id": "2", "type": "SSA (A^2)", "value": "0.3"
        }, {
         "biol_id": "2", "type": "ABSA (A^2)", "value": "0.2"
        }]}
        assign_metrics_to_assembly(assembly, mmcif_dict)
        self.assertEqual(assembly["delta_energy"], 0.5)
        self.assertEqual(assembly["surface_area"], 0.3)
        self.assertEqual(assembly["buried_surface_area"], 0.2)



class AssemblyTransformationsAssigningTests(TestCase):

    @patch("atomium.files.mmcif.get_transformation_ids")
    def test_can_assign_transformations(self, mock_get):
        mock_get.return_value = ["10", "11"]
        assembly = {"id": 2, "transformations": []}
        mmcif_dict = {"pdbx_struct_assembly_gen": [{
         "assembly_id": "1", "oper_expression": "A", "asym_id_list": "A,B"
        }, {
         "assembly_id": "3", "oper_expression": "B", "asym_id_list": "C,D"
        }, {
         "assembly_id": "2", "oper_expression": "C", "asym_id_list": "A,B,C,D"
        }], "pdbx_struct_oper_list": [{
         "id": "10",
         "matrix[1][1]": "0.309", "matrix[1][2]": "-0.5", "matrix[1][3]": "0.809",
         "vector[1]": "-193.90817",
         "matrix[2][1]": "-0.809", "matrix[2][2]": "0.3", "matrix[2][3]": "0.5",
         "vector[2]": "0.00000",
         "matrix[3][1]": "-0.5", "matrix[3][2]": "-0.809", "matrix[3][3]": "-0.3",
         "vector[3]": "507.65814"
        }, {
         "id": "11",
         "matrix[1][1]": "1.0", "matrix[1][2]": "0.0", "matrix[1][3]": "0.0",
         "vector[1]": "0.0",
         "matrix[2][1]": "0.0", "matrix[2][2]": "1.0", "matrix[2][3]": "0.0",
         "vector[2]": "0.0",
         "matrix[3][1]": "0.0", "matrix[3][2]": "0.0", "matrix[3][3]": "1.0",
         "vector[3]": "0.0"
        }, {"id": "12"}]}
        assign_transformations_to_assembly(assembly, mmcif_dict)
        mock_get.assert_called_with("C")
        self.assertEqual(assembly["transformations"], [{
         "chains": ["A", "B", "C", "D"],
         "vector": [-193.90817, 0.0, 507.65814],
         "matrix": [[0.309, -0.5, 0.809], [-0.809, 0.3, 0.5], [-0.5, -0.809, -0.3]]
        }, {
         "chains": ["A", "B", "C", "D"],
         "vector": [0.0, 0.0, 0.0],
         "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        }])



class TransformationIdGettingTests(TestCase):

    def test_can_get_simple_id(self):
        self.assertEqual(get_transformation_ids("1"), ["1"])


    def test_can_get_simple_ids(self):
        self.assertEqual(get_transformation_ids("1,2,3"), ["1", "2", "3"])


    def test_can_get_bracketed_ids(self):
        self.assertEqual(get_transformation_ids("(1,2,3)"), ["1", "2", "3"])


    def test_can_get_id_range(self):
        self.assertEqual(get_transformation_ids("(1-3)"), ["1", "2", "3"])


    def test_can_get_id_multiples(self):
        self.assertEqual(
         get_transformation_ids("(1-3)(A)(9,8)"),
         ["1", "2", "3", "A", "9", "8"]
        )



class AtomLineToAtomDictTests(TestCase):

    def setUp(self):
        self.d = {
         "group_PDB": "ATOM",
         "id": "26962",
         "type_symbol": "O",
         "label_atom_id": "O1",
         "label_alt_id": "A",
         "label_comp_id": "GLY",
         "label_asym_id": "B",
         "label_entity_id": "2",
         "label_seq_id": "13",
         "pdbx_PDB_ins_code": "C",
         "Cartn_x": "199.639",
         "Cartn_y": "91.034",
         "Cartn_z": "-25.211",
         "occupancy": "0.70",
         "B_iso_or_equiv": "11.21",
         "pdbx_formal_charge": "-2",
         "auth_seq_id": "13",
         "auth_comp_id": "TRP",
         "auth_asym_id": "B",
         "auth_atom_id": "O",
         "pdbx_PDB_model_num": "1"
        }


    def test_can_convert_empty_line_to_atom(self):
        for k in self.d: self.d[k] = "."
        atom = atom_line_to_atom_dict(self.d)
        self.assertEqual(atom, deepcopy(ATOM_DICT))


    def test_can_get_atom(self):
        atom = atom_line_to_atom_dict(self.d)
        self.assertEqual(atom, {
         "id": 26962, "name": "O1",
         "alt_loc": "A", "residue_name": "GLY",
         "chain_id": "B", "residue_id": 13, "residue_insert": "C",
         "x": 199.639, "y": 91.034, "z": -25.211,
         "occupancy": 0.7, "bfactor": 11.21, "anisotropy": [],
         "element": "O", "charge": -2, "polymer": True, "full_res_id": "B:13C"
        })


    def test_can_get_hetatm(self):
        self.d["label_seq_id"] = "."
        atom = atom_line_to_atom_dict(self.d)
        self.assertEqual(atom, {
         "id": 26962, "name": "O1",
         "alt_loc": "A", "residue_name": "GLY",
         "chain_id": "B", "residue_id": 13, "residue_insert": "C",
         "x": 199.639, "y": 91.034, "z": -25.211,
         "occupancy": 0.7, "bfactor": 11.21, "anisotropy": [],
         "element": "O", "charge": -2, "polymer": False, "full_res_id": "B:13C"
        })
