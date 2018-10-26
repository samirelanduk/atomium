from collections import deque
from datetime import date
import numpy as np
from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.mmcif import *

class MmcifStringToMmcifDictTests(TestCase):

    @patch("atomium.mmcif.consolidate_strings")
    @patch("atomium.mmcif.mmcif_lines_to_mmcif_blocks")
    @patch("atomium.mmcif.loop_block_to_list")
    @patch("atomium.mmcif.non_loop_block_to_list")
    @patch("atomium.mmcif.strip_quotes")
    def test_can_turn_mmcif_string_to_mmcif_dict(self, mock_strp, mock_nnlp, mock_loop, mock_block, mock_con):
        mock_con.return_value = "CONSOLIDATED"
        mock_block.return_value = [
         {"category": "A", "lines": ["1", "2"]},
         {"category": "B", "lines": ["loop_", "3"]},
         {"category": "C", "lines": ["loop_", "4"]},
         {"category": "D", "lines": ["5", "6"]},
        ]
        mock_nnlp.side_effect = ["CAT1", "CAT2"]
        mock_loop.side_effect = ["LOOP1", "LOOP2"]
        d = mmcif_string_to_mmcif_dict("1\n\n2\n3\n")
        mock_con.assert_called_with(deque(["1", "2", "3"]))
        mock_block.assert_called_with("CONSOLIDATED")
        mock_loop.assert_any_call({"category": "B", "lines": ["loop_", "3"]})
        mock_loop.assert_any_call({"category": "C", "lines": ["loop_", "4"]})
        mock_nnlp.assert_any_call({"category": "A", "lines": ["1", "2"]})
        mock_nnlp.assert_any_call({"category": "D", "lines": ["5", "6"]})
        mock_strp.assert_called_with({
         "A": "CAT1", "B": "LOOP1", "C": "LOOP2", "D": "CAT2"
        })
        self.assertEqual(d, {
         "A": "CAT1", "B": "LOOP1", "C": "LOOP2", "D": "CAT2"
        })



class StringConsolidationTests(TestCase):

    def test_consolidation_can_do_nothing(self):
        lines = consolidate_strings(deque(["line1", "line2", "line3"]))
        self.assertEqual(lines, deque(["line1", "line2", "line3"]))


    def test_can_consolidate_string(self):
        lines = consolidate_strings(deque(["line1", "line2", ";STRING", ";", "line3"]))
        self.assertEqual(lines, deque(["line1", "line2 \"STRING\"", "line3"]))


    def test_can_consolidate_multi_strings(self):
        lines = consolidate_strings(deque([
         "line1", "line2", "; STRING 1", "CONT CONT", ";", "line3", ";STRING2", ";", "line4"
        ]))
        self.assertEqual(lines, deque([
         "line1", "line2 \"STRING 1 CONT CONT\"", "line3 \"STRING2\"", "line4"
        ]))



class MmcifLinesToMmcifBlockTests(TestCase):

    def test_can_make_non_loop_blocks(self):
        blocks = mmcif_lines_to_mmcif_blocks(deque([
         "data_1", "_AAA.x 100", "_BBB.x 200", "_BBB.y 200", "_CCC.z 300"
        ]))
        self.assertEqual(blocks, [
         {"category": "AAA", "lines": ["_AAA.x 100"]},
         {"category": "BBB", "lines": ["_BBB.x 200", "_BBB.y 200"]},
         {"category": "CCC", "lines": ["_CCC.z 300"]}
        ])


    def test_can_make_loop_blocks(self):
        blocks = mmcif_lines_to_mmcif_blocks(deque([
         "data_1", "loop_", "_AAA.x", "100", "loop_", "_BBB.x", "_BBB.y", "11", "12", "13"
        ]))
        self.assertEqual(blocks, [
         {"category": "AAA", "lines": ["loop_", "_AAA.x", "100"]},
         {"category": "BBB", "lines": ["loop_", "_BBB.x", "_BBB.y", "11", "12", "13"]},
        ])


    def test_can_make_mixed_blocks(self):
        blocks = mmcif_lines_to_mmcif_blocks(deque([
         "data_1", "_AAA.x 100", "loop_", "_111.x", "100", "_BBB.x 200",
         "_BBB.y 200", "loop_", "_222.x", "_222.y", "11", "12", "13", "_CCC.z 300"
        ]))
        self.assertEqual(blocks, [
         {"category": "AAA", "lines": ["_AAA.x 100"]},
          {"category": "111", "lines": ["loop_", "_111.x", "100"]},
         {"category": "BBB", "lines": ["_BBB.x 200", "_BBB.y 200"]},
          {"category": "222", "lines": ["loop_", "_222.x", "_222.y", "11", "12", "13"]},
         {"category": "CCC", "lines": ["_CCC.z 300"]}
        ])



class NonLoopBlockToListTests(TestCase):

    def test_can_convert_non_loop_block_to_list(self):
        l = non_loop_block_to_list({"category": "AA", "lines": [
         "_AA.one XXX", "_AA.two 'YYY ZZZ'", "_AA.three 1.2.3"
        ]})
        self.assertEqual(l, [{
         "one": "XXX", "two": "'YYY ZZZ'", "three": "1.2.3"
        }])



class LoopBlockToListTests(TestCase):

    @patch("atomium.mmcif.split_values")
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


    @patch("atomium.mmcif.split_values")
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
         "A": [{
          "1": "20", "2": '"34"'
         }], "B": [{
          "3": "30", "4": '"34"'
         }, {
          "3": "50", "4": "'34 dfe'"
         }]
        }
        strip_quotes(d)
        self.assertEqual(d, {
         "A": [{
          "1": "20", "2": '34'
         }], "B": [{
          "3": "30", "4": '34'
         }, {
          "3": "50", "4": '34 dfe'
         }]
        })



class MmcifDictToDataDictTests(TestCase):

    @patch("atomium.mmcif.update_description_dict")
    @patch("atomium.mmcif.update_experiment_dict")
    @patch("atomium.mmcif.update_quality_dict")
    @patch("atomium.mmcif.update_geometry_dict")
    @patch("atomium.mmcif.update_models_list")
    def test_can_convert_mmcif_dict_to_data_dict(self, mock_md, mock_gm, mock_ql, mock_ex, mock_ds):
        mmcif_dict = {"A": "B"}
        d = mmcif_dict_to_data_dict(mmcif_dict)
        mock_ds.assert_called_with(mmcif_dict, d)
        mock_ex.assert_called_with(mmcif_dict, d)
        mock_ql.assert_called_with(mmcif_dict, d)
        mock_gm.assert_called_with(mmcif_dict, d)
        mock_md.assert_called_with(mmcif_dict, d)
        self.assertEqual(d, {
         "description": {
          "code": None, "title": None, "deposition_date": None,
          "classification": None, "keywords": [], "authors": []
         }, "experiment": {
          "technique": None, "source_organism": None, "expression_system": None
         }, "quality": {"resolution": None, "rvalue": None, "rfree": None},
         "geometry": {"assemblies": []}, "models": []
        })



class DescriptionDictionaryUpdatingTests(TestCase):

    @patch("atomium.mmcif.mmcif_to_data_transfer")
    def test_can_update_description_dictionary(self, mock_trans):
        m, d = {"M": 1}, {"D": 2}
        update_description_dict(m, d)
        mock_trans.assert_any_call(m, d, "description", "code", "entry", "id")
        mock_trans.assert_any_call(m, d, "description", "title", "struct", "title")
        mock_trans.assert_any_call(m, d, "description", "deposition_date", "pdbx_database_status", "recvd_initial_deposition_date", date=True)
        mock_trans.assert_any_call(m, d, "description", "classification", "struct_keywords", "pdbx_keywords")
        mock_trans.assert_any_call(m, d, "description", "keywords", "struct_keywords", "text", split=True)
        mock_trans.assert_any_call(m, d, "description", "authors", "audit_author", "name", multi=True)



class ExperimentDictionaryUpdatingTests(TestCase):

    @patch("atomium.mmcif.mmcif_to_data_transfer")
    def test_can_update_experiment_dictionary(self, mock_trans):
        m, d = {"M": 1}, {"D": 2}
        update_experiment_dict(m, d)
        mock_trans.assert_any_call(m, d, "experiment", "technique", "exptl", "method")
        mock_trans.assert_any_call(m, d, "experiment", "source_organism", "entity_src_gen", "pdbx_gene_src_scientific_name")
        mock_trans.assert_any_call(m, d, "experiment", "expression_system", "entity_src_gen", "pdbx_host_org_scientific_name")



class QualityDictionaryUpdatingTests(TestCase):

    @patch("atomium.mmcif.mmcif_to_data_transfer")
    def test_can_update_quality_dictionary(self, mock_trans):
        m, d = {"M": 1}, {"D": 2}
        update_quality_dict(m, d)
        mock_trans.assert_any_call(m, d, "quality", "resolution", "reflns", "d_resolution_high", func=float)
        mock_trans.assert_any_call(m, d, "quality", "rvalue", "refine", "ls_R_factor_R_work", func=float)
        mock_trans.assert_any_call(m, d, "quality", "rfree", "refine", "ls_R_factor_R_free", func=float)



class GeometryDictUpdatingTests(TestCase):

    def test_can_update_geometry_with_nothing(self):
        m, d = {}, {"geometry": {"assemblies": []}}
        update_geometry_dict(m, d)
        self.assertEqual(d, {"geometry": {"assemblies": []}})


    @patch("atomium.mmcif.assign_metrics_to_assembly")
    @patch("atomium.mmcif.assign_transformations_to_assembly")
    def test_can_add_assemblies_to_geometry(self, mock_trans, mock_ass):
        d = {"geometry": {"assemblies": []}}
        m = {"pdbx_struct_assembly": [{
         "id": "1", "method_details": "PISA",
        }, {
         "id": "2", "method_details": "?",
        }], "pdbx_struct_oper_list": [{
         "id": "2", "matrix[1][1]": "0.80901699", "matrix[1][2]": "0.30901699",
         "matrix[1][3]": "0.50000000", "vector[1]": "0.00000", "matrix[2][1]": "0.30901699",
         "matrix[2][2]": "0.50000000", "matrix[2][3]": "-0.80901699", "vector[2]": "0.00000",
         "matrix[3][1]": "-0.50000000", "matrix[3][2]": "0.80901699",
         "matrix[3][3]": "0.30901699", "vector[3]": "0.00000"
        }, {
         "id": "9", "matrix[1][1]": "0.5", "matrix[1][2]": "0.3",
         "matrix[1][3]": "0.6", "vector[1]": "0.2", "matrix[2][1]": "0.6",
         "matrix[2][2]": "0.5", "matrix[2][3]": "-0.5", "vector[2]": "0.0",
         "matrix[3][1]": "-0.1", "matrix[3][2]": "0.2",
         "matrix[3][3]": "0.2", "vector[3]": "0.6"
        }]}
        update_geometry_dict(m, d)
        self.assertEqual(d, {"geometry": {"assemblies": [{
         "id": 1, "software": "PISA", "delta_energy": None,
         "buried_surface_area": None, "surface_area": None, "transformations": []
        }, {
         "id": 2, "software": None, "delta_energy": None,
         "buried_surface_area": None, "surface_area": None, "transformations": []
        }]}})
        operations = {"2": [
         [0.80901699, 0.30901699, 0.5, 0.0], [0.30901699, 0.5, -0.80901699, 0.0],
         [-0.5, 0.80901699, 0.30901699, 0.0], [0, 0, 0, 1]
        ], "9": [
         [0.5, 0.3, 0.6, 0.2], [0.6, 0.5, -0.5, 0.0],
         [-0.1, 0.2, 0.2, 0.6], [0, 0, 0, 1]
        ]}
        for assembly in d["geometry"]["assemblies"]:
            mock_ass.assert_any_call(m, assembly)
            mock_trans.assert_any_call(m, operations, assembly)



class MetricsAssigningTests(TestCase):

    def setUp(self):
        self.assembly = {
         "id": 1, "software": "PISA", "delta_energy": None,
         "buried_surface_area": None, "surface_area": None, "transformations": []
        }


    def test_can_assign_nothing(self):
        assign_metrics_to_assembly({}, self.assembly)
        self.assertEqual(self.assembly, {
         "id": 1, "software": "PISA", "delta_energy": None,
         "buried_surface_area": None, "surface_area": None, "transformations": []
        })


    def test_can_assign_metrics_to_assembly(self):
        m = {"pdbx_struct_assembly_prop": [{
         "biol_id": "1", "type": "ABSA (A^2)", "value": "1720",
        }, {
         "biol_id": "1", "type": "MORE", "value": "-7",
        }, {
         "biol_id": "1", "type": "SSA (A^2)", "value": "3980",
        }, {
         "biol_id": "1", "type": "XXX", "value": "123",
        }, {
         "biol_id": "2", "type": "SSA (A^2)", "value": "9999",
        }]}
        assign_metrics_to_assembly(m, self.assembly)
        self.assertEqual(self.assembly, {
         "id": 1, "software": "PISA", "delta_energy": -7, "surface_area": 3980,
         "buried_surface_area": 1720, "transformations": []
        })



class TransformationAssigningTests(TestCase):

    def setUp(self):
        self.assembly = {
         "id": 2, "software": "PISA", "delta_energy": None,
         "buried_surface_area": None, "surface_area": None, "transformations": []
        }


    def test_can_assign_nothing(self):
        assign_transformations_to_assembly({}, {}, self.assembly)
        self.assertEqual(self.assembly, {
         "id": 2, "software": "PISA", "delta_energy": None,
         "buried_surface_area": None, "surface_area": None, "transformations": []
        })


    @patch("atomium.mmcif.get_operation_id_groups")
    @patch("atomium.mmcif.operation_id_groups_to_operations")
    def test_can_assign_transformations(self, mock_convert, mock_get):
        mock_get.return_value = [[1], [2]]
        mock_convert.return_value = [["abcd", "efgh", "ijkl", "mnop"], ["ABCD", "EFGH", "IJKL", "MNOP"]]
        m = {"pdbx_struct_assembly_gen": [{
         "assembly_id": "1", "oper_expression": "A", "asym_id_list": "A,B"
        }, {
         "assembly_id": "3", "oper_expression": "B", "asym_id_list": "C,D"
        }, {
         "assembly_id": "2", "oper_expression": "C", "asym_id_list": "A,B,C,D"
        }]}
        operations = {"1": [
         [0.80901699, 0.30901699, 0.5, 0.0], [0.30901699, 0.5, -0.80901699, 0.0],
         [-0.5, 0.80901699, 0.30901699, 0.0], [0, 0, 0, 1]
        ], "2": [
         [0.5, 0.3, 0.6, 0.2], [0.6, 0.5, -0.5, 0.0],
         [-0.1, 0.2, 0.2, 0.6], [0, 0, 0, 1]
        ]}
        assign_transformations_to_assembly(m, operations, self.assembly)
        mock_get.assert_called_with("C")
        mock_convert.assert_called_with(operations, [[1], [2]])
        self.assertEqual(self.assembly["transformations"], [{
         "chains": ["A", "B", "C", "D"],
         "vector": list("dhl"),
         "matrix": ["abc", "efg", "ijk"]
        }, {
         "chains": ["A", "B", "C", "D"],
         "vector": list("DHL"),
         "matrix": ["ABC", "EFG", "IJK"]
        }])



class OperationIdGroupsGettingTests(TestCase):

    def test_can_get_simple_id(self):
        self.assertEqual(get_operation_id_groups("1"), [["1"]])
        self.assertEqual(get_operation_id_groups("(1)"), [["1"]])


    def test_can_get_simple_ids(self):
        self.assertEqual(get_operation_id_groups("1,2,3"), [["1", "2", "3"]])
        self.assertEqual(get_operation_id_groups("(1,2,3)"), [["1", "2", "3"]])


    def test_can_get_id_range(self):
        self.assertEqual(get_operation_id_groups("1-3"), [["1", "2", "3"]])
        self.assertEqual(get_operation_id_groups("(1-3)"), [["1", "2", "3"]])


    def test_can_get_mixed_range(self):
        self.assertEqual(get_operation_id_groups("1-3,5,7"), [["1", "2", "3", "5", "7"]])
        self.assertEqual(get_operation_id_groups("(1-3,5,7)"), [["1", "2", "3", "5", "7"]])


    def test_can_get_id_multiples(self):
        self.assertEqual(
         get_operation_id_groups("(1,2)(3,4)"), [["1", "2"], ["3", "4"]]
        )
        self.assertEqual(
         get_operation_id_groups("(1-3)(3,4)"), [["1", "2", "3"], ["3", "4"]]
        )



class OperationIdGroupsToOperationsTests(TestCase):

    def test_can_covert_single_group(self):
        operations = {1: [[1, 2, 3, 4], [5, 6, 7, 8], [9, 1, 2, 3], [4, 5, 6, 7]], 2: []}
        self.assertEqual(
         operation_id_groups_to_operations(operations, [[1]]),
         [operations[1]]
        )


    def test_can_covert_multiple_groups(self):
        operations = {1: [[1, 2, 3, 4], [5, 6, 7, 8], [9, 1, 2, 3], [4, 5, 6, 7]], 2: []}
        self.assertEqual(
         operation_id_groups_to_operations(operations, [[1, 2]]),
         [operations[1], operations[2]]
        )


    def test_can_multiply_groups(self):
        operations = {
         1: [[1, 2, 3, 4], [5, 6, 7, 8], [9, 1, 2, 3], [4, 5, 6, 7]],
         2: [[1, 2, 3, 9], [5, 6, 7, 8], [9, 1, 9, 3], [4, 5, 6, 7]],
         3: [[1, 2, 3, 9], [5, 6, 7, 8], [9, 1, 2, 3], [4, 5, 6, 7]],
         4: [[1, 2, 3, 9], [5, 6, 5, 8], [9, 1, 2, 3], [4, 5, 6, 7]]
        }
        self.assertTrue(np.array_equal(
         operation_id_groups_to_operations(operations, [[1, 2], [3, 4]]), [
         np.array([[54, 37, 47, 62], [130, 93, 119, 170], [44, 41, 56, 116], [111, 79, 101, 143]]),
         np.array([[54, 37, 43, 62], [130, 93, 107, 170], [44, 41, 54, 116], [111, 79, 91, 143]]),
         np.array([[74, 62, 77, 97], [130, 93, 119, 170], [107, 48, 70, 137], [111, 79, 101, 143]]),
         np.array([[74, 62, 73, 97], [130, 93, 107, 170], [107, 48, 68, 137], [111, 79, 91, 143]])
        ]))



class ModelsListUpdatingTests(TestCase):

    @patch("atomium.mmcif.make_sequences")
    @patch("atomium.mmcif.make_aniso")
    @patch("atomium.mmcif.add_atom_to_polymer")
    @patch("atomium.mmcif.add_atom_to_non_polymer")
    @patch("atomium.mmcif.add_sequences_to_polymers")
    def test_can_update_one_model(self, mock_sq, mock_np, mock_p, mock_an, mock_ms):
        m = {"entity": [
         {"id": "1", "type": "polymer"}, {"id": "2", "type": "non-polymer"},
         {"id": "3", "type": "non-polymer"}, {"id": "4", "type": "water"},
         {"id": "5", "type": "water"}
        ], "struct_asym": [
         {"id": "A", "entity_id": "1"}, {"id": "B", "entity_id": "1"},
         {"id": "C", "entity_id": "2"}, {"id": "D", "entity_id": "3"},
         {"id": "E", "entity_id": "2"}, {"id": "F", "entity_id": "3"},
         {"id": "G", "entity_id": "4"}, {"id": "H", "entity_id": "5"}
        ], "atom_site": [
         {"pdbx_PDB_model_num": "1", "label_asym_id": n} for n in "ABCDEFGH"
        ]}
        d = {"models": []}
        model = {"polymer": {}, "non-polymer": {}, "water": {}}
        update_models_list(m, d)
        mock_ms.assert_called_with(m)
        mock_an.assert_called_with(m)
        mock_p.assert_any_call(m["atom_site"][0], mock_an.return_value, model)
        mock_p.assert_any_call(m["atom_site"][1], mock_an.return_value, model)
        mock_np.assert_any_call(m["atom_site"][2], mock_an.return_value, model, "non-polymer")
        mock_np.assert_any_call(m["atom_site"][3], mock_an.return_value, model, "non-polymer")
        mock_np.assert_any_call(m["atom_site"][4], mock_an.return_value, model, "non-polymer")
        mock_np.assert_any_call(m["atom_site"][5], mock_an.return_value, model, "non-polymer")
        mock_np.assert_any_call(m["atom_site"][6], mock_an.return_value, model, "water")
        mock_np.assert_any_call(m["atom_site"][7], mock_an.return_value, model, "water")
        self.assertEqual(d["models"], [model])


    @patch("atomium.mmcif.make_sequences")
    @patch("atomium.mmcif.make_aniso")
    @patch("atomium.mmcif.add_atom_to_polymer")
    @patch("atomium.mmcif.add_atom_to_non_polymer")
    @patch("atomium.mmcif.add_sequences_to_polymers")
    def test_can_update_multiple_model(self, mock_sq, mock_np, mock_p, mock_an, mock_ms):
        m = {"entity": [
         {"id": "1", "type": "polymer"}, {"id": "2", "type": "non-polymer"},
         {"id": "3", "type": "non-polymer"}, {"id": "4", "type": "water"},
         {"id": "5", "type": "water"}
        ], "struct_asym": [
         {"id": "A", "entity_id": "1"}, {"id": "B", "entity_id": "1"},
         {"id": "C", "entity_id": "2"}, {"id": "D", "entity_id": "3"},
         {"id": "E", "entity_id": "2"}, {"id": "F", "entity_id": "3"},
         {"id": "G", "entity_id": "4"}, {"id": "H", "entity_id": "5"}
        ], "atom_site": [
         {"pdbx_PDB_model_num": m, "label_asym_id": n} for m in "123" for n in "ABCDEFGH"
        ]}
        d = {"models": []}
        model = {"polymer": {}, "non-polymer": {}, "water": {}}
        update_models_list(m, d)
        mock_ms.assert_called_with(m)
        mock_an.assert_called_with(m)
        mock_p.assert_any_call(m["atom_site"][0], mock_an.return_value, model)
        mock_p.assert_any_call(m["atom_site"][1], mock_an.return_value, model)
        mock_np.assert_any_call(m["atom_site"][2], mock_an.return_value, model, "non-polymer")
        mock_np.assert_any_call(m["atom_site"][3], mock_an.return_value, model, "non-polymer")
        mock_np.assert_any_call(m["atom_site"][4], mock_an.return_value, model, "non-polymer")
        mock_np.assert_any_call(m["atom_site"][5], mock_an.return_value, model, "non-polymer")
        mock_np.assert_any_call(m["atom_site"][6], mock_an.return_value, model, "water")
        mock_np.assert_any_call(m["atom_site"][7], mock_an.return_value, model, "water")
        mock_sq.assert_called_with(model, m, {
         "A": "1", "B": "1","C": "2", "D": "3", "E": "2", "F": "3", "G": "4", "H": "5"
        })
        self.assertEqual(d["models"], [model, model, model])



class AnisoMakingTests(TestCase):

    def test_can_handle_no_aniso(self):
        self.assertEqual(make_aniso({}), {})


    def test_can_make_aniso(self):
        self.assertEqual(make_aniso({"atom_site_anisotrop": [{
         "id": "1", "U[1][1]": "0.98", "U[2][2]": "0.748", "U[3][3]": "0.4152",
         "U[1][2]": "-0.1159", "U[1][3]": "-0.0115", "U[2][3]": "-0.2655",
        }, {
         "id": "9", "U[1][1]": "0.909", "U[2][2]": "0.548", "U[3][3]": "0.4152",
         "U[1][2]": "-0.2159", "U[1][3]": "-0.115", "U[2][3]": "-0.2654",
        }]}), {
         1: [0.98, 0.748, 0.4152, -0.1159, -0.0115, -0.2655],
         9: [0.909, 0.548, 0.4152, -0.2159, -0.115, -0.2654]
        })



class AtomToPolymerAddingTests(TestCase):

    @patch("atomium.mmcif.make_residue_id")
    @patch("atomium.mmcif.atom_dict_to_atom_dict")
    def test_can_add_atom_to_existing_residue(self, mock_at, mock_id):
        atom = {"id": "10", "auth_asym_id": "C"}
        mock_id.return_value = "A.10"
        model = {"polymer": {"C": {"residues": {"A.10": {"atoms": {}}}}}}
        add_atom_to_polymer(atom, "ANISO", model)
        mock_id.assert_called_with(atom)
        mock_at.assert_called_with(atom, "ANISO")
        self.assertEqual(model, {"polymer": {"C": {
         "residues": {"A.10": {"atoms": {10: mock_at.return_value}}}
        }}})


    @patch("atomium.mmcif.make_residue_id")
    @patch("atomium.mmcif.atom_dict_to_atom_dict")
    def test_can_add_atom_to_existing_chain(self, mock_at, mock_id):
        self.maxDiff = None
        atom = {"id": "10", "auth_asym_id": "C", "auth_comp_id": "VAL"}
        mock_id.return_value = "A.10"
        model = {"polymer": {"C": {"residues": {}}}}
        add_atom_to_polymer(atom, "ANISO", model)
        mock_id.assert_called_with(atom)
        mock_at.assert_called_with(atom, "ANISO")
        self.assertEqual(model, {"polymer": {"C": {
         "residues": {"A.10": {"atoms": {10: mock_at.return_value}, "name": "VAL", "number": 1}},
        }}})


    @patch("atomium.mmcif.make_residue_id")
    @patch("atomium.mmcif.atom_dict_to_atom_dict")
    def test_can_add_atom_to_non_existing_chain(self, mock_at, mock_id):
        self.maxDiff = None
        atom = {"id": "10", "auth_asym_id": "C", "auth_comp_id": "VAL", "label_asym_id": "F"}
        mock_id.return_value = "A.10"
        model = {"polymer": {}}
        add_atom_to_polymer(atom, "ANISO", model)
        mock_id.assert_called_with(atom)
        mock_at.assert_called_with(atom, "ANISO")
        self.assertEqual(model, {"polymer": {"C": {"internal_id": "F",
         "residues": {"A.10": {"atoms": {10: mock_at.return_value}, "name": "VAL", "number": 1}},
        }}})



class AtomToNonPolymerAddingTests(TestCase):

    @patch("atomium.mmcif.make_residue_id")
    @patch("atomium.mmcif.atom_dict_to_atom_dict")
    def test_can_add_atom_to_existing_molecule(self, mock_at, mock_id):
        atom = {"id": "10", "auth_asym_id": "C"}
        mock_id.return_value = "A.10"
        model = {"TYPE": {"A.10": {"atoms": {}}}}
        add_atom_to_non_polymer(atom, "ANISO", model, "TYPE")
        mock_id.assert_called_with(atom)
        mock_at.assert_called_with(atom, "ANISO")
        self.assertEqual(
         model, {"TYPE": {"A.10": {"atoms": {10: mock_at.return_value}}}}
        )


    @patch("atomium.mmcif.make_residue_id")
    @patch("atomium.mmcif.atom_dict_to_atom_dict")
    def test_can_add_atom_to_non_existing_molecule(self, mock_at, mock_id):
        atom = {"id": "10", "auth_asym_id": "C", "auth_comp_id": "XYZ", "label_asym_id": "A"}
        mock_id.return_value = "A.10"
        model = {"TYPE": {}}
        add_atom_to_non_polymer(atom, "ANISO", model, "TYPE")
        mock_id.assert_called_with(atom)
        mock_at.assert_called_with(atom, "ANISO")
        self.assertEqual(
         model, {"TYPE": {"A.10": {
          "internal_id": "A", "polymer": "C", "name": "XYZ", "atoms": {10: mock_at.return_value}
         }}}
        )



class ResidueIdMakingTests(TestCase):

    def test_can_make_residue_id(self):
        self.assertEqual(make_residue_id({
         "auth_asym_id": "A", "auth_seq_id": 100, "pdbx_PDB_ins_code": "B"
        }), "A.100B")


    def test_can_make_residue_id_with_missing_info(self):
        self.assertEqual(make_residue_id({
         "auth_asym_id": "A", "auth_seq_id": 100, "pdbx_PDB_ins_code": "?"
        }), "A.100")
        self.assertEqual(make_residue_id({
         "auth_asym_id": "A", "auth_seq_id": 100, "pdbx_PDB_ins_code": "."
        }), "A.100")



class SequencesToPolymerAdditionTests(TestCase):

    @patch("atomium.mmcif.make_sequences")
    def test_can_add_no_sequences(self, mock_make):
        mock_make.return_value = {}
        model = {"polymer": {"A": {"internal_id": "T"}, "B": {"internal_id": "X"}}, "non-polmer": {}}
        add_sequences_to_polymers(model, {1: 2}, {"A": "1", "B": "2"})
        self.assertEqual(model, {"polymer": {
         "A": {"internal_id": "T", "sequence": ""}, "B": {"internal_id": "X", "sequence": ""}
        }, "non-polmer": {}})


    @patch("atomium.mmcif.make_sequences")
    def test_can_add_sequences(self, mock_make):
        mock_make.return_value = {"1": "ABC", "2": "DEFG"}
        model = {"polymer": {"A": {"internal_id": "T"}, "B": {"internal_id": "X"}}, "non-polmer": {}}
        add_sequences_to_polymers(model, {1: 2}, {"T": "1", "X": "2"})
        self.assertEqual(model, {"polymer": {
         "A": {"internal_id": "T", "sequence": "ABC"}, "B": {"internal_id": "X", "sequence": "DEFG"}
        }, "non-polmer": {}})



class SequenceMakingTests(TestCase):

    def test_can_make_no_sequences(self):
        self.assertEqual(make_sequences({}), {})
        self.assertEqual(make_sequences({
         "entity": [{"type": "non-polymer"}, {"type": "water"}]
        }), {})
        self.assertEqual(make_sequences({
         "entity": [{"id": "1", "type": "polymer"}]
        }), {"1": ""})


    def test_can_make_sequences(self):
        self.assertEqual(make_sequences({
         "entity": [{"id": "1", "type": "polymer"}, {"id": "2", "type": "polymer"}],
         "entity_poly_seq": [
          {"entity_id": "1", "mon_id": "MET"}, {"entity_id": "1", "mon_id": "XXX"},
          {"entity_id": "2", "mon_id": "DA"}, {"entity_id": "2", "mon_id": "DT"}, {"entity_id": "2", "mon_id": "DC"}
         ]
        }), {"1": "MX", "2": "ATC"})



class AtomDictToAtomDictTests(TestCase):

    def setUp(self):
        self.d = {
         "group_PDB": "ATOM", "id": "26962", "type_symbol": "O",
         "label_atom_id": "O1", "label_alt_id": "A", "label_comp_id": "GLY",
         "label_asym_id": "B", "label_entity_id": "2", "label_seq_id": "13",
         "pdbx_PDB_ins_code": "C", "Cartn_x": "199.639", "Cartn_y": "91.034",
         "Cartn_z": "-25.211", "occupancy": "0.70", "B_iso_or_equiv": "11.21",
         "pdbx_formal_charge": "-2", "auth_seq_id": "13", "auth_comp_id": "TRP",
         "auth_asym_id": "B", "auth_atom_id": "O", "pdbx_PDB_model_num": "1"
        }


    def test_can_get_atom(self):
        atom = atom_dict_to_atom_dict(self.d, {})
        self.assertEqual(atom, {
         "name": "O1", "alt_loc": "A",
         "x": 199.639, "y": 91.034, "z": -25.211,
         "occupancy": 0.7, "bvalue": 11.21, "anisotropy": [0] * 6,
         "element": "O", "charge": -2
        })


    def test_can_get_atom_with_anisotropy(self):
        self.d["label_alt_id"] = "."
        atom = atom_dict_to_atom_dict(self.d, {1: 100, 26962: 200})
        self.assertEqual(atom, {
         "name": "O1", "alt_loc": None,
         "x": 199.639, "y": 91.034, "z": -25.211,
         "occupancy": 0.7, "bvalue": 11.21, "anisotropy": 200,
         "element": "O", "charge": -2
        })



class MmcifDictTransferTests(TestCase):

    def setUp(self):
        self.d = {"A": {1: None, 2: None, 3: None}, "B": {4: None, 5: None}}
        self.m = {"M": [{10: "ten", 11: "2018-09-17"}], "N": [{12: "tw"}, {12: "tw2"}]}


    def test_can_do_nothing(self):
        mmcif_to_data_transfer(self.m, self.d, "A", 1, "X", 10)
        self.assertEqual(self.d["A"][1], None)
        mmcif_to_data_transfer(self.m, self.d, "A", 1, "M", 100)
        self.assertEqual(self.d["A"][1], None)


    def test_can_transfer_from_mmcif_to_data_dict(self):
        mmcif_to_data_transfer(self.m, self.d, "A", 1, "M", 10)
        self.assertEqual(self.d["A"][1], "ten")
        mmcif_to_data_transfer(self.m, self.d, "B", 5, "N", 12)
        self.assertEqual(self.d["B"][5], "tw")


    def test_can_transfer_from_mmcif_to_data_dict_multi(self):
        mmcif_to_data_transfer(self.m, self.d, "B", 5, "N", 12, multi=True)
        self.assertEqual(self.d["B"][5], ["tw", "tw2"])


    def test_can_transfer_from_mmcif_to_data_dict_date(self):
        mmcif_to_data_transfer(self.m, self.d, "B", 5, "M", 11, date=True)
        self.assertEqual(self.d["B"][5], date(2018, 9, 17))


    def test_can_transfer_from_mmcif_to_data_dict_split(self):
        self.m["M"][0][11] = "5,6, 7"
        mmcif_to_data_transfer(self.m, self.d, "B", 5, "M", 11, split=True)
        self.assertEqual(self.d["B"][5], ["5", "6", "7"])


    def test_can_transfer_from_mmcif_to_data_dict_func(self):
        self.m["M"][0][11] = "102"
        mmcif_to_data_transfer(self.m, self.d, "B", 5, "M", 11, func=int)
        self.assertEqual(self.d["B"][5], 102)
