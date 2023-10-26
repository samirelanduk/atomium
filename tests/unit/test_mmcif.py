from unittest import TestCase
from unittest.mock import patch
from atomium.mmcif import *

class MmcifStringToDictTests(TestCase):

    @patch("atomium.mmcif.get_sections_from_filestring")
    @patch("atomium.mmcif.add_loop_section")
    @patch("atomium.mmcif.add_non_loop_section")
    def test_can_get_dict(self, mock_non, mock_loop, mock_sections):
        mock_sections.return_value = [[""], ["loop_"], [""]]
        d = mmcif_string_to_mmcif_dict("filestring\n\n")
        mock_sections.assert_called_with("filestring")
        mock_loop.assert_called_with(["loop_"], d)
        mock_non.assert_any_call([""], d)
        mock_non.assert_any_call([""], d)



class SectionsFromFilestringTests(TestCase):

    def test_can_get_sections(self):
        filestring = ("data_2BFB\n"
        "# \n"
        "_entry.id   2BFB\n" 
        "# \n"
        "loop_\n"
        "_database_2.database_id\n"
        "_database_2.database_code\n"
        "PDB   2BFB\n"
        "PDBE  EBI-21880\n"
        "# ")
        sections = get_sections_from_filestring(filestring)
        self.assertEqual(sections, [[
            "_entry.id   2BFB"
        ], [
            "loop_", "_database_2.database_id", "_database_2.database_code",
            "PDB   2BFB", "PDBE  EBI-21880"
        ]])



class NonLoopSectionTests(TestCase):

    def test_single_row(self):
        mmcif = {}
        add_non_loop_section(["_entry.id   2BFB"], mmcif)
        self.assertEqual(mmcif["entry"], [{"id": "2BFB"}])
    

    def test_multi_row(self):
        mmcif = {}
        add_non_loop_section([
            "_audit_conform.dict_name       mmcif_pdbx.dic",
            "_audit_conform.dict_version    5.347"
        ], mmcif)
        self.assertEqual(mmcif["audit_conform"], [
            {"dict_name": "mmcif_pdbx.dic", "dict_version": "5.347"}
        ])
    

    def test_blank_lines_row(self):
        mmcif = {}
        add_non_loop_section([
            "_audit_conform.dict_name       mmcif_pdbx.dic",
            "\n ",
            "_audit_conform.dict_version    5.347"
        ], mmcif)
        self.assertEqual(mmcif["audit_conform"], [
            {"dict_name": "mmcif_pdbx.dic", "dict_version": "5.347"}
        ])
    

    def test_spaces_in_values(self):
        mmcif = {}
        add_non_loop_section([
            "_refine.pdbx_refine_id 'X-RAY  DIFFRACTION'",
            "_refine.name       \"XANTHOSINE-5'-MONOPHOSPHATE\""
        ], mmcif)
        self.assertEqual(mmcif["refine"], [
            {"pdbx_refine_id": "X-RAY  DIFFRACTION", "name": "XANTHOSINE-5'-MONOPHOSPHATE"}
        ])
    

    def test_quote_strings(self):
        mmcif = {}
        add_non_loop_section(["_atoms.name1 \"N1'\"", "_atoms.name2 'N1\"'"], mmcif)
        self.assertEqual(mmcif["atoms"], [{"name1": "N1'", "name2": "N1\""}])
    

    def test_value_on_next_line_with_quotes(self):
        mmcif = {}
        add_non_loop_section([
            "_cat.key1 1234",
            "_cat.key2", "'5678'",
            "_cat.key3 ABCD",
            "_cat.key4", "\"EFGHIJKL\"",
        ], mmcif)
        self.assertEqual(mmcif, {"cat": [
            {"key1": "1234", "key2": "5678", "key3": "ABCD", "key4": "EFGHIJKL"}
        ]})
    

    def test_value_on_next_line_without_quotes(self):
        mmcif = {}
        add_non_loop_section([
            "_cat.key1 1234",
            "_cat.key2", "5678",
            "_cat.key3 ABCD",
            "_cat.key4", "\"EFGHIJKL\"",
        ], mmcif)
        self.assertEqual(mmcif, {"cat": [
            {"key1": "1234", "key2": "5678", "key3": "ABCD", "key4": "EFGHIJKL"}
        ]})


    @patch("atomium.mmcif.get_semicolon_value")
    def test_semicolon_value_on_next_lines(self, mock_semi):
        mmcif = {}
        mock_semi.side_effect = lambda d: d.popleft() and "VAL"
        add_non_loop_section([
            "_cat.key1 1234",
            "_cat.key2", ";5678", ";",
            "_cat.key3 ABCD",
        ], mmcif)
        self.assertEqual(mmcif, {"cat": [
            {"key1": "1234", "key2": "VAL", "key3": "ABCD"}
        ]})



class LoopSectionTests(TestCase):

    def setUp(self):
        self.patch = patch("atomium.mmcif.get_line_values")
        self.mock_values = self.patch.start()


    def tearDown(self):
        self.patch.stop()    


    def test_basic_table(self):
        mmcif = {}
        self.mock_values.side_effect = [["1", "2"], ["3", "4"], ["5", "6"]]
        add_loop_section([
            "loop_", "_db2.id", "_db2.code",
            "PDB   2BFB", "PDBE  EBI-21880", "DB3 XXX"
        ], mmcif)
        self.assertEqual(mmcif, {"db2": [
            {"id": "1", "code": "2"},
            {"id": "3", "code": "4"},
            {"id": "5", "code": "6"}
        ]})
    

    def test_blank_lines_row(self):
        mmcif = {}
        self.mock_values.side_effect = [["1", "2"], ["3", "4"], ["5", "6"]]
        add_loop_section([
            "loop_", "_db2.id", "_db2.code",
            "PDB   2BFB", " ", "PDBE  EBI-21880", "", "DB3 XXX"
        ], mmcif)
        self.assertEqual(mmcif, {"db2": [
            {"id": "1", "code": "2"},
            {"id": "3", "code": "4"},
            {"id": "5", "code": "6"}
        ]})
        

    def test_single_line_quote_values(self):
        mmcif = {}
        self.mock_values.side_effect = [["1", "2", "3"], ["4"], ["5"], ["6"], ["7"], ["8"], ["9"]]
        add_loop_section([
            "loop_", "_cat.id", "_cat.name", "_cat.age",
            "20 Mike 2", "30", "'Jane Smith'", "3", "40", "\"Lee Bee\"", "4",
        ], mmcif)
        self.assertEqual(mmcif, {"cat": [
            {"id": "1", "name": "2", "age": "3"},
            {"id": "4", "name": "5", "age": "6"},
            {"id": "7", "name": "8", "age": "9"},
        ]})
    

    @patch("atomium.mmcif.get_semicolon_value")
    def test_semicolon_values(self, mock_semi):
        mmcif = {}
        mock_semi.side_effect = lambda d: d.popleft() and "VAL"
        self.mock_values.side_effect = [["1", "2", "3"], ["4"], ["6"]]
        add_loop_section([
            "loop_", "_cat.id", "_cat.name", "_cat.age",
            "20 Mike 3",
            "30", ";Roger Mortimer", ";", "9",
        ], mmcif)
        self.assertEqual(mmcif, {"cat": [
            {"id": "1", "name": "2", "age": "3"},
            {"id": "4", "name": "VAL", "age": "6"},
        ]})



class LineValuesTests(TestCase):

    def test_can_get_normal_values(self):
        self.assertEqual(get_line_values("AAA BBB CCC"), ["AAA", "BBB", "CCC"])
    

    def test_can_get_single_quote_values(self):
        self.assertEqual(get_line_values("AAA 'B \" B' CCC"), ["AAA", "B \" B", "CCC"])
    

    def test_can_get_double_quote_values(self):
        self.assertEqual(get_line_values("AAA \"B ' B\" CCC"), ["AAA", "B ' B", "CCC"])
    

    def test_can_handle_multiple_spaces(self):
        self.assertEqual(get_line_values("AAA 'BX    B' CCC"), ["AAA", "BX    B", "CCC"])



class SemicolonValueTests(TestCase):

    def test_can_get_single_line_semicolon_value(self):
        lines = deque([";A single value", ";", "Next line"])
        self.assertEqual(get_semicolon_value(lines), "A single value")
        self.assertEqual(lines, deque([";", "Next line"]))
    

    def test_can_get_multi_line_semicolon_value_with_spaces(self):
        lines = deque(["; A single value", "Another", ";", "Next line"])
        self.assertEqual(get_semicolon_value(lines), "A single value\nAnother")
        self.assertEqual(lines, deque([";", "Next line"]))
    

    def test_can_get_multi_line_semicolon_value_without_spaces(self):
        lines = deque([";ABCDE", "FGEH", ";", "Next line"])
        self.assertEqual(get_semicolon_value(lines), "ABCDEFGEH")
        self.assertEqual(lines, deque([";", "Next line"]))



class MmcifDictToMmcifFilestringTests(TestCase):

    @patch("atomium.mmcif.get_non_loop_lines")
    @patch("atomium.mmcif.get_loop_lines")
    def test_can_save_mmcif_dict(self, mock_loop, mock_non):
        mock_loop.side_effect = [[], ["loop1", "loop2"]]
        mock_non.return_value = ["nonloop1", "nonloop2"]
        mmcif = {"cat1": [], "cat2": [1], "cat3": [1, 2]}
        filestring = mmcif_dict_to_mmcif_filestring(mmcif)
        mock_loop.assert_any_call("cat1", [])
        mock_loop.assert_any_call("cat3", [1, 2])
        mock_non.assert_called_with("cat2", [1])
        output = [
            "data_XXXX",
            "#", "nonloop1", "nonloop2", "#", "loop1", "loop2", "#"
        ]
        self.assertEqual(filestring, "\n".join(output))
    

    @patch("atomium.mmcif.get_non_loop_lines")
    @patch("atomium.mmcif.get_loop_lines")
    def test_can_save_mmcif_dict_with_entry_id(self, mock_loop, mock_non):
        mock_loop.side_effect = [[], ["loop1", "loop2"]]
        mock_non.return_value = ["nonloop1", "nonloop2"]
        mmcif = {"cat1": [], "entry": [{"id": ".."}], "cat3": [1, 2]}
        filestring = mmcif_dict_to_mmcif_filestring(mmcif)
        mock_loop.assert_any_call("cat1", [])
        mock_loop.assert_any_call("cat3", [1, 2])
        mock_non.assert_called_with("entry", [{"id": ".."}])
        output = [
            "data_..",
            "#", "nonloop1", "nonloop2", "#", "loop1", "loop2", "#"
        ]
        self.assertEqual(filestring, "\n".join(output))



class NonLoopLineGenerationTests(TestCase):

    @patch("atomium.mmcif.format_value")
    def test_lines_from_normal_values(self, mock_format):
        mock_format.side_effect = lambda l: l
        lines = get_non_loop_lines("cat1", [{"attr1": "xxx", "mod_attr2": "yyy"}])
        self.assertEqual(lines, [
            "_cat1.attr1     xxx",
            "_cat1.mod_attr2 yyy"
        ])
        mock_format.assert_any_call("xxx")
        mock_format.assert_any_call("yyy")
    

    @patch("atomium.mmcif.format_value")
    def test_lines_from_list_values(self, mock_format):
        mock_format.side_effect = ["xxx", [";line1", "line2", ";"], "zzz"]
        lines = get_non_loop_lines("cat1", [
            {"attr1": "xxx", "mod_attr2": "yyy", "attr_2": "zzz"}
        ])
        self.assertEqual(lines, [
            "_cat1.attr1     xxx",
            "_cat1.mod_attr2",
            ";line1", "line2", ";",
            "_cat1.attr_2    zzz",
        ])
        mock_format.assert_any_call("xxx")
        mock_format.assert_any_call("yyy")
        mock_format.assert_any_call("zzz")



class LoopLineGenerationTests(TestCase):

    def test_can_handle_empty_section(self):
        self.assertEqual(get_loop_lines("cat", []), [])


    @patch("atomium.mmcif.format_value")
    def test_lines_from_normal_values(self, mock_format):
        mock_format.side_effect = lambda l: l
        category = [
            {"attr1": "1111", "attr2": "2", "attr3": "3"},
            {"attr1": "444", "attr2": "555", "attr3": "66"},
        ]
        lines = get_loop_lines("cat1", category)
        self.assertEqual(lines, [
            "loop_", "_cat1.attr1", "_cat1.attr2", "_cat1.attr3",
            "1111 2   3 ", "444  555 66"
        ])
        for val in ["1111", "2", "3", "444", "555", "66"]:
            mock_format.assert_any_call(val)
        
    
    @patch("atomium.mmcif.format_value")
    def test_semicolon_value_at_start(self, mock_format):
        mock_format.side_effect = [[";line1", "line2", ";"], "2", "3", [";line3", "line4", ";"], "555", "66"]
        category = [
            {"attr1": "1111", "attr2": "2", "attr3": "3"},
            {"attr1": "444", "attr2": "555", "attr3": "66"},
        ]
        lines = get_loop_lines("cat1", category)
        self.assertEqual(lines, [
            "loop_", "_cat1.attr1", "_cat1.attr2", "_cat1.attr3",
            ";line1", "line2", ";", "2   3 ",
            ";line3", "line4", ";", "555 66"
        ])
        for val in ["1111", "2", "3", "444", "555", "66"]:
            mock_format.assert_any_call(val)
    

    @patch("atomium.mmcif.format_value")
    def test_semicolon_values_in_middle(self, mock_format):
        mock_format.side_effect = ["1111", [";line1", "line2", ";"], "3", "444", [";line3", "line4", ";"], "66"]
        category = [
            {"attr1": "1111", "attr2": "2", "attr3": "3"},
            {"attr1": "444", "attr2": "555", "attr3": "66"},
        ]
        lines = get_loop_lines("cat1", category)
        self.assertEqual(lines, [
            "loop_", "_cat1.attr1", "_cat1.attr2", "_cat1.attr3",
            "1111", ";line1", "line2", ";", "3 ",
            "444 ", ";line3", "line4", ";", "66"
        ])
        for val in ["1111", "2", "3", "444", "555", "66"]:
            mock_format.assert_any_call(val)
    

    @patch("atomium.mmcif.format_value")
    def test_semicolon_values_at_end(self, mock_format):
        mock_format.side_effect = ["1111", "2", [";line1", "line2", ";"], "444", "555", [";line3", "line4", ";"]]
        category = [
            {"attr1": "1111", "attr2": "2", "attr3": "3"},
            {"attr1": "444", "attr2": "555", "attr3": "66"},
        ]
        lines = get_loop_lines("cat1", category)
        self.assertEqual(lines, [
            "loop_", "_cat1.attr1", "_cat1.attr2", "_cat1.attr3",
            "1111 2  ", ";line1", "line2", ";",
            "444  555", ";line3", "line4", ";",
        ])
        for val in ["1111", "2", "3", "444", "555", "66"]:
            mock_format.assert_any_call(val)
    

    @patch("atomium.mmcif.format_value")
    def test_semicolon_values_at_different_locations(self, mock_format):
        mock_format.side_effect = [[";line1", "line2", ";"], "2", "3", "444", [";line3", "line4", ";"], "66"]
        category = [
            {"attr1": "1111", "attr2": "2", "attr3": "3"},
            {"attr1": "444", "attr2": "555", "attr3": "66"},
        ]
        lines = get_loop_lines("cat1", category)
        self.assertEqual(lines, [
            "loop_", "_cat1.attr1", "_cat1.attr2", "_cat1.attr3",
            ";line1", "line2", ";", "2 3 ",
            "444", ";line3", "line4", ";", "66"
        ])
        for val in ["1111", "2", "3", "444", "555", "66"]:
            mock_format.assert_any_call(val)



class ValueFormattingTests(TestCase):

    def test_can_return_unchanged(self):
        self.assertEqual(format_value("xxx"), "xxx")
        self.assertEqual(format_value("-ABC;"), "-ABC;")
    

    def test_can_escape_quotes(self):
        self.assertEqual(format_value("x\"xx"), "'x\"xx'")
        self.assertEqual(format_value("x'xx"), "\"x'xx\"")
    

    def test_can_handle_spaces(self):
        self.assertEqual(format_value("x x x"), "'x x x'")
        self.assertEqual(format_value("x x\"x"), "'x x\"x'")
        self.assertEqual(format_value("x x'x"), "\"x x'x\"")
    

    def test_can_handle_line_breaks(self):
        self.assertEqual(format_value("xx\nx\nxy"), [";xx", "x", "xy", ";"])
    

    def test_can_handle_multiple_quotes(self):
        self.assertEqual(format_value("x\"x'x"), [";x\"x'x", ";"])