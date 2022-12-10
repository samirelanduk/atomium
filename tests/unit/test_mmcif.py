from unittest import TestCase
from unittest.mock import patch
from atomium.mmcif import *

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
    

    def test_spaces_in_values(self):
        mmcif = {}
        add_non_loop_section([
            "_refine.pdbx_refine_id 'X-RAY DIFFRACTION'",
            "_refine.name       \"XANTHOSINE-5'-MONOPHOSPHATE\""
        ], mmcif)
        self.assertEqual(mmcif["refine"], [
            {"pdbx_refine_id": "X-RAY DIFFRACTION", "name": "XANTHOSINE-5'-MONOPHOSPHATE"}
        ])
    

    def test_quote_strings(self):
        mmcif = {}
        add_non_loop_section(["_atoms.name1 \"N1'\"", "_atoms.name2 'N1\"'"], mmcif)
        self.assertEqual(mmcif["atoms"], [{"name1": "N1'", "name2": "N1\""}])
    

    def test_value_on_next_line(self):
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