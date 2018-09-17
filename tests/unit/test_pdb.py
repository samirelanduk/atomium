from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.pdb import *

class PdbStringToPdbDictTests(TestCase):

    @patch("atomium.pdb.update_dict")
    def test_can_convert_pdb_string_to_pdb_dict_one_model(self, mock_up):
        filestring = "\n".join([
         "HEAD   LINE ", "TITLE  L1     ", "TITLE  L2", ""
         "REMARK 1  X", "REMARK 2  Y", "REMARK 2  Z",
         "ATOM   1", "ATOM   2", "HETATM 1",
         "CONECT 10", "CONECT 20"
        ])
        mock_up.side_effect = lambda d, k, v: d.update({k: v})
        d = pdb_string_to_pdb_dict(filestring)
        self.assertEqual(d, {
         "HEAD": "HEAD   LINE", "TITLE": "TITLE  L2",
         "REMARK": {"1": "REMARK 1  X", "2": "REMARK 2  Z"},
         "MODEL": [["ATOM   1", "ATOM   2", "HETATM 1"]],
         "CONECT": "CONECT 20"
        })
        mock_up.assert_any_call(d, "HEAD", "HEAD   LINE")
        mock_up.assert_any_call(d, "TITLE", "TITLE  L1")
        mock_up.assert_any_call(d, "TITLE", "TITLE  L2")
        mock_up.assert_any_call(d["REMARK"], "1", "REMARK 1  X")
        mock_up.assert_any_call(d["REMARK"], "2", "REMARK 2  Y")
        mock_up.assert_any_call(d["REMARK"], "2", "REMARK 2  Z")
        mock_up.assert_any_call(d, "CONECT", "CONECT 10")
        mock_up.assert_any_call(d, "CONECT", "CONECT 20")


    @patch("atomium.pdb.update_dict")
    def test_can_convert_pdb_string_to_pdb_dict_multi_model(self, mock_up):
        filestring = "\n".join([
         "HEAD   LINE ", "TITLE  L1     ", "TITLE  L2",
         "MODEL    1", "ATOM   1", "HETATM 1", "ENDMDL",
         "MODEL    2", "ATOM   2", "HETATM 2", "ENDMDL",
         "MODEL    3", "ATOM   3", "HETATM 3", "ENDMDL",
         "CONECT 10", "CONECT 20"
        ])
        mock_up.side_effect = lambda d, k, v: d.update({k: v})
        d = pdb_string_to_pdb_dict(filestring)
        self.assertEqual(d, {
         "HEAD": "HEAD   LINE", "TITLE": "TITLE  L2",
         "MODEL": [["ATOM   1", "HETATM 1"], ["ATOM   2", "HETATM 2"], ["ATOM   3", "HETATM 3"]],
         "CONECT": "CONECT 20"
        })
        mock_up.assert_any_call(d, "HEAD", "HEAD   LINE")
        mock_up.assert_any_call(d, "TITLE", "TITLE  L1")
        mock_up.assert_any_call(d, "TITLE", "TITLE  L2")
        mock_up.assert_any_call(d, "CONECT", "CONECT 10")
        mock_up.assert_any_call(d, "CONECT", "CONECT 20")



class DictUpdatingTests(TestCase):

    def test_can_add_to_list(self):
        d = {"a": [1], "b": 2}
        update_dict(d, "a", 5)
        self.assertEqual(d, {"a": [1, 5], "b": 2})


    def test_can_create_list(self):
        d = {"a": [1], "b": 2}
        update_dict(d, "c", 5)
        self.assertEqual(d, {"a": [1], "b": 2, "c": [5]})
