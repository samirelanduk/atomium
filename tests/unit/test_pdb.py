from datetime import date
from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.pdb import *

class PdbStringToPdbDictTests(TestCase):

    @patch("atomium.pdb.update_dict")
    def test_can_convert_pdb_string_to_pdb_dict_one_model(self, mock_up):
        ring = "\n".join([
         "HEAD   LINE ", "TITLE  L1     ", "TITLE  L2", ""
         "REMARK 1  X", "REMARK 2  Y", "REMARK 2  Z",
         "ATOM   1", "ATOM   2", "HETATM 1",
         "CONECT 10", "CONECT 20"
        ])
        mock_up.side_effect = lambda d, k, v: d.update({k: v})
        d = pdb_string_to_pdb_dict(ring)
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
        ring = "\n".join([
         "HEAD   LINE ", "TITLE  L1     ", "TITLE  L2",
         "MODEL    1", "ATOM   1", "HETATM 1", "ENDMDL",
         "MODEL    2", "ATOM   2", "HETATM 2", "ENDMDL",
         "MODEL    3", "ATOM   3", "HETATM 3", "ENDMDL",
         "CONECT 10", "CONECT 20"
        ])
        mock_up.side_effect = lambda d, k, v: d.update({k: v})
        d = pdb_string_to_pdb_dict(ring)
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



class PdbDictToDataDictTests(TestCase):

    @patch("atomium.pdb.update_description_dict")
    @patch("atomium.pdb.update_experiment_dict")
    @patch("atomium.pdb.update_quality_dict")
    def test_can_convert_pdb_dict_to_data_dict(self, mock_ql, mock_ex, mock_ds):
        pdb_dict = {"A": "B"}
        d = pdb_dict_to_data_dict(pdb_dict)
        mock_ds.assert_called_with(pdb_dict, d)
        mock_ex.assert_called_with(pdb_dict, d)
        mock_ql.assert_called_with(pdb_dict, d)
        self.assertEqual(d, {
         "description": {
          "code": None, "title": None, "deposition_date": None,
          "classification": None, "keywords": [], "authors": []
         }, "experiment": {
          "technique": None, "source_organism": None, "expression_system": None
         }, "quality": {"resolution": None, "rvalue": None, "rfree": None}
        })



class DescriptionDictUpdatingTests(TestCase):

    @patch("atomium.pdb.extract_header")
    @patch("atomium.pdb.extract_title")
    @patch("atomium.pdb.extract_keywords")
    @patch("atomium.pdb.extract_authors")
    def test_can_update_description_dict(self, mock_aut, mock_key, mock_tit, mock_hed):
        d = {"description": "dict"}
        pdb_dict = {"PDB": "DICT"}
        update_description_dict(pdb_dict, d)
        mock_hed.assert_called_with(pdb_dict, "dict")
        mock_tit.assert_called_with(pdb_dict, "dict")
        mock_key.assert_called_with(pdb_dict, "dict")
        mock_aut.assert_called_with(pdb_dict, "dict")



class ExperimentDictUpdatingTests(TestCase):

    @patch("atomium.pdb.extract_technique")
    @patch("atomium.pdb.extract_source")
    def test_can_update_experiment_dict(self, mock_src, mock_tech):
        d = {"experiment": "dict"}
        pdb_dict = {"PDB": "DICT"}
        update_experiment_dict(pdb_dict, d)
        mock_src.assert_called_with(pdb_dict, "dict")
        mock_tech.assert_called_with(pdb_dict, "dict")



class QualityDictUpdatingTests(TestCase):

    @patch("atomium.pdb.extract_resolution_remark")
    @patch("atomium.pdb.extract_rvalue_remark")
    def test_can_update_quality_dict(self, mock_rfac, mock_res):
        d = {"quality": "dict"}
        pdb_dict = {"PDB": "DICT"}
        update_quality_dict(pdb_dict, d)
        mock_res.assert_called_with(pdb_dict, "dict")
        mock_rfac.assert_called_with(pdb_dict, "dict")



class HeaderExtractionTests(TestCase):

    def setUp(self):
        self.d = {"code": None, "deposition_date": None, "classification": None}


    def test_can_extract_no_header(self):
        extract_header({"TITLE": ["1"]}, self.d)
        self.assertEqual(self.d, {"code": None, "deposition_date": None, "classification": None})


    def test_can_extract_empty_header(self):
        extract_header({"HEADER": [" " * 74]}, self.d)
        self.assertEqual(self.d, {"code": None, "deposition_date": None, "classification": None})


    def test_can_extract_header(self):
        extract_header({"HEADER": [
         "HEADER    UNKNOWN FUNCTION" + " " * 24 + "21-AUG-17   6AR7" + " " * 14
        ]}, self.d)
        self.assertEqual(self.d, {
         "code": "6AR7", "deposition_date": date(2017, 8, 21), "classification": "UNKNOWN FUNCTION"
        })



class TitleExtractionTests(TestCase):

    def test_missing_title_extraction(self):
        d = {"title": None}
        extract_title({"HEADER": ["1"]}, d)
        self.assertEqual(d, {"title": None})


    @patch("atomium.pdb.merge_lines")
    def test_title_extraction(self, mock_merge):
        d = {"title": None}
        mock_merge.return_value = "TITLE TITLE TITLE"
        extract_title({"TITLE": ["TITLE     L1", "TITLE    2 L2"]}, d)
        mock_merge.assert_called_with(["TITLE     L1", "TITLE    2 L2"], 10)
        self.assertEqual(d["title"], "TITLE TITLE TITLE")



class KeywordExtractionTests(TestCase):

    def test_missing_keyword_extraction(self):
        d = {"keywords": []}
        extract_keywords({"HEADER": ["1"]}, d)
        self.assertEqual(d, {"keywords": []})


    @patch("atomium.pdb.merge_lines")
    def test_keywords_extraction(self, mock_merge):
        d = {"keywords": []}
        mock_merge.return_value = "KEY1, KEY2, KEY3"
        extract_keywords({"KEYWDS": ["KEY     L1", "KEY    2 L2"]}, d)
        mock_merge.assert_called_with(["KEY     L1", "KEY    2 L2"], 10)
        self.assertEqual(d["keywords"], ["KEY1", "KEY2", "KEY3"])



class AuthorExtractionTests(TestCase):

    def test_missing_author_extraction(self):
        d = {"authors": []}
        extract_authors({"HEADER": ["1"]}, d)
        self.assertEqual(d, {"authors": []})


    @patch("atomium.pdb.merge_lines")
    def test_authors_extraction(self, mock_merge):
        d = {"authors": []}
        mock_merge.return_value = "AT1, AT2, AT3"
        extract_authors({"AUTHOR": ["AT     L1", "AT    2 L2"]}, d)
        mock_merge.assert_called_with(["AT     L1", "AT    2 L2"], 10)
        self.assertEqual(d["authors"], ["AT1", "AT2", "AT3"])



class TechniqueExtractionTests(TestCase):

    def test_missing_technique_extraction(self):
        d = {"technique": None}
        extract_technique({"HEADER": ["1"]}, d)
        self.assertEqual(d, {"technique": None})


    def test_empty_technique_extraction(self):
        d = {"technique": None}
        extract_technique({"EXPDTA": ["     "]}, d)
        self.assertEqual(d, {"technique": None})


    def test_technique_extraction(self):
        d = {"technique": None}
        extract_technique({"EXPDTA": ["EXPDTA    X-RAY DIFFRACTION       "]}, d)
        self.assertEqual(d, {"technique": "X-RAY DIFFRACTION"})



class SourceExtractionTests(TestCase):

    def test_missing_source_extraction(self):
        d = {"source_organism": None, "expression_system": None}
        extract_source({"HEADER": ["1"]}, d)
        self.assertEqual(d, {"source_organism": None, "expression_system": None})


    @patch("atomium.pdb.merge_lines")
    def test_empty_source_extraction(self, mock_merge):
        mock_merge.return_value = "JYVGBHUBBGYBHKJNHBK"
        d = {"source_organism": None, "expression_system": None}
        extract_source({"SOURCE": ["1", "2"]}, d)
        self.assertEqual(d, {"source_organism": None, "expression_system": None})
        mock_merge.assert_called_with(["1", "2"], 10)


    @patch("atomium.pdb.merge_lines")
    def test_empty_source_extraction(self, mock_merge):
        mock_merge.return_value = (
         "MOL_ID: 1;"
         " ORGANISM_SCIENTIFIC: METHANOTHERMOBACTER"
         " THERMAUTOTROPHICUS STR. DELTA H;"
         " ORGANISM_TAXID: 187420;"
         " STRAIN: DELTA H;"
         " EXPRESSION_SYSTEM: ESCHERICHIA COLI;"
         " EXPRESSION_SYSTEM_TAXID: 562;"
         " EXPRESSION_SYSTEM_PLASMID: PET15B"
        )
        d = {"source_organism": None, "expression_system": None}
        extract_source({"SOURCE": ["1", "2"]}, d)
        self.assertEqual(d, {
         "source_organism": "METHANOTHERMOBACTER THERMAUTOTROPHICUS STR. DELTA H",
         "expression_system": "ESCHERICHIA COLI"
        })
        mock_merge.assert_called_with(["1", "2"], 10)



class ResolutionExtractionTests(TestCase):

    def setUp(self):
        self.remark_lines = {
         "1": ["", "BLAH BLAH"],
         "2": ["", "REMARK 2   RESOLUTION.    1.90 ANGSTROMS."],
         "24": ["", "BLAH BLAH"]
        }


    def test_missing_remarks_extraction(self):
        d = {"resolution": None}
        del self.remark_lines["2"]
        extract_resolution_remark({"REMARK": self.remark_lines}, d)
        self.assertEqual(d, {"resolution": None})
        extract_resolution_remark({"ABC": []}, d)
        self.assertEqual(d, {"resolution": None})


    def test_empty_resolution_extraction(self):
        d = {"resolution": None}
        self.remark_lines["2"][1] = "REMARK 2   RESOLUTION. NOT APPLICABLE."
        extract_resolution_remark({"REMARK": self.remark_lines}, d)
        self.assertEqual(d, {"resolution": None})


    def test_resolution_extraction(self):
        d = {"resolution": None}
        extract_resolution_remark({"REMARK": self.remark_lines}, d)
        self.assertEqual(d, {"resolution": 1.9})



class RvalueExtractionTests(TestCase):

    def setUp(self):
        self.remark_lines = {
         "1": ["", "BLAH BLAH"],
         "3": [
          "REMARK 3     CROSS-VALIDATION METHOD          : THROUGHOUT",
          "REMARK 3     FREE R VALUE TEST SET SELECTION  : RANDOM",
          "REMARK 3     R VALUE            (WORKING SET) : 0.193",
          "REMARK 3     FREE R VALUE                     : 0.229",
          "REMARK 3     FREE R VALUE TEST SET SIZE   (%) : 4.900",
          "REMARK 3     FREE R VALUE TEST SET COUNT      : 1583",
          "REMARK 3     BIN R VALUE           (WORKING SET) : 0.2340 "
         ],
         "24": ["", "BLAH BLAH"]
        }


    def test_missing_rvalue_extraction(self):
        d = {"rvalue": None, "rfree": None}
        del self.remark_lines["3"]
        extract_resolution_remark({"REMARK": self.remark_lines}, d)
        self.assertEqual(d, {"rvalue": None, "rfree": None})
        extract_rvalue_remark({"ABC": []}, d)
        self.assertEqual(d, {"rvalue": None, "rfree": None})


    def test_empty_rvalue_extraction(self):
        d = {"rvalue": None, "rfree": None}
        self.remark_lines["3"][2] = self.remark_lines["3"][2][:-7]
        self.remark_lines["3"][3] = self.remark_lines["3"][3][:-7]
        self.remark_lines["3"].pop()
        extract_rvalue_remark({"REMARK": self.remark_lines}, d)
        self.assertEqual(d, {"rvalue": None, "rfree": None})


    def test_rvalue_extraction(self):
        d = {"rvalue": None, "rfree": None}
        extract_rvalue_remark({"REMARK": self.remark_lines}, d)
        self.assertEqual(d, {"rvalue": 0.193, "rfree": 0.229})



class LineMergingTests(TestCase):

    def setUp(self):
        self.lines = ["0123456789 ", "abcdefghij ", "0123456789 "]
        self.punc_lines = ["0123, 456789 ", "abcd  efghij ", "0123; 456789 "]


    def test_can_merge_lines(self):
        self.assertEqual(
         merge_lines(self.lines, 5),
         "56789 fghij 56789"
        )
        self.assertEqual(
         merge_lines(self.lines, 8),
         "89 ij 89"
        )


    def test_can_vary_join(self):
        self.assertEqual(
         merge_lines(self.lines, 5, join=""),
         "56789fghij56789"
        )
        self.assertEqual(
         merge_lines(self.lines, 8, join="."),
         "89.ij.89"
        )
