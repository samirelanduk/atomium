from datetime import datetime
from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.files.pdbstring2pdbdict import *

class PdbStringConversionTest(TestCase):

    def setUp(self):
        self.pdb_dict = {}
        self.lines = ["line1", "line2", "line3"]



class PdbStringToPdbDictTests(TestCase):

    @patch("atomium.files.utilities.string_to_lines")
    @patch("atomium.files.pdbstring2pdbdict.extract_annotation")
    @patch("atomium.files.pdbstring2pdbdict.extract_structure")
    def test_can_convert_pdb_string_to_dict(self, mock_struc, mock_ann, mock_lines):
        mock_lines.return_value = ["line1", "line2"]
        pdb_dict = pdb_string_to_pdb_dict("filestring")
        mock_lines.assert_called_with("filestring", width=80)
        mock_ann.assert_called_with({}, ["line1", "line2"])
        mock_struc.assert_called_with({}, ["line1", "line2"])
        self.assertEqual(pdb_dict, {})



class AnnotationExtractionTests(PdbStringConversionTest):
    @patch("atomium.files.pdbstring2pdbdict.extract_header")
    @patch("atomium.files.pdbstring2pdbdict.extract_title")
    @patch("atomium.files.pdbstring2pdbdict.extract_resolution")
    @patch("atomium.files.pdbstring2pdbdict.extract_rfactor")
    @patch("atomium.files.pdbstring2pdbdict.extract_source")
    @patch("atomium.files.pdbstring2pdbdict.extract_technique")
    @patch("atomium.files.pdbstring2pdbdict.extract_keywords")
    @patch("atomium.files.pdbstring2pdbdict.extract_sequence")
    def test_can_extract_header(self, mock_seq, mock_key, mock_tech, mock_source, mock_rfac, mock_res, mock_title, mock_header):
        extract_annotation(self.pdb_dict, self.lines)
        mock_title.assert_called_with(self.pdb_dict, self.lines)
        mock_header.assert_called_with(self.pdb_dict, self.lines)
        mock_res.assert_called_with(self.pdb_dict, self.lines)
        mock_rfac.assert_called_with(self.pdb_dict, self.lines)
        mock_source.assert_called_with(self.pdb_dict, self.lines)
        mock_tech.assert_called_with(self.pdb_dict, self.lines)
        mock_key.assert_called_with(self.pdb_dict, self.lines)
        mock_seq.assert_called_with(self.pdb_dict, self.lines)



class HeaderExtractionTests(PdbStringConversionTest):

    @patch("atomium.files.pdbstring2pdbdict.get_line")
    def test_empty_header_extraction(self, mock_line):
        mock_line.return_value = "HEADER".ljust(80)
        extract_header(self.pdb_dict, self.lines)
        mock_line.assert_called_with("HEADER", self.lines)
        self.assertEqual(self.pdb_dict["deposition_date"], None)
        self.assertEqual(self.pdb_dict["code"], None)
        self.assertEqual(self.pdb_dict["classification"], None)


    @patch("atomium.files.pdbstring2pdbdict.get_line")
    def test_missing_header_extraction(self, mock_line):
        mock_line.return_value = None
        extract_header(self.pdb_dict, self.lines)
        mock_line.assert_called_with("HEADER", self.lines)
        self.assertEqual(self.pdb_dict["deposition_date"], None)
        self.assertEqual(self.pdb_dict["code"], None)
        self.assertEqual(self.pdb_dict["classification"], None)


    @patch("atomium.files.pdbstring2pdbdict.get_line")
    def test_header_extraction(self, mock_line):
        mock_line.return_value = (
         "HEADER    UNKNOWN FUNCTION" + " " * 24 + "21-AUG-17   6AR7" + " " * 14
        )
        extract_header(self.pdb_dict, self.lines)
        mock_line.assert_called_with("HEADER", self.lines)
        self.assertEqual(
         self.pdb_dict["deposition_date"], datetime(2017, 8, 21).date()
        )
        self.assertEqual(self.pdb_dict["code"], "6AR7")
        self.assertEqual(self.pdb_dict["classification"], "UNKNOWN FUNCTION")



class TitleExtractionTests(PdbStringConversionTest):

    @patch("atomium.files.pdbstring2pdbdict.get_lines")
    def test_missing_title_extraction(self, mock_lines):
        mock_lines.return_value = []
        extract_title(self.pdb_dict, self.lines)
        mock_lines.assert_any_call("TITLE", self.lines)
        self.assertEqual(self.pdb_dict["title"], None)


    @patch("atomium.files.pdbstring2pdbdict.get_lines")
    @patch("atomium.files.pdbstring2pdbdict.merge_lines")
    def test_title_extraction(self, mock_merge, mock_lines):
        mock_lines.return_value = [
         "TITLE     L1".ljust(80), "TITLE    2 L2".ljust(80)
        ]
        mock_merge.return_value = "TITLE TITLE TITLE"
        extract_title(self.pdb_dict, self.lines)
        mock_lines.assert_any_call("TITLE", self.lines)
        mock_merge.assert_called_with(mock_lines.return_value, 10)
        self.assertEqual(self.pdb_dict["title"], "TITLE TITLE TITLE")



class ResolutionExtractionTests(PdbStringConversionTest):

    def setUp(self):
        PdbStringConversionTest.setUp(self)
        self.remark_lines = [
         "REMARK   1",
         "REMARK   1 BLAH BLAH.",
         "REMARK   2",
         "REMARK   2 RESOLUTION.    1.90 ANGSTROMS.",
         "REMARK  24",
         "REMARK  24 BLAH BLAH."
        ]

    @patch("atomium.files.pdbstring2pdbdict.get_lines")
    def test_empty_resolution_extraction(self, mock_lines):
        self.remark_lines[1] = "REMARK   2 RESOLUTION. NOT APPLICABLE."
        mock_lines.return_value = self.remark_lines
        extract_resolution(self.pdb_dict, self.lines)
        mock_lines.assert_any_call("REMARK", self.lines)
        self.assertEqual(self.pdb_dict["resolution"], None)


    @patch("atomium.files.pdbstring2pdbdict.get_lines")
    def test_missing_remarks_extraction(self, mock_lines):
        mock_lines.return_value = []
        extract_resolution(self.pdb_dict, self.lines)
        mock_lines.assert_any_call("REMARK", self.lines)
        self.assertEqual(self.pdb_dict["resolution"], None)


    @patch("atomium.files.pdbstring2pdbdict.get_lines")
    def test_resolution_extraction(self, mock_lines):
        mock_lines.return_value = self.remark_lines
        extract_resolution(self.pdb_dict, self.lines)
        mock_lines.assert_any_call("REMARK", self.lines)
        self.assertEqual(self.pdb_dict["resolution"], 1.9)



class RfactorExtractionTests(PdbStringConversionTest):

    def setUp(self):
        PdbStringConversionTest.setUp(self)
        self.remark_lines = [
         "REMARK   1",
         "REMARK   1 BLAH BLAH.",
         "REMARK   3   CROSS-VALIDATION METHOD          : THROUGHOUT",
         "REMARK   3   FREE R VALUE TEST SET SELECTION  : RANDOM",
         "REMARK   3   R VALUE            (WORKING SET) : 0.193",
         "REMARK   3   FREE R VALUE                     : 0.229",
         "REMARK   3   FREE R VALUE TEST SET SIZE   (%) : 4.900",
         "REMARK   3   FREE R VALUE TEST SET COUNT      : 1583",
         "REMARK  24",
         "REMARK  24 BLAH BLAH."
        ]


    @patch("atomium.files.pdbstring2pdbdict.get_lines")
    def test_empty_rfactor_extraction(self, mock_lines):
        self.remark_lines[4] = "REMARK   3   R VALUE             : 0.193"
        mock_lines.return_value = self.remark_lines
        extract_rfactor(self.pdb_dict, self.lines)
        mock_lines.assert_any_call("REMARK", self.lines)
        self.assertEqual(self.pdb_dict["rfactor"], None)


    @patch("atomium.files.pdbstring2pdbdict.get_lines")
    def test_missing_remarks_extraction(self, mock_lines):
        mock_lines.return_value = []
        extract_rfactor(self.pdb_dict, self.lines)
        mock_lines.assert_any_call("REMARK", self.lines)
        self.assertEqual(self.pdb_dict["rfactor"], None)


    @patch("atomium.files.pdbstring2pdbdict.get_lines")
    def test_rfactor_extraction(self, mock_lines):
        mock_lines.return_value = self.remark_lines
        extract_rfactor(self.pdb_dict, self.lines)
        mock_lines.assert_any_call("REMARK", self.lines)
        self.assertEqual(self.pdb_dict["rfactor"], 0.193)
        self.assertEqual(self.pdb_dict["rfree"], 0.229)
        self.assertEqual(self.pdb_dict["rcount"], 1583)



class SourceExtractionTests(PdbStringConversionTest):

    @patch("atomium.files.pdbstring2pdbdict.get_lines")
    def test_no_source_extraction(self, mock_lines):
        mock_lines.return_value = []
        extract_source(self.pdb_dict, self.lines)
        mock_lines.assert_any_call("SOURCE", self.lines)
        self.assertEqual(self.pdb_dict["organism"], None)
        self.assertEqual(self.pdb_dict["expression_system"], None)


    @patch("atomium.files.pdbstring2pdbdict.get_lines")
    def test_missing_source_extraction(self, mock_lines):
        mock_lines.return_value = self.lines
        extract_source(self.pdb_dict, self.lines)
        mock_lines.assert_any_call("SOURCE", self.lines)
        self.assertEqual(self.pdb_dict["organism"], None)
        self.assertEqual(self.pdb_dict["expression_system"], None)


    @patch("atomium.files.pdbstring2pdbdict.get_lines")
    @patch("atomium.files.pdbstring2pdbdict.merge_lines")
    def test_organism_source_extraction(self, mock_merge,  mock_lines):
        mock_lines.return_value = self.lines
        mock_merge.return_value = (
         "MOL_ID: 1;"
         " ORGANISM_SCIENTIFIC: METHANOTHERMOBACTER"
         " THERMAUTOTROPHICUS STR. DELTA H;"
         " ORGANISM_TAXID: 187420;"
         " STRAIN: DELTA H;"
        )
        extract_source(self.pdb_dict, self.lines)
        mock_lines.assert_any_call("SOURCE", self.lines)
        mock_merge.assert_called_with(self.lines, 10)
        self.assertEqual(
         self.pdb_dict["organism"],
         "METHANOTHERMOBACTER THERMAUTOTROPHICUS STR. DELTA H"
        )
        self.assertEqual(self.pdb_dict["expression_system"], None)


    @patch("atomium.files.pdbstring2pdbdict.get_lines")
    @patch("atomium.files.pdbstring2pdbdict.merge_lines")
    def test_expression_source_extraction(self, mock_merge,  mock_lines):
        mock_lines.return_value = self.lines
        mock_merge.return_value = (
         " STRAIN: DELTA H;"
         " EXPRESSION_SYSTEM: ESCHERICHIA COLI;"
         " EXPRESSION_SYSTEM_TAXID: 562;"
         " EXPRESSION_SYSTEM_PLASMID: PET15B"
        )
        extract_source(self.pdb_dict, self.lines)
        mock_lines.assert_any_call("SOURCE", self.lines)
        mock_merge.assert_called_with(self.lines, 10)
        self.assertEqual(
         self.pdb_dict["expression_system"],
         "ESCHERICHIA COLI"
        )
        self.assertEqual(self.pdb_dict["organism"], None)



class TechniqueExtractionTests(PdbStringConversionTest):

    @patch("atomium.files.pdbstring2pdbdict.get_line")
    def test_empty_technique_extraction(self, mock_line):
        mock_line.return_value = "EXPDTA".ljust(80)
        extract_technique(self.pdb_dict, self.lines)
        mock_line.assert_called_with("EXPDTA", self.lines)
        self.assertEqual(self.pdb_dict["technique"], None)


    @patch("atomium.files.pdbstring2pdbdict.get_line")
    def test_missing_tecnhique_extraction(self, mock_line):
        mock_line.return_value = None
        extract_technique(self.pdb_dict, self.lines)
        mock_line.assert_called_with("EXPDTA", self.lines)
        self.assertEqual(self.pdb_dict["technique"], None)


    @patch("atomium.files.pdbstring2pdbdict.get_line")
    def test_technique_extraction(self, mock_line):
        mock_line.return_value = "EXPDTA    X-RAY DIFFRACTION".ljust(80)
        extract_technique(self.pdb_dict, self.lines)
        mock_line.assert_called_with("EXPDTA", self.lines)
        self.assertEqual(self.pdb_dict["technique"], "X-RAY DIFFRACTION")



class KeywordExtractionTests(PdbStringConversionTest):

    @patch("atomium.files.pdbstring2pdbdict.get_lines")
    @patch("atomium.files.pdbstring2pdbdict.merge_lines")
    def test_empty_keyword_extraction(self, mock_merge, mock_lines):
        mock_lines.return_value = []
        extract_keywords(self.pdb_dict, self.lines)
        mock_lines.assert_any_call("KEYWDS", self.lines)
        self.assertEqual(self.pdb_dict["keywords"], [])


    @patch("atomium.files.pdbstring2pdbdict.get_lines")
    @patch("atomium.files.pdbstring2pdbdict.merge_lines")
    def test_one_line_keyword_extraction(self, mock_merge, mock_lines):
        mock_lines.return_value = ["KEYWDS    TIM BARREL, LYASE"]
        mock_merge.return_value = "TIM BARREL, LYASE"
        extract_keywords(self.pdb_dict, self.lines)
        mock_lines.assert_any_call("KEYWDS", self.lines)
        mock_merge.assert_any_call(["KEYWDS    TIM BARREL, LYASE"], 10)
        self.assertEqual(self.pdb_dict["keywords"], ["TIM BARREL", "LYASE"])



class SequenceExtractionTests(PdbStringConversionTest):

    @patch("atomium.files.pdbstring2pdbdict.get_lines")
    @patch("atomium.files.pdbstring2pdbdict.merge_lines")
    def test_empty_sequence_extraction(self, mock_merge, mock_lines):
        mock_lines.return_value = []
        extract_sequence(self.pdb_dict, self.lines)
        mock_lines.assert_any_call("SEQRES", self.lines)
        self.assertEqual(self.pdb_dict["sequences"], {})


    @patch("atomium.files.pdbstring2pdbdict.get_lines")
    @patch("atomium.files.pdbstring2pdbdict.merge_lines")
    def test_sequence_extraction(self, mock_merge, mock_lines):
        mock_lines.return_value = ["line1", "line2"]
        mock_merge.return_value = "C 15 JHY FRT C 15 FDR D YTR"
        extract_sequence(self.pdb_dict, self.lines)
        mock_lines.assert_any_call("SEQRES", self.lines)
        mock_merge.assert_called_with(["line1", "line2"], 10)
        self.assertEqual(self.pdb_dict["sequences"], {"C": ["JHY", "FRT", "FDR"], "D": ["YTR"]})



        mock_merge.return_value = (
         "MOL_ID: 1;"
         " ORGANISM_SCIENTIFIC: METHANOTHERMOBACTER"
         " THERMAUTOTROPHICUS STR. DELTA H;"
         " ORGANISM_TAXID: 187420;"
         " STRAIN: DELTA H;"
        )



class StructureExtractionTests(TestCase):

    def setUp(self):
        self.pdb_dict = {}
        self.lines = [
         "model1", "atom1", "an1", "hetam1", "model2", "atom2", "an2", "hetam2", "con1", "con2"
        ]


    @patch("atomium.files.pdbstring2pdbdict.get_lines")
    @patch("atomium.files.pdbstring2pdbdict.extract_connections")
    @patch("atomium.files.pdbstring2pdbdict.lines_to_model")
    def test_can_extract_structure_one_model(self, mock_model, mock_con, mock_lines):
        mock_lines.side_effect = [
         [], ["atom1", "atom2"], ["hetatm1", "hetatm2"], ["an1", "an2"], ["con1", "con2"]
        ]
        mock_model.return_value = {"model": "1"}
        extract_structure(self.pdb_dict, self.lines)
        mock_lines.assert_any_call("MODEL", self.lines, number=True)
        mock_lines.assert_any_call("ATOM", self.lines, number=False)
        mock_lines.assert_any_call("HETATM", self.lines, number=False)
        mock_lines.assert_any_call("ANISOU", self.lines, number=False)
        mock_lines.assert_any_call("CONECT", self.lines)
        mock_model.assert_called_with(["atom1", "atom2"], ["hetatm1", "hetatm2"], ["an1", "an2"])
        mock_con.assert_called_with(self.pdb_dict, ["con1", "con2"])
        self.assertEqual(self.pdb_dict["models"], [{"model": "1"}])


    @patch("atomium.files.pdbstring2pdbdict.get_lines")
    @patch("atomium.files.pdbstring2pdbdict.extract_connections")
    @patch("atomium.files.pdbstring2pdbdict.lines_to_model")
    def test_can_extract_structure_multiple_models(self, mock_model, mock_con, mock_lines):
        mock_lines.side_effect = [
         [[self.lines[0], 0], [self.lines[4], 4]], [[self.lines[1], 1], [self.lines[5], 5]],
         [[self.lines[3], 3], [self.lines[7], 7]], [[self.lines[2], 2], [self.lines[6], 6]],
         [self.lines[8], self.lines[9]],
        ]
        mock_model.side_effect = [{"model": "1"}, {"model": "2"}]
        extract_structure(self.pdb_dict, self.lines)
        mock_lines.assert_any_call("MODEL", self.lines, number=True)
        mock_lines.assert_any_call("ATOM", self.lines, number=True)
        mock_lines.assert_any_call("HETATM", self.lines, number=True)
        mock_lines.assert_any_call("CONECT", self.lines)
        mock_model.assert_any_call([self.lines[1]], [self.lines[3]], [self.lines[2]])
        mock_model.assert_any_call([self.lines[5]], [self.lines[7]], [self.lines[6]])
        mock_con.assert_called_with(self.pdb_dict, self.lines[8:10])
        self.assertEqual(self.pdb_dict["models"], [{"model": "1"}, {"model": "2"}])



class LinesToModelTests(TestCase):

    @patch("atomium.files.pdbstring2pdbdict.atom_line_to_atom_dict")
    @patch("atomium.files.pdbstring2pdbdict.assign_anisou")
    @patch("atomium.files.pdbstring2pdbdict.atoms_to_chains")
    def test_can_convert_lines_to_model(self, mock_chain, mock_an, mock_dict):
        mock_dict.side_effect = [
         {"a": 1}, {"a": 2}, {"a": 3}, {"a": 4},
         {"h": 1}, {"h": 2}, {"h": 3}, {"h": 4}
        ]
        mock_chain.return_value = [{"c": 1}, {"c": 2}]
        model = lines_to_model(
         ["a1", "a2", "a3", "a4"], ["h1", "h2", "h3", "h4"], ["an1", "an2"]
        )
        for char in ["a", "h"]:
            for num in ["1", "2", "3", "4"]:
                mock_dict.assert_any_call(char + num)
        mock_an.assert_called_with(
         ["an1", "an2"],
         [{"a": 1}, {"a": 2}, {"a": 3}, {"a": 4}],
         [{"h": 1}, {"h": 2}, {"h": 3}, {"h": 4}]
        )
        mock_chain.assert_called_with(
         [{"a": 1}, {"a": 2}, {"a": 3}, {"a": 4}],
         [{"h": 1}, {"h": 2}, {"h": 3}, {"h": 4}]
        )
        self.assertEqual(model, {"chains": [{"c": 1}, {"c": 2}]})



class AtomLineToAtomDictTests(TestCase):

    def test_can_convert_empty_line_to_atom(self):
        atom = atom_line_to_atom_dict("ATOM".ljust(80))
        self.assertEqual(atom, {
         "atom_id": None, "atom_name": None,
         "alt_loc": None, "residue_name": None,
         "chain_id": "", "residue_id": None, "insert_code": "", "full_id": "",
         "x": None, "y": None, "z": None,
         "occupancy": 1, "temp_factor": None,
         "element": None, "charge": 0
        })


    def test_can_convert_full_line_to_atom(self):
        atom = atom_line_to_atom_dict(
         "ATOM    107  N1 AGLY B  13C     " +
         "12.681  37.302 -25.211 0.70  15.56           N2-"
        )
        self.assertEqual(atom, {
         "atom_id": 107, "atom_name": "N1",
         "alt_loc": "A", "residue_name": "GLY",
         "chain_id": "B", "residue_id": 13, "insert_code": "C", "full_id": "B13C",
         "x": 12.681, "y": 37.302, "z": -25.211,
         "occupancy": 0.7, "temp_factor": 15.56,
         "element": "N", "charge": -2
        })


    def test_can_handle_numeric_residue_names(self):
        atom = atom_line_to_atom_dict(
         "ATOM    107  N1 A062 B  13C     " +
         "12.681  37.302 -25.211 0.70  15.56           N2-"
        )
        self.assertEqual(atom["residue_name"], "062")



class AtomsToResiduesTests(TestCase):

    def test_can_convert_atoms_to_residues(self):
        atoms = [
         {"full_id": "B10", "residue_name": "GLY"},
         {"full_id": "B10", "residue_name": "GLY"},
         {"full_id": "A10", "residue_name": "VAL"},
         {"full_id": "A10", "residue_name": "VAL"},
         {"full_id": "A10", "residue_name": "VAL"},
        ]
        residues = atoms_to_residues(atoms)
        self.assertEqual(residues, [{
         "id": "B10", "name": "GLY", "atoms": atoms[:2]
        }, {
         "id": "A10", "name": "VAL", "atoms": atoms[2:]
        }])



class AtomsToChainsTests(TestCase):

    @patch("atomium.files.pdbstring2pdbdict.atoms_to_residues")
    def test_can_convert_atoms_to_chains(self, mock_res):
        atoms = [
         {"full_id": "A10", "residue_name": "GLY", "chain_id": "A"},
         {"full_id": "A11", "residue_name": "PRO", "chain_id": "A"},
         {"full_id": "B10", "residue_name": "VAL", "chain_id": "B"},
         {"full_id": "B10", "residue_name": "VAL", "chain_id": "B"},
         {"full_id": "B11", "residue_name": "LYS", "chain_id": "B"},
        ]
        heteroatoms = [
         {"full_id": "A1000", "residue_name": "AAA", "chain_id": "A"},
         {"full_id": "A1100", "residue_name": "AAA", "chain_id": "A"},
         {"full_id": "B1000", "residue_name": "XXX", "chain_id": "B"},
         {"full_id": "B1100", "residue_name": "XXX", "chain_id": "B"},
        ]
        mock_res.side_effect = [[{"r": 1}], [{"r": 2}], [{"r": 3}], [{"r": 4}]]
        chains = atoms_to_chains(atoms, heteroatoms)
        mock_res.assert_any_call(atoms[:2])
        mock_res.assert_any_call(atoms[2:])
        self.assertEqual(chains, [{
         "chain_id": "A", "residues": [{"r": 1}], "ligands": [{"r": 2}]
        }, {
         "chain_id": "B", "residues": [{"r": 3}], "ligands": [{"r": 4}]
        }])



class AnisouAssingingTests(TestCase):

    def test_can_assign_anisou(self):
        atoms = [{"atom_id": 107}, {"atom_id": 108}]
        heteroatoms = [{"atom_id": 109}, {"atom_id": 110}]
        assign_anisou([
         "ANISOU  107  N   GLY A  13     2406   1892   1614    198    519   -328",
         "ANISOU  110  O   GLY A  13     3837   2505   1611    164   -121    189"
        ], atoms, heteroatoms)
        self.assertEqual(atoms, [
         {"atom_id": 107, "anisotropy": [0.2406, 0.1892, 0.1614, 0.0198, 0.0519, -0.0328]},
         {"atom_id": 108, "anisotropy": [0, 0, 0, 0, 0, 0]}
        ])
        self.assertEqual(heteroatoms, [
         {"atom_id": 109, "anisotropy": [0, 0, 0, 0, 0, 0]},
         {"atom_id": 110, "anisotropy": [0.3837, 0.2505, 0.1611, 0.0164, -0.0121, 0.0189]}
        ])



class ConnectionExtractionTests(TestCase):

    def test_extract_connections(self):
        pdb_dict = {}
        conect_lines = [
         "CONECT 1179  746 1184 1195 1203".ljust(80),
         "CONECT 1179 1211 1222".ljust(80),
         "CONECT 1221  544 1017 1020 1022".ljust(80)
        ]
        extract_connections(pdb_dict, conect_lines)
        self.assertEqual(pdb_dict, {"connections": [
         {"atom": 1179, "bond_to": [746, 1184, 1195, 1203, 1211, 1222]},
         {"atom": 1221, "bond_to": [544, 1017, 1020, 1022]}
        ]})




class GetLineTests(TestCase):

    def setUp(self):
        self.lines = ["AAA   X", "AAA   Y", "BBBBBBX"]


    def test_can_get_line(self):
        self.assertEqual(get_line("BBBBBB", self.lines), "BBBBBBX")


    def test_can_get_stripped_line(self):
        self.assertEqual(get_line("AAA", self.lines), "AAA   X")


    def test_can_get_line_with_number(self):
        self.assertEqual(
         get_line("BBBBBB", self.lines, number=True), ["BBBBBBX", 2]
        )


    def test_can_get_none(self):
        self.assertIsNone(get_line("AA", self.lines))



class GetLinesTests(TestCase):

    def setUp(self):
        self.lines = ["AAA   X", "AAA   Y", "BBBBBBX"]


    def test_can_get_lines(self):
        self.assertEqual(get_lines("BBBBBB", self.lines), ["BBBBBBX"])


    def test_can_get_stripped_lines(self):
        self.assertEqual(get_lines("AAA", self.lines), ["AAA   X", "AAA   Y"])


    def test_can_get_lines_with_numbers(self):
        self.assertEqual(
         get_lines("AAA", self.lines, number=True),
         [["AAA   X", 0], ["AAA   Y", 1]]
        )


    def test_can_get_no_lines(self):
        self.assertEqual(get_lines("AA", self.lines), [])



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
