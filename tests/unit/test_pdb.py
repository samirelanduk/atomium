from unittest import TestCase
from unittest.mock import patch
from atomium.pdb import *

class PdbStringToDictTests(TestCase):

    @patch("atomium.pdb.parse_metadata")
    @patch("atomium.pdb.get_entities")
    @patch("atomium.pdb.build_structure_categories")
    @patch("atomium.pdb.parse_annotation")
    @patch("atomium.pdb.update_auth_and_label")
    def test_can_get_dict(self, mock_auth, mock_ann, mock_build, mock_get, mock_parse):
        mock_parse.return_value = {"mmcif": 1}
        mock_get.return_value = ["polymers", "nonpolymers"]
        d = pdb_string_to_mmcif_dict("filestring\n\n")
        self.assertEqual(d, {"mmcif": 1})
        mock_parse.assert_called_with("filestring\n\n")
        mock_get.assert_called_with("filestring\n\n")
        mock_build.assert_called_with("filestring\n\n", "polymers", "nonpolymers", {"mmcif": 1})
        mock_ann.assert_called_with("filestring\n\n", {"mmcif": 1})
        mock_auth.assert_called_with({"mmcif": 1})



class ParseMetadataTests(TestCase):

    @patch("atomium.pdb.parse_header")
    @patch("atomium.pdb.parse_obslte")
    @patch("atomium.pdb.parse_title")
    @patch("atomium.pdb.parse_split")
    @patch("atomium.pdb.parse_caveat")
    @patch("atomium.pdb.parse_keywds")
    @patch("atomium.pdb.parse_expdta")
    @patch("atomium.pdb.parse_mdltyp")
    @patch("atomium.pdb.parse_author")
    @patch("atomium.pdb.parse_revdat")
    @patch("atomium.pdb.parse_sprsde")
    @patch("atomium.pdb.parse_jrnl")
    @patch("atomium.pdb.parse_remarks")
    @patch("atomium.pdb.parse_cryst1")
    @patch("atomium.pdb.parse_origx")
    @patch("atomium.pdb.parse_scale")
    @patch("atomium.pdb.parse_mtrix")
    def test_can_parse_metadata(self, *mocks):
        mmcif = parse_metadata("filestring")
        self.assertEqual(mmcif, {})
        for mock in mocks:
            mock.assert_called_with("filestring", {})



class EntityGettingTests(TestCase):

    @patch("atomium.pdb.parse_compnd_and_source")
    @patch("atomium.pdb.parse_dbref")
    @patch("atomium.pdb.parse_seqadv")
    @patch("atomium.pdb.parse_seqres")
    @patch("atomium.pdb.parse_modres")
    @patch("atomium.pdb.parse_het")
    @patch("atomium.pdb.parse_hetnam")
    @patch("atomium.pdb.parse_hetsyn")
    @patch("atomium.pdb.parse_formul")
    @patch("atomium.pdb.update_entities_from_atoms")
    @patch("atomium.pdb.finalize_entities")
    @patch("atomium.pdb.finalize_polymers")
    @patch("atomium.pdb.add_molecules_to_entities")
    def test_can_get_entity(self, mock_add, mock_finpol, mock_finent, mock_upent, mock_formul, mock_hetsyn, mock_hetnam, mock_het, mock_modres, mock_seqres, mock_seqadv, mock_dbref, mock_cs):
        mock_cs.return_value = ["pe1", "pe2"]
        mock_dbref.return_value = ["p1", "p2", "p3"]
        mock_het.return_value = [["ne1", "ne2"], ["n1", "n2", "n3"]]
        polymers, non_polymers = get_entities("filestring")
        mock_cs.assert_called_with("filestring")
        mock_dbref.assert_called_with("filestring")
        mock_seqadv.assert_called_with("filestring", ["p1", "p2", "p3"])
        mock_seqres.assert_called_with("filestring", ["p1", "p2", "p3"])
        mock_modres.assert_called_with("filestring", ["p1", "p2", "p3"])
        mock_het.assert_called_with("filestring")
        mock_hetnam.assert_called_with("filestring", ["ne1", "ne2"])
        mock_hetsyn.assert_called_with("filestring", ["ne1", "ne2"])
        mock_formul.assert_called_with("filestring", ["ne1", "ne2"])
        mock_upent.assert_called_with(
            "filestring", ["pe1", "pe2"], ["p1", "p2", "p3"],
            ["ne1", "ne2"], ["n1", "n2", "n3"]
        )
        mock_finent.assert_called_with(["pe1", "pe2"], ["ne1", "ne2"])
        mock_finpol.assert_called_with(["p1", "p2", "p3"])
        mock_add.assert_called_with(
            ["pe1", "pe2"], ["p1", "p2", "p3"], ["ne1", "ne2"], ["n1", "n2", "n3"]
        )
        self.assertEqual(polymers, ["pe1", "pe2"])
        self.assertEqual(non_polymers, ["ne1", "ne2"])



class StructureCategoryBuildingTests(TestCase):

    @patch("atomium.pdb.build_entity_category")
    @patch("atomium.pdb.build_entity_name_com")
    @patch("atomium.pdb.build_entity_poly")
    @patch("atomium.pdb.build_entity_poly_seq")
    @patch("atomium.pdb.build_struct_ref")
    @patch("atomium.pdb.build_struct_ref_seq")
    @patch("atomium.pdb.build_struct_ref_seq_dif")
    @patch("atomium.pdb.build_pdbx_struct_mod_residue")
    @patch("atomium.pdb.build_pdbx_entity_nonpoly")
    @patch("atomium.pdb.build_chem_comp")
    @patch("atomium.pdb.build_atom_type")
    @patch("atomium.pdb.build_atom_site")
    @patch("atomium.pdb.build_atom_site_anisotrop")
    @patch("atomium.pdb.update_atom_ids")
    @patch("atomium.pdb.build_struct_asym")
    def test_can_build_structure_categories(self, *mocks):
        mocks = mocks[::-1]
        build_structure_categories("filestring", "polymers", "nonpolymers", {"mmcif": 1})
        mocks[0].assert_called_with("polymers", "nonpolymers", {"mmcif": 1})
        for mock in mocks[1:8]:
            mock.assert_called_with("polymers", {"mmcif": 1})
        for mock in mocks[8:10]:
            mock.assert_called_with("nonpolymers", {"mmcif": 1})
        mocks[10].assert_called_with("filestring", {"mmcif": 1})
        mocks[11].assert_called_with("filestring", "polymers", "nonpolymers", {"mmcif": 1})
        mocks[12].assert_called_with("filestring", {"mmcif": 1})
        mocks[13].assert_called_with({"mmcif": 1})
        mocks[14].assert_called_with({"mmcif": 1})
    


class HeaderParsingTests(TestCase):

    @patch("atomium.pdb.pdb_date_to_mmcif_date")
    def test_can_handle_no_header(self, mock_date):
        mmcif = {}
        mock_date.return_value = None
        parse_header("", mmcif)
        self.assertEqual(mmcif, {
            "entry": [{"id": "XXXX"}],
            "pdbx_database_status": [{"status_code": "REL", "entry_id": "XXXX", "recvd_initial_deposition_date": "?"}],
            "struct_keywords": [{"entry_id": "XXXX", "pdbx_keywords": "?", "text": "?"}]
        })
    

    @patch("atomium.pdb.pdb_date_to_mmcif_date")
    def test_can_handle_empty_header(self, mock_date):
        mmcif = {}
        mock_date.return_value = None
        parse_header("HEADER  ", mmcif)
        self.assertEqual(mmcif, {
            "entry": [{"id": "XXXX"}],
            "pdbx_database_status": [{"status_code": "REL", "entry_id": "XXXX", "recvd_initial_deposition_date": "?"}],
            "struct_keywords": [{"entry_id": "XXXX", "pdbx_keywords": "?", "text": "?"}]
        })


    @patch("atomium.pdb.pdb_date_to_mmcif_date")
    def test_can_parse_keywords(self, mock_date):
        mmcif = {}
        mock_date.return_value = None
        parse_header("HEADER    LYASE        ", mmcif)
        self.assertEqual(mmcif, {
            "entry": [{"id": "XXXX"}],
            "pdbx_database_status": [{"status_code": "REL", "entry_id": "XXXX", "recvd_initial_deposition_date": "?"}],
            "struct_keywords": [{"entry_id": "XXXX", "pdbx_keywords": "LYASE", "text": "?"}]
        })


    @patch("atomium.pdb.pdb_date_to_mmcif_date")
    def test_can_parse_date(self, mock_date):
        mmcif = {}
        mock_date.return_value = "2002-05-06"
        parse_header("HEADER                                            06-MAY-02  ", mmcif)
        self.assertEqual(mmcif, {
            "entry": [{"id": "XXXX"}],
            "pdbx_database_status": [{"status_code": "REL", "entry_id": "XXXX", "recvd_initial_deposition_date": "2002-05-06"}],
            "struct_keywords": [{"entry_id": "XXXX", "pdbx_keywords": "?", "text": "?"}]
        })
        mock_date.assert_called_with("06-MAY-02")


    @patch("atomium.pdb.pdb_date_to_mmcif_date")
    def test_can_parse_code(self, mock_date):
        mmcif = {}
        mock_date.return_value = None
        parse_header("HEADER                                                        1LOL", mmcif)
        self.assertEqual(mmcif, {
            "entry": [{"id": "1LOL"}],
            "pdbx_database_status": [{"status_code": "REL", "entry_id": "1LOL", "recvd_initial_deposition_date": "?"}],
            "struct_keywords": [{"entry_id": "1LOL", "pdbx_keywords": "?", "text": "?"}]
        })
    

    @patch("atomium.pdb.pdb_date_to_mmcif_date")
    def test_can_parse_full_header(self, mock_date):
        mmcif = {}
        mock_date.return_value = "2002-05-06"
        parse_header("HEADER    LYASE                                   06-MAY-02   1LOL", mmcif)
        self.assertEqual(mmcif, {
            "entry": [{"id": "1LOL"}],
            "pdbx_database_status": [{"status_code": "REL", "entry_id": "1LOL", "recvd_initial_deposition_date": "2002-05-06"}],
            "struct_keywords": [{"entry_id": "1LOL", "pdbx_keywords": "LYASE", "text": "?"}]
        })
        mock_date.assert_called_with("06-MAY-02")



class ObslteParsingTests(TestCase):

    def test_can_handle_no_obslte(self):
        mmcif = {}
        parse_obslte("", mmcif)
        self.assertEqual(mmcif, {})
    

    @patch("atomium.pdb.pdb_date_to_mmcif_date")
    def test_can_parse_single_code(self, mock_date):
        mmcif = {}
        mock_date.return_value = "1994-01-21"
        parse_obslte("OBSLTE     31-JAN-94 1MBP      2MBP", mmcif)
        self.assertEqual(mmcif, {"pdbx_database_PDB_obs_spr": [{
            "id": "OBSLTE", "details": "?", "date": "1994-01-21",
            "pdb_id": "2MBP", "replace_pdb_id": "1MBP"
        }]})
        mock_date.assert_called_with("31-JAN-94")
    

    @patch("atomium.pdb.pdb_date_to_mmcif_date")
    def test_can_parse_multiple_codes(self, mock_date):
        mmcif = {}
        mock_date.return_value = "1994-01-21"
        parse_obslte("OBSLTE     31-JAN-94 1MBP      2MBP      3MBP", mmcif)
        self.assertEqual(mmcif, {"pdbx_database_PDB_obs_spr": [{
            "id": "OBSLTE", "details": "?", "date": "1994-01-21",
            "pdb_id": "2MBP 3MBP", "replace_pdb_id": "1MBP"
        }]})
        mock_date.assert_called_with("31-JAN-94")



class TitleParsingTests(TestCase):

    def test_can_handle_no_title(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        parse_title("", mmcif)
        self.assertEqual(mmcif, {
            "entry": [{"id": "1XXX"}],
            "struct": [{
                "entry_id": "1XXX", "title": "?", "pdbx_descriptor": "?",
                "pdbx_model_details": "?", "pdbx_CASP_flag": "?",
                "pdbx_model_type_details": "?"
            }]
        })
    

    def test_can_handle_single_line_title(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        parse_title("TITLE     RHIZOPUSPEPSIN COMPLEXED WITH REDUCED PEPTIDE INHIBITOR   ", mmcif)
        self.assertEqual(mmcif, {
            "entry": [{"id": "1XXX"}],
            "struct": [{
                "entry_id": "1XXX", "pdbx_descriptor": "?",
                "pdbx_model_details": "?", "pdbx_CASP_flag": "?",
                "pdbx_model_type_details": "?",
                "title": "RHIZOPUSPEPSIN COMPLEXED WITH REDUCED PEPTIDE INHIBITOR"
            }]
        })
    

    def test_can_handle_multi_line_title(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        parse_title(
            "TITLE     STRUCTURE OF THE TRANSFORMED MONOCLINIC LYSOZYME BY         \n"
            "TITLE    2 CONTROLLED DEHYDRATION  ",
            mmcif
        )
        self.assertEqual(mmcif, {
            "entry": [{"id": "1XXX"}],
            "struct": [{
                "entry_id": "1XXX", "pdbx_descriptor": "?",
                "pdbx_model_details": "?", "pdbx_CASP_flag": "?",
                "pdbx_model_type_details": "?",
                "title": "STRUCTURE OF THE TRANSFORMED MONOCLINIC LYSOZYME BY CONTROLLED DEHYDRATION"
            }]
        })



class SplitParsingTests(TestCase):

    def test_can_handle_no_split(self):
        mmcif = {}
        parse_split("", mmcif)
        self.assertEqual(mmcif, {})
    

    def test_can_handle_split(self):
        mmcif = {}
        parse_split("SPLIT      1VOQ 1VOR 1VOS 1VOU 1VOV", mmcif)
        self.assertEqual(mmcif, {"pdbx_database_related": [
            {"db_name": "PDB", "db_id": "1VOQ", "content_type": "split", "details": "Split 1"},
            {"db_name": "PDB", "db_id": "1VOR", "content_type": "split", "details": "Split 2"},
            {"db_name": "PDB", "db_id": "1VOS", "content_type": "split", "details": "Split 3"},
            {"db_name": "PDB", "db_id": "1VOU", "content_type": "split", "details": "Split 4"},
            {"db_name": "PDB", "db_id": "1VOV", "content_type": "split", "details": "Split 5"},
        ]})
    

    def test_can_handle_multi_line_split(self):
        mmcif = {}
        parse_split("SPLIT      1VOQ 1VOR 1VOS   \nSPLIT      1VOU 1VOV", mmcif)
        self.assertEqual(mmcif, {"pdbx_database_related": [
            {"db_name": "PDB", "db_id": "1VOQ", "content_type": "split", "details": "Split 1"},
            {"db_name": "PDB", "db_id": "1VOR", "content_type": "split", "details": "Split 2"},
            {"db_name": "PDB", "db_id": "1VOS", "content_type": "split", "details": "Split 3"},
            {"db_name": "PDB", "db_id": "1VOU", "content_type": "split", "details": "Split 4"},
            {"db_name": "PDB", "db_id": "1VOV", "content_type": "split", "details": "Split 5"},
        ]})



class CaveatParsingTests(TestCase):

    def test_can_handle_no_caveat(self):
        mmcif = {}
        parse_caveat("", mmcif)
        self.assertEqual(mmcif, {})
    

    def test_can_handle_caveat(self):
        mmcif = {}
        parse_caveat("CAVEAT     1CK8    THERE ARE CHIRALITY ERRORS IN C-ALPHA CENTERS       ", mmcif)
        self.assertEqual(mmcif, {"database_PDB_caveat": [
            {"id": "1", "text": "THERE ARE CHIRALITY ERRORS IN C-ALPHA CENTERS"},
        ]})
    

    def test_can_handle_multi_line_caveat(self):
        mmcif = {}
        parse_caveat(
            "CAVEAT     4HHB    THR A 137 HAS WRONG CHIRALITY AT ATOM CB THR B 12 HAS WRONG  \n"
            "CAVEAT   2 4HHB    CHIRALITY AT ATOM CB THR B 50 HAS WRONG CHIRALITY AT ATOM   \n"
            "CAVEAT   3 4HHB    CB",
            mmcif
        )
        self.assertEqual(mmcif, {"database_PDB_caveat": [
            {"id": "1", "text": "THR A 137 HAS WRONG CHIRALITY AT ATOM CB THR B 12 HAS WRONG CHIRALITY AT ATOM CB THR B 50 HAS WRONG CHIRALITY AT ATOM CB"},
        ]})



class KeywdsParsingTests(TestCase):

    def test_can_handle_no_keywds(self):
        mmcif = {"struct_keywords": [{"text": "?"}]}
        parse_keywds("", mmcif)
        self.assertEqual(mmcif, {"struct_keywords": [{"text": "?"}]})
    

    def test_can_handle_keywds(self):
        mmcif = {"struct_keywords": [{"text": "?"}]}
        parse_keywds("KEYWDS    LYASE,  TRICARBOXYLIC ACID CYCLE       ", mmcif)
        self.assertEqual(mmcif, {"struct_keywords": [
            {"text": "LYASE,  TRICARBOXYLIC ACID CYCLE"},
        ]})
    

    def test_can_handle_multi_line_keywds(self):
        mmcif = {"struct_keywords": [{"text": "?"}]}
        parse_keywds(
            "KEYWDS    LYASE,  TRICARBOXYLIC ACID CYCLE, MITOCHONDRION, OXIDATIVE \n"
            "KEYWDS   2 METABOLISM",
            mmcif
        )
        self.assertEqual(mmcif, {"struct_keywords": [
            {"text": "LYASE,  TRICARBOXYLIC ACID CYCLE, MITOCHONDRION, OXIDATIVE METABOLISM"},
        ]})



class ExpdtaParsingTests(TestCase):

    def test_can_handle_no_expdta(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        parse_expdta("", mmcif)
        self.assertEqual(mmcif, {"entry": [{"id": "1XXX"}]})
    

    def test_can_parse_expdta(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        parse_expdta("EXPDTA    NEUTRON DIFFRACTION  ", mmcif)
        self.assertEqual(mmcif, {
            "entry": [{"id": "1XXX"}],
            "exptl": [{"entry_id": "1XXX", "method": "NEUTRON DIFFRACTION"}]
        })
    

    def test_can_parse_expdta(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        parse_expdta("EXPDTA    NEUTRON DIFFRACTION; X-RAY DIFFRACTION  ", mmcif)
        self.assertEqual(mmcif, {
            "entry": [{"id": "1XXX"}],
            "exptl": [
                {"entry_id": "1XXX", "method": "NEUTRON DIFFRACTION"},
                {"entry_id": "1XXX", "method": "X-RAY DIFFRACTION"}
            ]
        })



class MdltypParsingTests(TestCase):

    def test_can_handle_no_mdltyp(self):
        mmcif = {"struct": [{"pdbx_model_type_details": "?"}]}
        parse_mdltyp("", mmcif)
        self.assertEqual(mmcif, {"struct": [{"pdbx_model_type_details": "?"}]})
    

    def test_can_parse_mdltyp(self):
        mmcif = {"struct": [{"pdbx_model_type_details": "?"}]}
        parse_mdltyp("MDLTYP    MINIMIZED AVERAGE", mmcif)
        self.assertEqual(mmcif, {
            "struct": [{"pdbx_model_type_details": "MINIMIZED AVERAGE"}]
        })


    def test_can_parse_mdltyp_with_asyms(self):
        mmcif = {"struct": [{"pdbx_model_type_details": "?"}]}
        parse_mdltyp("MDLTYP    MINIMIZED AVERAGE ; OTHER ; CA ATOMS ONLY, CHAIN A, B", mmcif)
        self.assertEqual(mmcif, {
            "struct": [{"pdbx_model_type_details": "MINIMIZED AVERAGE ; OTHER"}],
            "pdbx_coordinate_model": [
                {"asym_id": "A", "type": "CA ATOMS ONLY"},
                {"asym_id": "B", "type": "CA ATOMS ONLY"},
            ]
        })
        

    def test_can_parse_multi_line_mdltyp_with_asyms(self):
        mmcif = {"struct": [{"pdbx_model_type_details": "?"}]}
        parse_mdltyp(
            "MDLTYP    CA ATOMS ONLY, CHAIN A, B, C, D, E, F, G, H, I, J, K ; P ATOMS ONLY,\n"
            "MDLTYP   2 CHAIN X, Y, Z",
            mmcif
        )
        self.assertEqual(mmcif, {
            "struct": [{"pdbx_model_type_details": "?"}],
            "pdbx_coordinate_model": [
                {"asym_id": "A", "type": "CA ATOMS ONLY"},
                {"asym_id": "B", "type": "CA ATOMS ONLY"},
                {"asym_id": "C", "type": "CA ATOMS ONLY"},
                {"asym_id": "D", "type": "CA ATOMS ONLY"},
                {"asym_id": "E", "type": "CA ATOMS ONLY"},
                {"asym_id": "F", "type": "CA ATOMS ONLY"},
                {"asym_id": "G", "type": "CA ATOMS ONLY"},
                {"asym_id": "H", "type": "CA ATOMS ONLY"},
                {"asym_id": "I", "type": "CA ATOMS ONLY"},
                {"asym_id": "J", "type": "CA ATOMS ONLY"},
                {"asym_id": "K", "type": "CA ATOMS ONLY"},
                {"asym_id": "X", "type": "P ATOMS ONLY"},
                {"asym_id": "Y", "type": "P ATOMS ONLY"},
                {"asym_id": "Z", "type": "P ATOMS ONLY"},
            ]
        })



class AuthorParsingTests(TestCase):

    def test_can_handle_no_author(self):
        mmcif = {}
        parse_author("", mmcif)
        self.assertEqual(mmcif, {})


    @patch("atomium.pdb.pdb_names_to_mmcif_names")
    def test_can_parse_author(self, mock_names):
        mock_names.return_value = ["Berry, M.B"]
        mmcif = {}
        parse_author("AUTHOR    M.B.BERRY  ", mmcif)
        self.assertEqual(mmcif, {
            "audit_author": [{"name": "Berry, M.B", "pdbx_ordinal": "1"}]
        })
        mock_names.assert_called_with(["M.B.BERRY"])
    

    @patch("atomium.pdb.pdb_names_to_mmcif_names")
    def test_can_parse_multi_line_author(self, mock_names):
        mock_names.return_value = ["Berry, M.B", "Phillips, G.N"]
        mmcif = {}
        parse_author("AUTHOR    M.B.BERRY  \nAUTHOR   2 G.N.PHILLIPS", mmcif)
        self.assertEqual(mmcif, {
            "audit_author": [
                {"name": "Berry, M.B", "pdbx_ordinal": "1"},
                {"name": "Phillips, G.N", "pdbx_ordinal": "2"},
            ]
        })
        mock_names.assert_called_with(["M.B.BERRY", "G.N.PHILLIPS"])



class RevdatParsingTests(TestCase):

    def test_can_handle_no_revdat(self):
        mmcif = {}
        parse_revdat("", mmcif)
        self.assertEqual(mmcif, {})
    

    @patch("atomium.pdb.pdb_date_to_mmcif_date")
    def test_can_parse_revdat(self, mock_date):
        mmcif = {}
        mock_date.side_effect = ["1984-07-17", "2011-07-13","2020-06-17",  "2021-03-31"]
        parse_revdat(
            "REVDAT   4   31-MAR-21 4HHB    1       REMARK ATOM\n"
            "REVDAT   3   17-JUN-20 4HHB    1       CAVEAT COMPND SOURCE DBREF\n"
            "REVDAT   3 2                   1       ATOM\n"
            "REVDAT   2   13-JUL-11 4HHB    1       VERSN\n"
            "REVDAT   1   17-JUL-84 4HHB    0\n",
            mmcif
        )
        self.assertEqual(mmcif, {
            "pdbx_audit_revision_history": [{
                "ordinal": "1", "data_content_type": "Structure model",
                "major_revision": "1", "minor_revision": "0",
                "revision_date": "1984-07-17"
            }, {
                "ordinal": "2", "data_content_type": "Structure model",
                "major_revision": "2", "minor_revision": "0",
                "revision_date": "2011-07-13"
            }, {
                "ordinal": "3", "data_content_type": "Structure model",
                "major_revision": "3", "minor_revision": "0",
                "revision_date": "2020-06-17"
            }, {
                "ordinal": "4", "data_content_type": "Structure model",
                "major_revision": "4", "minor_revision": "0",
                "revision_date": "2021-03-31"
            }]
        })
        mock_date.assert_any_call("17-JUL-84")
        mock_date.assert_any_call("13-JUL-11")
        mock_date.assert_any_call("17-JUN-20")
        mock_date.assert_any_call("31-MAR-21")



class SprsdeParsingTests(TestCase):

    def test_can_handle_no_sprsde(self):
        mmcif = {}
        parse_sprsde("", mmcif)
        self.assertEqual(mmcif, {})
    

    @patch("atomium.pdb.pdb_date_to_mmcif_date")
    def test_can_parse_single_code(self, mock_date):
        mmcif = {}
        mock_date.return_value = "1994-01-21"
        parse_sprsde("SPRSDE     31-JAN-94 4HHB      1HHB", mmcif)
        self.assertEqual(mmcif, {"pdbx_database_PDB_obs_spr": [{
            "id": "SPRSDE", "details": "?", "date": "1994-01-21",
            "pdb_id": "4HHB", "replace_pdb_id": "1HHB"
        }]})
        mock_date.assert_called_with("31-JAN-94")
    

    @patch("atomium.pdb.pdb_date_to_mmcif_date")
    def test_can_parse_multiple_codes(self, mock_date):
        mmcif = {"pdbx_database_PDB_obs_spr": [1]}
        mock_date.return_value = "1994-01-21"
        parse_sprsde("SPRSDE     31-JAN-94 1GDJ      1LH4 2LH4", mmcif)
        self.assertEqual(mmcif, {"pdbx_database_PDB_obs_spr": [{
            "id": "SPRSDE", "details": "?", "date": "1994-01-21",
            "pdb_id": "1GDJ", "replace_pdb_id": "1LH4 2LH4"
        }, 1]})
        mock_date.assert_called_with("31-JAN-94")



class JrnlParsingTests(TestCase):

    def setUp(self):
        self.patch1 = patch("atomium.pdb.pdb_names_to_mmcif_names")
        self.patch2 = patch("atomium.pdb.parse_journal_title")
        self.patch3 = patch("atomium.pdb.parse_journal_references")
        self.patch4 = patch("atomium.pdb.parse_journal_ids")
        self.mock_names = self.patch1.start()
        self.mock_title = self.patch2.start()
        self.mock_refs = self.patch3.start()
        self.mock_ids = self.patch4.start()
        self.mock_names.return_value = ["name1", "name2"]
    

    def tearDown(self):
        self.patch1.stop()
        self.patch2.stop()
        self.patch3.stop()
        self.patch4.stop()


    def test_can_handle_no_jrnl(self):
        mmcif = {}
        parse_jrnl("", mmcif)
        self.assertEqual(mmcif, {})
    

    def test_can_parse_jrnl_auth(self):
        mmcif = {}
        self.mock_names.side_effect = [["name1", "name2"], None]
        lines = [
            "JRNL        AUTH   M.HYVONEN,M.J.MACIAS,M.NILGES,H.OSCHKINAT,   ",
            "JRNL        AUTH 2 M.SARASTE,M.WILMANNS ",
            "ATOM    123",
        ]
        parse_jrnl("\n".join(lines), mmcif)
        self.assertEqual(mmcif, {
            "citation_author": [
                {"citation_id": "primary", "name": "name1", "pdbx_ordinal": "1"},
                {"citation_id": "primary", "name": "name2", "pdbx_ordinal": "2"},
            ],
            "citation": [{
                "id": "primary", "title": "?", "journal_abbrev": "?", "journal_volume": "?",
                "page_first": "?", "page_last": "?", "year": "?", "journal_id_ASTM": "?",
                "country": "?", "journal_id_ISSN": "?", "journal_id_CSD": "?",
                "book_publisher": "?", "pdbx_database_id_PubMed": "?", "pdbx_database_id_DOI": "?"
            }]
        })
        self.mock_names.assert_any_call(["M.HYVONEN,M.J.MACIAS,M.NILGES,H.OSCHKINAT,", "M.SARASTE,M.WILMANNS"])
        self.mock_names.assert_any_call([])
        self.mock_title.assert_called_with(lines[:2], mmcif)
        self.mock_refs.assert_called_with(lines[:2], mmcif)
        self.mock_ids.assert_called_with(lines[:2], mmcif)
    

    def test_can_parse_jrnl_edit(self):
        mmcif = {}
        self.mock_names.side_effect = [None, ["name1", "name2"]]
        lines = [
            "JRNL        EDIT   M.HYVONEN,M.J.MACIAS,M.NILGES,H.OSCHKINAT,   ",
            "JRNL        EDIT 2 M.SARASTE,M.WILMANNS ",
            "ATOM    123",
        ]
        parse_jrnl("\n".join(lines), mmcif)
        self.assertEqual(mmcif, {
            "citation_editor": [
                {"citation_id": "primary", "name": "name1", "pdbx_ordinal": "1"},
                {"citation_id": "primary", "name": "name2", "pdbx_ordinal": "2"},
            ],
            "citation": [{
                "id": "primary", "title": "?", "journal_abbrev": "?", "journal_volume": "?",
                "page_first": "?", "page_last": "?", "year": "?", "journal_id_ASTM": "?",
                "country": "?", "journal_id_ISSN": "?", "journal_id_CSD": "?",
                "book_publisher": "?", "pdbx_database_id_PubMed": "?", "pdbx_database_id_DOI": "?"
            }]
        })
        self.mock_names.assert_any_call(["M.HYVONEN,M.J.MACIAS,M.NILGES,H.OSCHKINAT,", "M.SARASTE,M.WILMANNS"])
        self.mock_names.assert_any_call([])
        self.mock_title.assert_called_with(lines[:2], mmcif)
        self.mock_refs.assert_called_with(lines[:2], mmcif)
        self.mock_ids.assert_called_with(lines[:2], mmcif)
    

    def test_can_parse_jrnl_full(self):
        mmcif = {}
        self.mock_names.side_effect = [["name1", "name2"], ["name3", "name4"]]
        lines = [
            "JRNL        AUTH   M.HYVONEN,M.J.MACIAS,M.NILGES,H.OSCHKINAT,   ",
            "JRNL        AUTH 2 M.SARASTE,M.WILMANNS ",
            "JRNL        EDIT   L.HYVONEN,M.J.MACIAS,M.NILGES,H.OSCHKINAT,   ",
            "JRNL        EDIT 2 L.SARASTE,M.WILMANNS ",
            "ATOM    123",
        ]
        parse_jrnl("\n".join(lines), mmcif)
        self.assertEqual(mmcif, {
            "citation_author": [
                {"citation_id": "primary", "name": "name1", "pdbx_ordinal": "1"},
                {"citation_id": "primary", "name": "name2", "pdbx_ordinal": "2"},
            ],
            "citation_editor": [
                {"citation_id": "primary", "name": "name3", "pdbx_ordinal": "1"},
                {"citation_id": "primary", "name": "name4", "pdbx_ordinal": "2"},
            ],
            "citation": [{
                "id": "primary", "title": "?", "journal_abbrev": "?", "journal_volume": "?",
                "page_first": "?", "page_last": "?", "year": "?", "journal_id_ASTM": "?",
                "country": "?", "journal_id_ISSN": "?", "journal_id_CSD": "?",
                "book_publisher": "?", "pdbx_database_id_PubMed": "?", "pdbx_database_id_DOI": "?"
            }]
        })
        self.mock_names.assert_any_call(["M.HYVONEN,M.J.MACIAS,M.NILGES,H.OSCHKINAT,", "M.SARASTE,M.WILMANNS"])
        self.mock_names.assert_any_call(["L.HYVONEN,M.J.MACIAS,M.NILGES,H.OSCHKINAT,", "L.SARASTE,M.WILMANNS"])
        self.mock_title.assert_called_with(lines[:4], mmcif)
        self.mock_refs.assert_called_with(lines[:4], mmcif)
        self.mock_ids.assert_called_with(lines[:4], mmcif)



class JournalTitleParsingTests(TestCase):

    def test_can_handle_no_title(self):
        mmcif = {"citation": [{"title": "?"}]}
        parse_journal_title(["JRNL    "], mmcif)
        self.assertEqual(mmcif, {"citation": [{"title": "?"}]})
    

    def test_can_parse_title(self):
        mmcif = {"citation": [{"title": "?"}]}
        parse_journal_title([
            "JRNL        TITL   STRUCTURE OF THE BINDING SITE FOR INOSITOL    ",
            "JRNL        TITL 2 PHOSPHATES IN A PH DOMAIN.",
            "JRNL        AUTH NAME.",
        ], mmcif)
        self.assertEqual(mmcif, {"citation": [{"title": "STRUCTURE OF THE BINDING SITE FOR INOSITOL PHOSPHATES IN A PH DOMAIN."}]})



class JournalReferencesParsingTests(TestCase):

    def test_can_handle_no_ref(self):
        mmcif = {}
        parse_journal_references(["JRNL    "], mmcif)
        self.assertEqual(mmcif, {})
    

    def test_can_parse_ref(self):
        mmcif = {"citation": [{}]}
        parse_journal_references([
            "JRNL        REF    PROC. NATL. ACAD. SCI.        V. 114  7037 2017 ",
            "JRNL        REF  2 U.S.A",
            "JRNL        AUTH NAME.",
        ], mmcif)
        self.assertEqual(mmcif, {"citation": [{
            "journal_abbrev": "Proc. Natl. Acad. Sci. U.S.A",
            "journal_volume": "114",
            "page_first": "7037",
            "year": "2017",
        }]})
    

    def test_can_parse_refn(self):
        mmcif = {"citation": [{}]}
        parse_journal_references([
            "JRNL        REFN                   ESSN 1091-6490 ",
            "JRNL        AUTH NAME.",
        ], mmcif)
        self.assertEqual(mmcif, {"citation": [{"journal_id_ISSN": "1091-6490",}]})
    

    def test_can_parse_publ(self):
        mmcif = {"citation": [{}]}
        parse_journal_references([
            "JRNL        PUBL   NATIONAL BIOMEDICAL RESEARCH FOUNDATION, SILVER ",
            "JRNL        PUBL 2 SPRING,MD.       ",
            "JRNL        AUTH NAME.",
        ], mmcif)
        self.assertEqual(mmcif, {"citation": [{"book_publisher": "NATIONAL BIOMEDICAL RESEARCH FOUNDATION, SILVER SPRING,MD."}]})
    

    def test_can_parse_full_references(self):
        mmcif = {"citation": [{}]}
        parse_journal_references([
            "JRNL        REF    PROC. NATL. ACAD. SCI.        V. 114  7037 2017 ",
            "JRNL        REF  2 U.S.A",
            "JRNL        REFN                   ESSN 1091-6490 ",
            "JRNL        PUBL   NATIONAL BIOMEDICAL RESEARCH FOUNDATION, SILVER ",
            "JRNL        PUBL 2 SPRING,MD.       ",
            "JRNL        AUTH NAME.",
        ], mmcif)
        self.assertEqual(mmcif, {"citation": [{
            "journal_abbrev": "Proc. Natl. Acad. Sci. U.S.A",
            "journal_volume": "114",
            "page_first": "7037",
            "year": "2017",
            "journal_id_ISSN": "1091-6490",
            "book_publisher": "NATIONAL BIOMEDICAL RESEARCH FOUNDATION, SILVER SPRING,MD."
        }]})



class JournalIdsParsingTests(TestCase):

    def test_can_handle_no_ids(self):
        mmcif = {}
        parse_journal_ids(["JRNL    "], mmcif)
        self.assertEqual(mmcif, {})
    

    def test_can_parse_pmid(self):
        mmcif = {"citation": [{}]}
        parse_journal_ids([
            "JRNL        PMID   6726807 ",
            "JRNL        AUTH NAME.",
        ], mmcif)
        self.assertEqual(mmcif, {"citation": [{"pdbx_database_id_PubMed": "6726807"}]})
    
    
    def test_can_parse_doi(self):
        mmcif = {"citation": [{}]}
        parse_journal_ids([
            "JRNL        DOI    10.1016/0022-2836(84)90472-8 ",
            "JRNL        AUTH NAME.",
        ], mmcif)
        self.assertEqual(mmcif, {"citation": [{
            "pdbx_database_id_DOI": "10.1016/0022-2836(84)90472-8"
        }]})
    

    def test_can_parse_full_ids(self):
        mmcif = {"citation": [{}]}
        parse_journal_ids([
            "JRNL        PMID   6726807 ",
            "JRNL        DOI    10.1016/0022-2836(84)90472-8 ",
            "JRNL        AUTH NAME.",
        ], mmcif)
        self.assertEqual(mmcif, {"citation": [{
            "pdbx_database_id_PubMed": "6726807",
            "pdbx_database_id_DOI": "10.1016/0022-2836(84)90472-8"
        }]})



class RemarkParsingTests(TestCase):

    @patch("atomium.pdb.parse_remark_2")
    @patch("atomium.pdb.parse_remark_3")
    @patch("atomium.pdb.parse_remark_350")
    @patch("atomium.pdb.parse_remark_800")
    def test_can_parse_remark_records(self, *mocks):
        parse_remarks("filestring", {"mmicf": 1})
        for mock in mocks:
            mock.assert_called_with("filestring", {"mmicf": 1})



class Remark2ParsingTests(TestCase):

    def test_can_handle_no_remark_2(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        parse_remark_2("", mmcif)
        self.assertEqual(mmcif, {"entry": [{"id": "1XXX"}]})
    

    def test_can_get_resolution(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        parse_remark_2("REMARK   2   \nREMARK   2 RESOLUTION.    1.74 ANGSTROMS.", mmcif)
        self.assertEqual(mmcif, {"entry": [{"id": "1XXX"}], "reflns": [{
            "entry_id": "1XXX", "observed_criterion_sigma_I": "?", "observed_criterion_sigma_F": "?",
            "d_resolution_low": "?", "d_resolution_high": "1.74", "number_obs": "?", "number_all": "?",
            "percent_possible_obs": "?", "pdbx_Rmerge_I_obs": "?", "pdbx_Rsym_value": "?", "pdbx_netI_over_sigmaI": "?",
            "B_iso_Wilson_estimate": "?", "pdbx_redundancy": "?", "R_free_details": "?", "limit_h_max": "?",
            "limit_h_min": "?", "limit_k_max": "?", "limit_k_min": "?", "limit_l_max": "?", "limit_l_min": "?",
            "observed_criterion_F_max": "?", "observed_criterion_F_min": "?", "pdbx_chi_squared": "?",
            "pdbx_scaling_rejects": "?", "pdbx_diffrn_id": "?", "pdbx_ordinal": "?"
        }]})
    

    def test_can_get_na_resolution(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        parse_remark_2("REMARK   2   \nREMARK   2 RESOLUTION.   NOT APPLICABLE. ", mmcif)
        self.assertEqual(mmcif, {"entry": [{"id": "1XXX"}]})



class Remark3ParsingTests(TestCase):

    def test_can_handle_no_remark_3(self):
        mmcif = {"reflns": [{}]}
        parse_remark_3("", mmcif)
        self.assertEqual(mmcif, {"reflns": [{}]})
    

    @patch("atomium.pdb.update_reflns_from_remark_3")
    @patch("atomium.pdb.update_refine_from_remark_3")
    def test_can_parse_remark_3(self, mock_refine, mock_reflns):
        mmcif = {"reflns": [{}]}
        parse_remark_3("HEADER \nREMARK   3   \nREMARK   3 XXX.\nATOM", mmcif)
        mock_refine.assert_called_with(["REMARK   3   ", "REMARK   3 XXX."], {"reflns": [{}]})
        mock_reflns.assert_called_with(["REMARK   3   ", "REMARK   3 XXX."], {"reflns": [{}]})



class ReflnsFromRemark3Tests(TestCase):

    def test_nothing_to_get(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        update_reflns_from_remark_3([
            "REMARK   3  B VALUES.    ",
            "REMARK   3   OVERALL ANISOTROPIC B VALUE.       "
        ], mmcif)
        self.assertEqual(mmcif, {"entry": [{"id": "1XXX"}]})
    

    def test_null_value(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        update_reflns_from_remark_3([
            "REMARK   3  B VALUES.    ",
            "REMARK   3   FROM WILSON PLOT           (A**2) : NULL   ",
            "REMARK   3   MEAN B VALUE      (OVERALL, A**2) : 24.20    ",
            "REMARK   3   OVERALL ANISOTROPIC B VALUE.       "
        ], mmcif)
        self.assertEqual(mmcif, {"entry": [{"id": "1XXX"}]})
    

    def test_can_create_table(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        update_reflns_from_remark_3([
            "REMARK   3  B VALUES.    ",
            "REMARK   3   FROM WILSON PLOT           (A**2) : 11.20   ",
            "REMARK   3   MEAN B VALUE      (OVERALL, A**2) : 24.20    ",
            "REMARK   3   OVERALL ANISOTROPIC B VALUE.       "
        ], mmcif)
        self.assertEqual(mmcif, {
            "entry": [{"id": "1XXX"}], "reflns": [{
                "entry_id": "1XXX", "observed_criterion_sigma_I": "?", "observed_criterion_sigma_F": "?",
                "d_resolution_low": "?", "d_resolution_high": "?", "number_obs": "?", "number_all": "?",
                "percent_possible_obs": "?", "pdbx_Rmerge_I_obs": "?", "pdbx_Rsym_value": "?", "pdbx_netI_over_sigmaI": "?",
                "B_iso_Wilson_estimate": "11.20", "pdbx_redundancy": "?", "R_free_details": "?", "limit_h_max": "?",
                "limit_h_min": "?", "limit_k_max": "?", "limit_k_min": "?", "limit_l_max": "?", "limit_l_min": "?",
                "observed_criterion_F_max": "?", "observed_criterion_F_min": "?", "pdbx_chi_squared": "?",
                "pdbx_scaling_rejects": "?", "pdbx_diffrn_id": "?", "pdbx_ordinal": "?"
            }]
        })


    def test_can_update_table(self):
        mmcif = {"entry": [{"id": "1XXX"}], "reflns": [{}]}
        update_reflns_from_remark_3([
            "REMARK   3  B VALUES.    ",
            "REMARK   3   FROM WILSON PLOT           (A**2) : 11.20   ",
            "REMARK   3   MEAN B VALUE      (OVERALL, A**2) : 24.20    ",
            "REMARK   3   OVERALL ANISOTROPIC B VALUE.       "
        ], mmcif)
        self.assertEqual(mmcif, {
            "entry": [{"id": "1XXX"}], "reflns": [{"B_iso_Wilson_estimate": "11.20"}]
        })



class RefineFromRemark3Tests(TestCase):

    def test_nothing_to_get(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        update_refine_from_remark_3([
            "REMARK   3  B VALUES.    ",
            "REMARK   3   OVERALL ANISOTROPIC B VALUE.       "
        ], mmcif)
        self.assertEqual(mmcif, {"entry": [{"id": "1XXX"}]})
    

    def test_can_get_values(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        lines = [
            "REMARK   3                                                      ",
            "REMARK   3  REFINEMENT TARGET : ENGH & HUBER                    ",
            "REMARK   3                                                      ",
            "REMARK   3  DATA USED IN REFINEMENT.                            ",
            "REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 1.90           ",
            "REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 27.07          ",
            "REMARK   3   DATA CUTOFF HIGH         (ABS(F)) : 679650.230     ",
            "REMARK   3   DATA CUTOFF LOW          (ABS(F)) : 0.0000         ",
            "REMARK   3   COMPLETENESS (WORKING+TEST)   (%) : 97.2           ",
            "REMARK   3   NUMBER OF REFLECTIONS             : 32092          ",
            "REMARK   3                                                      ",
            "REMARK   3  FIT TO DATA USED IN REFINEMENT.                     ",
            "REMARK   3   CROSS-VALIDATION METHOD          : THROUGHOUT      ",
            "REMARK   3   FREE R VALUE TEST SET SELECTION  : RANDOM          ",
            "REMARK   3   R VALUE            (WORKING SET) : 0.193           ",
            "REMARK   3   FREE R VALUE                     : 0.229           ",
            "REMARK   3   FREE R VALUE TEST SET SIZE   (%) : 4.900           ",
            "REMARK   3   FREE R VALUE TEST SET COUNT      : 1583            ",
            "REMARK   3   ESTIMATED ERROR OF FREE R VALUE  : 0.006           ",
            "REMARK   3                                                      ",
            "REMARK   3  FIT IN THE HIGHEST RESOLUTION BIN.                  ",
            "REMARK   3   BIN RESOLUTION RANGE HIGH       (A) : 1.90         ",
            "REMARK   3   BIN RESOLUTION RANGE LOW        (A) : 2.02         ",
            "REMARK   3   BIN COMPLETENESS (WORKING+TEST) (%) : 96.40        ",
            "REMARK   3   REFLECTIONS IN BIN    (WORKING SET) : 4989         ",
            "REMARK   3   BIN R VALUE           (WORKING SET) : 0.2340       ",
            "REMARK   3   BIN FREE R VALUE                    : 0.2620       ",
            "REMARK   3   BIN FREE R VALUE TEST SET SIZE  (%) : 5.30         ",
            "REMARK   3   BIN FREE R VALUE TEST SET COUNT     : 279          ",
            "REMARK   3   ESTIMATED ERROR OF BIN FREE R VALUE : 0.016        ",
            "REMARK   3                                                      ",
            "REMARK   3  B VALUES.                                           ",
            "REMARK   3   FROM WILSON PLOT           (A**2) : 11.20          ",
            "REMARK   3   MEAN B VALUE      (OVERALL, A**2) : 24.20          ",
            "REMARK   3   OVERALL ANISOTROPIC B VALUE.                       ",
            "REMARK   3    B11 (A**2) : 7.19000                              ",
            "REMARK   3    B22 (A**2) : -3.85000                             ",
            "REMARK   3    B33 (A**2) : -3.34000                             ",
            "REMARK   3    B12 (A**2) : 0.00000                              ",
            "REMARK   3    B13 (A**2) : 3.48000                              ",
            "REMARK   3    B23 (A**2) : 0.00000                              ",
            "REMARK   3                                                      ",
            "REMARK   3  ISOTROPIC THERMAL MODEL : RESTRAINED                ",
            "REMARK   3                                                      ",
            "REMARK   3  BULK SOLVENT MODELING.                              ",
            "REMARK   3   METHOD USED : FLAT MODEL                           ",
            "REMARK   3   KSOL        : 0.39                                 ",
            "REMARK   3   BSOL        : 54.02                                ",
        ]
        update_refine_from_remark_3(lines, mmcif)
        self.assertEqual(mmcif, {
            "entry": [{"id": "1XXX"}], "refine": [{
                "entry_id": "1XXX", "ls_number_reflns_obs": "32092", "ls_number_reflns_all": "?",
                "pdbx_ls_sigma_I": "?", "pdbx_ls_sigma_F": "?", "pdbx_data_cutoff_high_absF": "679650.230",
                "pdbx_data_cutoff_low_absF": "0.0000", "ls_d_res_low": "27.07", "ls_d_res_high": "1.90",
                "ls_percent_reflns_obs": "97.2", "ls_R_factor_obs": "0.193", "ls_R_factor_all": "0.193",
                "ls_R_factor_R_work": "0.193", "ls_R_factor_R_free": "0.229", "ls_R_factor_R_free_error": "0.006",
                "ls_R_factor_R_free_error_details": "?", "ls_percent_reflns_R_free": "4.900",
                "ls_number_reflns_R_free": "1583", "ls_number_parameters": "?", "ls_number_restraints": "?",
                "occupancy_min": "?", "occupancy_max": "?", "correlation_coeff_Fo_to_Fc": "?",
                "correlation_coeff_Fo_to_Fc_free": "?", "B_iso_mean": "24.20", "aniso_B[1][1]": "7.19000",
                "aniso_B[2][2]": "-3.85000", "aniso_B[3][3]": "-3.34000", "aniso_B[1][2]": "0.00000",
                "aniso_B[1][3]": "3.48000", "aniso_B[2][3]": "0.00000", "solvent_model_details": "FLAT MODEL",
                "solvent_model_param_ksol": "0.39", "solvent_model_param_bsol": "54.02",
                "pdbx_solvent_vdw_probe_radii": "?", "pdbx_solvent_ion_probe_radii": "?",
                "pdbx_solvent_shrinkage_radii": "?", "pdbx_ls_cross_valid_method": "THROUGHOUT",
                "details": "?", "pdbx_starting_model": "?", "pdbx_method_to_determine_struct": "?",
                "pdbx_isotropic_thermal_model": "RESTRAINED", "pdbx_stereochemistry_target_values": "ENGH & HUBER",
                "pdbx_stereochem_target_val_spec_case": "?", "pdbx_R_Free_selection_details": "RANDOM",
                "pdbx_overall_ESU_R_Free": "?", "overall_SU_B": "?", "ls_redundancy_reflns_obs": "?",
                "B_iso_min": "?", "B_iso_max": "?", "overall_SU_R_Cruickshank_DPI": "?", "overall_SU_R_free": "?",
                "overall_SU_ML": "?", "pdbx_overall_ESU_R": "?", "pdbx_data_cutoff_high_rms_absF": "679650.230",
                "pdbx_refine_id": "?", "pdbx_overall_phase_error": "?", "ls_wR_factor_R_free": "?",
                "ls_wR_factor_R_work": "?", "overall_FOM_free_R_set": "?", "overall_FOM_work_R_set": "?",
                "pdbx_diffrn_id": "?", "pdbx_TLS_residual_ADP_flag": "?", "pdbx_overall_SU_R_free_Cruickshank_DPI": "?",
                "pdbx_overall_SU_R_Blow_DPI": "?", "pdbx_overall_SU_R_free_Blow_DPI": "?"
            }]
        })


    def test_nothing_to_get(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        lines = [
            "REMARK   3                                                 ",
            "REMARK   3  REFINEMENT TARGET : NULL                       ",
            "REMARK   3                                                 ",
            "REMARK   3  DATA USED IN REFINEMENT.                       ",
            "REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : NULL      ",
            "REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : NULL      ",
            "REMARK   3   DATA CUTOFF HIGH         (ABS(F)) : NULL      ",
            "REMARK   3   DATA CUTOFF LOW          (ABS(F)) : NULL      ",
            "REMARK   3   COMPLETENESS (WORKING+TEST)   (%) : NULL      ",
            "REMARK   3   NUMBER OF REFLECTIONS             : NULL      ",
            "REMARK   3                                                 ",
            "REMARK   3  FIT TO DATA USED IN REFINEMENT.                ",
            "REMARK   3   CROSS-VALIDATION METHOD          : NULL       ",
            "REMARK   3   FREE R VALUE TEST SET SELECTION  : NULL       ",
            "REMARK   3   R VALUE            (WORKING SET) : NULL       ",
            "REMARK   3   FREE R VALUE                     : NULL       ",
            "REMARK   3   FREE R VALUE TEST SET SIZE   (%) : NULL       ",
            "REMARK   3   FREE R VALUE TEST SET COUNT      : NULL       ",
            "REMARK   3   ESTIMATED ERROR OF FREE R VALUE  : NULL       ",
            "REMARK   3                                                 ",
            "REMARK   3  FIT IN THE HIGHEST RESOLUTION BIN.             ",
            "REMARK   3   BIN RESOLUTION RANGE HIGH       (A) : 1.90    ",
            "REMARK   3   BIN RESOLUTION RANGE LOW        (A) : 2.02    ",
            "REMARK   3   BIN COMPLETENESS (WORKING+TEST) (%) : 96.40   ",
            "REMARK   3   REFLECTIONS IN BIN    (WORKING SET) : 4989    ",
            "REMARK   3   BIN R VALUE           (WORKING SET) : 0.2340  ",
            "REMARK   3   BIN FREE R VALUE                    : 0.2620  ",
            "REMARK   3   BIN FREE R VALUE TEST SET SIZE  (%) : 5.30    ",
            "REMARK   3   BIN FREE R VALUE TEST SET COUNT     : 279     ",
            "REMARK   3   ESTIMATED ERROR OF BIN FREE R VALUE : 0.016   ",
            "REMARK   3                                                 ",
            "REMARK   3  B VALUES.                                      ",
            "REMARK   3   FROM WILSON PLOT           (A**2) : 11.20     ",
            "REMARK   3   MEAN B VALUE      (OVERALL, A**2) : NULL      ",
            "REMARK   3   OVERALL ANISOTROPIC B VALUE.                  ",
            "REMARK   3    B11 (A**2) : NULL                            ",
            "REMARK   3    B22 (A**2) : NULL                            ",
            "REMARK   3    B33 (A**2) : NULL                            ",
            "REMARK   3    B12 (A**2) : NULL                            ",
            "REMARK   3    B13 (A**2) : NULL                            ",
            "REMARK   3    B23 (A**2) : NULL                            ",
            "REMARK   3                                                 ",
            "REMARK   3  ISOTROPIC THERMAL MODEL : NULL                 ",
            "REMARK   3                                                 ",
            "REMARK   3  BULK SOLVENT MODELING.                         ",
            "REMARK   3   METHOD USED : NULL                            ",
            "REMARK   3   KSOL        : NULL                            ",
            "REMARK   3   BSOL        : NULL                            ",
        ]
        update_refine_from_remark_3(lines, mmcif)
        self.assertEqual(mmcif, {
            "entry": [{"id": "1XXX"}]
        })



class Remark350ParsingTests(TestCase):

    def test_can_handle_no_remark_350(self):
        mmcif = {}
        parse_remark_350("", mmcif)
        self.assertEqual(mmcif, {})
    

    @patch("atomium.pdb.parse_assembly")
    def test_can_parse_remark_350(self, mock_parse):
        lines = (
            "REMARK 350                                                           \n"
            "REMARK 350 GIVEN BELOW.  BOTH NON-CRYSTALLOGRAPHIC AND               \n"
            "REMARK 350 CRYSTALLOGRAPHIC OPERATIONS ARE GIVEN.                    \n"
            "REMARK 350                                                           \n"
            "REMARK 350 BIOMOLECULE: 1                                            \n"
            "REMARK 350 SOFTWARE USED: PISA                                       \n"
            "REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000 \n"
            "REMARK 350                                                           \n"
            "REMARK 350 BIOMOLECULE: 2                                            \n"
            "REMARK 350 APPLY THE FOLLOWING TO CHAINS: C, D                       \n"
            "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000 \n"
            "REMARK 350                                                           \n"
            "REMARK 350 BIOMOLECULE: 3                                            \n"
            "REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: DIMERIC                \n"
            "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000 \n"
        )
        mmcif = {}
        mock_parse.side_effect = lambda l, m: m["pdbx_struct_assembly"].append(1)
        parse_remark_350(lines, mmcif)
        mock_parse.assert_any_call([
            "REMARK 350 BIOMOLECULE: 1                                            ",
            "REMARK 350 SOFTWARE USED: PISA                                       ",
            "REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000 ",
            "REMARK 350                                                           ",
        ], mmcif)
        mock_parse.assert_any_call([
            "REMARK 350 BIOMOLECULE: 2                                            ",
            "REMARK 350 APPLY THE FOLLOWING TO CHAINS: C, D                       ",
            "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000 ",
            "REMARK 350                                                           ",
        ], mmcif)
        mock_parse.assert_any_call([
            "REMARK 350 BIOMOLECULE: 3                                            ",
            "REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: DIMERIC                ",
            "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000 ",
        ], mmcif)
        self.assertEqual(mmcif, {"pdbx_struct_assembly": [1, 1, 1]})



class AssemblyParsingTests(TestCase):

    @patch("atomium.pdb.create_assembly_dict")
    @patch("atomium.pdb.add_assembly_info_to_mmcif")
    @patch("atomium.pdb.add_assembly_gen_to_mmcif")
    def test_can_parse_assembly(self, mock_gen, mock_info, mock_create):
        mmcif = {"mmcif": 1}
        mock_create.return_value = {"gens": [1, 2, 3]}
        parse_assembly(["line1", "line2"], mmcif)
        mock_create.assert_called_with(["line1", "line2"])
        mock_info.assert_called_with({"gens": [1, 2, 3]}, {"mmcif": 1})
        mock_gen.assert_any_call(1, {"mmcif": 1})
        mock_gen.assert_any_call(2, {"mmcif": 1})
        mock_gen.assert_any_call(3, {"mmcif": 1})



class AssemblyDictCreationTests(TestCase):

    def test_can_parse_assembly(self):
        lines = [
            "REMARK 350 BIOMOLECULE: 1                                            ",
            "REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: TRIMERIC               ",
            "REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE: TRIMERIC        ",
            "REMARK 350 SOFTWARE USED: PISA                                       ",
            "REMARK 350 TOTAL BURIED SURFACE AREA: 2660 ANGSTROM**2               ",
            "REMARK 350 SURFACE AREA OF THE COMPLEX: 10680 ANGSTROM**2            ",
            "REMARK 350 CHANGE IN SOLVENT FREE ENERGY: -11.0 KCAL/MOL             ",
            "REMARK 350 APPLY THE FOLLOWING TO CHAINS: D,                         ",
            "REMARK 350                    AND CHAINS: J, K,                      ",
            "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000       42.38700 ",
            "REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000 ",
            "REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000 ",
            "REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, B                       ",
            "REMARK 350   BIOMT1   2  1.000000  0.000000  0.000000        4.00000 ",
            "REMARK 350   BIOMT2   2  0.000000  1.000000  0.000000        5.00000 ",
            "REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        6.00000 ",
            "REMARK 350   BIOMT1   2  1.000000  0.000000  0.000000        0.00000 ",
            "REMARK 350   BIOMT2   2  0.000000  99.00000  0.000000        0.00000 ",
            "REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000 ",
        ]
        assembly = create_assembly_dict(lines)
        self.assertEqual(assembly, {
            "gens": [{
                "chains": ["D", "J", "K"], "transformations": [{
                    "matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                    "vector": [42.387, 0, 0]
                }]
            }, {
                "chains": ["A", "B"], "transformations": [{
                    "matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                    "vector": [4, 5, 6]
                }, {
                    "matrix": [[1, 0, 0], [0, 99, 0], [0, 0, 1]],
                    "vector": [0, 0, 0]
                }]
            }],
            "software": "PISA", "ABSA (A^2)": "2660",
            "SSA (A^2)": "10680", "MORE": "-11.0"
        })
    

    def test_can_handle_no_meta_data(self):
        lines = [
            "REMARK 350 BIOMOLECULE: 1                                            ",
            "REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: TRIMERIC               ",
            "REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE: TRIMERIC        ",
            "REMARK 350 APPLY THE FOLLOWING TO CHAINS: A,                         ",
            "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000       42.38700 ",
            "REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000 ",
            "REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000 ",
        ]
        assembly = create_assembly_dict(lines)
        self.assertEqual(assembly, {
            "gens": [{
                "chains": ["A"], "transformations": [{
                    "matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                    "vector": [42.387, 0, 0]
                }]
            }]
        })



class AssemblyInfoTests(TestCase):

    def test_can_handle_meta_info(self):
        mmcif = {"pdbx_struct_assembly": [], "pdbx_struct_assembly_prop": []}
        assembly = {
            "software": "PISA", "ABSA (A^2)": "2660",
            "SSA (A^2)": "10680", "MORE": "-11.0"
        }
        add_assembly_info_to_mmcif(assembly, mmcif)
        self.assertEqual(mmcif, {
            "pdbx_struct_assembly": [{
                "id": "1", "details": "?", "method_details": "PISA",
                "oligomeric_details": "?", "oligomeric_count": "?"
            }],
            "pdbx_struct_assembly_prop": [{
                "biol_id": "1", "type": "ABSA (A^2)", "value": "2660", "details": "?"
            }, {
                "biol_id": "1", "type": "MORE", "value": "-11.0", "details": "?"
            }, {
                "biol_id": "1", "type": "SSA (A^2)", "value": "10680", "details": "?"
            }]
        })
    

    def test_can_handle_no_meta_info(self):
        mmcif = {"pdbx_struct_assembly": [{}], "pdbx_struct_assembly_prop": []}
        assembly = {}
        add_assembly_info_to_mmcif(assembly, mmcif)
        self.assertEqual(mmcif, {
            "pdbx_struct_assembly": [{}, {
                "id": "2", "details": "?", "method_details": "?",
                "oligomeric_details": "?", "oligomeric_count": "?"
            }],
            "pdbx_struct_assembly_prop": []
        })



class AssemblyGenTests(TestCase):

    def test_simple_gen(self):
        gen = {"chains": ["A"], "transformations": [{
            "matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            "vector": [42.387, 0, 0]
        }]}
        mmcif = {
            "pdbx_struct_assembly": [{"id": "1"}, {"id": "2"}],
            "pdbx_struct_assembly_gen": [], "pdbx_struct_oper_list": [],
        }
        add_assembly_gen_to_mmcif(gen, mmcif)
        self.assertEqual(mmcif, {
            "pdbx_struct_assembly": [{"id": "1"}, {"id": "2"}],
            "pdbx_struct_assembly_gen": [{
                "assembly_id": "2", "oper_expression": "1", "asym_id_list": "A"
            }],
            "pdbx_struct_oper_list": [{
                "id": "1", "type": "?", "type": "?", "name": "?",
                "symmetry_operation": "x,y,z", "matrix[1][1]": "1", "matrix[1][2]": "0",
                "matrix[1][3]": "0", "vector[1]": "42.387", "matrix[2][1]": "0",
                "matrix[2][2]": "1", "matrix[2][3]": "0", "vector[2]": "0",
                "matrix[3][1]": "0", "matrix[3][2]": "0", "matrix[3][3]": "1",
                "vector[3]": "0"
            }]
        })
    

    def test_complex_gen(self):
        gen = {"chains": ["A", "B", "C"], "transformations": [{
            "matrix": [[2, 0, 0], [0, 3, 0], [0, 0, 4]],
            "vector": [1, 0, 5]
        }, {
            "matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            "vector": [42.387, 0, 0]
        }, {
            "matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            "vector": [1, 1, 1]
        }]}
        mmcif = {
            "pdbx_struct_assembly": [{"id": "1"}, {"id": "2"}],
            "pdbx_struct_assembly_gen": [],
            "pdbx_struct_oper_list": [{
                "id": "1", "type": "?", "type": "?", "name": "?",
                "symmetry_operation": "x,y,z", "matrix[1][1]": "1", "matrix[1][2]": "0",
                "matrix[1][3]": "0", "vector[1]": "42.387", "matrix[2][1]": "0",
                "matrix[2][2]": "1", "matrix[2][3]": "0", "vector[2]": "0",
                "matrix[3][1]": "0", "matrix[3][2]": "0", "matrix[3][3]": "1",
                "vector[3]": "0"
            }]
        }
        add_assembly_gen_to_mmcif(gen, mmcif)
        self.assertEqual(mmcif, {
            "pdbx_struct_assembly": [{"id": "1"}, {"id": "2"}],
            "pdbx_struct_assembly_gen": [{
                "assembly_id": "2", "oper_expression": "2,1,3", "asym_id_list": "A,B,C"
            }],
            "pdbx_struct_oper_list": [{
                "id": "1", "type": "?", "type": "?", "name": "?",
                "symmetry_operation": "x,y,z", "matrix[1][1]": "1",
                "matrix[1][2]": "0", "matrix[1][3]": "0",
                "vector[1]": "42.387", "matrix[2][1]": "0",
                "matrix[2][2]": "1", "matrix[2][3]": "0",
                "vector[2]": "0", "matrix[3][1]": "0",
                "matrix[3][2]": "0", "matrix[3][3]": "1",
                "vector[3]": "0"
            }, {
                "id": "2", "type": "?", "type": "?", "name": "?",
                "symmetry_operation": "x,y,z", "matrix[1][1]": "2",
                "matrix[1][2]": "0", "matrix[1][3]": "0",
                "vector[1]": "1", "matrix[2][1]": "0",
                "matrix[2][2]": "3", "matrix[2][3]": "0",
                "vector[2]": "0", "matrix[3][1]": "0",
                "matrix[3][2]": "0", "matrix[3][3]": "4",
                "vector[3]": "5"
            }, {
                "id": "3", "type": "?", "type": "?", "name": "?",
                "symmetry_operation": "x,y,z", "matrix[1][1]": "1",
                "matrix[1][2]": "0", "matrix[1][3]": "0",
                "vector[1]": "1", "matrix[2][1]": "0",
                "matrix[2][2]": "1", "matrix[2][3]": "0",
                "vector[2]": "1", "matrix[3][1]": "0",
                "matrix[3][2]": "0", "matrix[3][3]": "1",
                "vector[3]": "1"
            }]
        })



class Remark465ParsingTests(TestCase):

    def test_can_handle_no_remark_465(self):
        mmcif = {}
        parse_remark_465("", mmcif)
        self.assertEqual(mmcif, {})
    

    def test_can_handle_empty_remark_465(self):
        filestring = (
            "REMARK 465                                                          \n"
            "REMARK 465 MISSING  RESIDUES                                        \n"
            "REMARK 465 THE FOLLOWING  RESIDUES WERE NOT LOCATED IN THE          \n"
            "REMARK 465 EXPERIMENT.  (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN  \n"
            "REMARK 465 IDENTIFIER;  SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)   \n"
            "REMARK 465                                                          \n"
        )
        mmcif = {}
        parse_remark_465(filestring, mmcif)
        self.assertEqual(mmcif, {})


    def test_can_parse_remark_465(self):
        filestring = (
            "REMARK 465                                                          \n"             
            "REMARK 465 MISSING  RESIDUES                                        \n"             
            "REMARK 465 THE FOLLOWING  RESIDUES WERE NOT LOCATED IN THE          \n"             
            "REMARK 465 EXPERIMENT.  (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN  \n"             
            "REMARK 465 IDENTIFIER;  SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)   \n"             
            "REMARK 465                                                          \n"            
            "REMARK 465   M RES C SSSEQI                                         \n"            
            "REMARK 465     ARG A    46                                          \n"            
            "REMARK 465     GLY A    47                                          \n"            
            "REMARK 465     ALA A    48                                          \n"            
            "REMARK 465     ARG A    49                                          \n"            
            "REMARK 465     MET A    49A     "
        )
        mmcif = {}
        parse_remark_465(filestring, mmcif)
        self.assertEqual(mmcif, {"pdbx_unobs_or_zero_occ_residues": [{
            "id": "1", "PDB_model_num": "1", "polymer_flag": "Y", "occupancy_flag": "1",
            "auth_asym_id": "A", "auth_comp_id": "ARG", "auth_seq_id": "46", "PDB_ins_code": "?",
            "label_asym_id": "A", "label_comp_id": "ARG", "label_seq_id": "46",
        }, {
            "id": "2", "PDB_model_num": "1", "polymer_flag": "Y", "occupancy_flag": "1",
            "auth_asym_id": "A", "auth_comp_id": "GLY", "auth_seq_id": "47", "PDB_ins_code": "?",
            "label_asym_id": "A", "label_comp_id": "GLY", "label_seq_id": "47",
        }, {
            "id": "3", "PDB_model_num": "1", "polymer_flag": "Y", "occupancy_flag": "1",
            "auth_asym_id": "A", "auth_comp_id": "ALA", "auth_seq_id": "48", "PDB_ins_code": "?",
            "label_asym_id": "A", "label_comp_id": "ALA", "label_seq_id": "48",
        }, {
            "id": "4", "PDB_model_num": "1", "polymer_flag": "Y", "occupancy_flag": "1",
            "auth_asym_id": "A", "auth_comp_id": "ARG", "auth_seq_id": "49", "PDB_ins_code": "?",
            "label_asym_id": "A", "label_comp_id": "ARG", "label_seq_id": "49",
        }, {
            "id": "5", "PDB_model_num": "1", "polymer_flag": "Y", "occupancy_flag": "1",
            "auth_asym_id": "A", "auth_comp_id": "MET", "auth_seq_id": "49", "PDB_ins_code": "A",
            "label_asym_id": "A", "label_comp_id": "MET", "label_seq_id": "49",
        }]})



class Remark800ParsingTests(TestCase):

    def test_can_handle_no_remark_800(self):
        mmcif = {}
        parse_remark_800("", mmcif)
        self.assertEqual(mmcif, {})
    

    def test_can_handle_empty_remark_800(self):
        mmcif = {}
        filestring = (
            "REMARK 800                                                      \n"
            "REMARK 800 SITE                                                 \n"
        )
        parse_remark_800(filestring, mmcif)
        self.assertEqual(mmcif, {})


    def test_can_parse_remark_800(self):
        mmcif = {}
        filestring = (
            "REMARK 800                                                      \n"
            "REMARK 800 SITE                                                 \n"
            "REMARK 800 SITE_IDENTIFIER: AC1                                 \n"             
            "REMARK 800 EVIDENCE_CODE: SOFTWARE                              \n"             
            "REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE BU2 A 5001\n"             
            "REMARK 800 SITE_IDENTIFIER: AC2                                 \n"             
            "REMARK 800 EVIDENCE_CODE: SOFTWARE                              \n"             
            "REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE BU2 B 5002\n"
        )
        parse_remark_800(filestring, mmcif)
        self.assertEqual(mmcif, {"struct_site": [{
            "id": "AC1", "pdbx_evidence_code": "Software", "pdbx_auth_asym_id": "?",
            "pdbx_auth_comp_id": "?", "pdbx_auth_seq_id": "?", "pdbx_auth_ins_code": "?",
            "pdbx_num_residues": "?",  "details": "BINDING SITE FOR RESIDUE BU2 A 5001"
        }, {
            "id": "AC2", "pdbx_evidence_code": "Software", "pdbx_auth_asym_id": "?",
            "pdbx_auth_comp_id": "?", "pdbx_auth_seq_id": "?", "pdbx_auth_ins_code": "?",
            "pdbx_num_residues": "?",  "details": "BINDING SITE FOR RESIDUE BU2 B 5002"
        }]})



class Cryst1ParsingTests(TestCase):

    def test_can_handle_no_cryst1(self):
        mmcif = {}
        parse_cryst1("", mmcif)
        self.assertEqual(mmcif, {})
    

    def test_can_parse_cryst1(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        filestring = "HEADER\nCRYST1   52.000   58.600   61.900  90.00  94.28  80.00 P 21 21 21    8"
        parse_cryst1(filestring, mmcif)
        self.assertEqual(mmcif, {
            "entry": [{"id": "1XXX"}],
            "symmetry": [{
                "entry_id": "1XXX", "space_group_name_H-M": "P 21 21 21",
                "pdbx_full_space_group_name_H-M": "?", "cell_setting": "?",
                "Int_Tables_number": "?"
            }],
            "cell": [{
                "entry_id": "1XXX", "length_a": "52.000", "length_b": "58.600", "length_c": "61.900",
                "angle_alpha": "90.00", "angle_beta": "94.28", "angle_gamma": "80.00", "Z_pdb": "8",
                "pdbx_unique_axis": "?"
            }]
        })



class OrigxParsingTests(TestCase):

    def test_can_handle_no_origxn(self):
        mmcif = {}
        parse_origx("", mmcif)
        self.assertEqual(mmcif, {})
    

    def test_can_parse_origx(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        filestring = (
            "ORIGX1      0.963457  0.136613  0.230424       16.61000          \n"     
            "ORIGX2     -0.158977  0.983924  0.081383       13.72000          \n"     
            "ORIGX3     -0.215598 -0.115048  0.969683       37.65000          \n"
        )
        parse_origx(filestring, mmcif)
        self.assertEqual(mmcif, {
            "entry": [{"id": "1XXX"}],
            "database_PDB_matrix": [{
                "entry_id": "1XXX", "origx[1][1]": "0.963457", "origx[1][2]": "0.136613",
                "origx[1][3]": "0.230424", "origx[2][1]": "-0.158977", "origx[2][2]": "0.983924",
                "origx[2][3]": "0.081383", "origx[3][1]": "-0.215598", "origx[3][2]": "-0.115048", 
                "origx[3][3]": "0.969683", "origx_vector[1]": "16.61000", "origx_vector[2]": "13.72000",
                "origx_vector[3]": "37.65000", 
            }]
        })



class ScaleParsingTests(TestCase):

    def test_can_handle_no_scale(self):
        mmcif = {}
        parse_scale("", mmcif)
        self.assertEqual(mmcif, {})
    

    def test_can_parse_scale(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        filestring = (
            "SCALE1      0.963457  0.136613  0.230424       16.61000          \n"     
            "SCALE2     -0.158977  0.983924  0.081383       13.72000          \n"     
            "SCALE3     -0.215598 -0.115048  0.969683       37.65000          \n"
        )
        parse_scale(filestring, mmcif)
        self.assertEqual(mmcif, {
            "entry": [{"id": "1XXX"}],
            "atom_sites": [{
                "entry_id": "1XXX", "fract_transf_matrix[1][1]": "0.963457", "fract_transf_matrix[1][2]": "0.136613",
                "fract_transf_matrix[1][3]": "0.230424", "fract_transf_matrix[2][1]": "-0.158977",
                "fract_transf_matrix[2][2]": "0.983924", "fract_transf_matrix[2][3]": "0.081383",
                "fract_transf_matrix[3][1]": "-0.215598", "fract_transf_matrix[3][2]": "-0.115048", 
                "fract_transf_matrix[3][3]": "0.969683", "fract_transf_vector[1]": "16.61000",
                "fract_transf_vector[2]": "13.72000", "fract_transf_vector[3]": "37.65000", 
            }]
        })



class MtrixParsingTests(TestCase):

    def test_can_handle_no_mtrix(self):
        mmcif = {}
        parse_mtrix("", mmcif)
        self.assertEqual(mmcif, {})
    

    def test_can_parse_single_mtrix(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        filestring = (
            "MTRIX1   1 -1.000000  0.000000  0.000000        0.00000    1        \n"  
            "MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1        \n"  
            "MTRIX3   1  0.000000  0.000000 -1.000000        0.00000    1        \n"
        )
        parse_mtrix(filestring, mmcif)
        self.assertEqual(mmcif, {
            "entry": [{"id": "1XXX"}],
            "struct_ncs_oper": [{
                "id": "1", "code": "given", "details": "?",
                "matrix[1][1]": "-1.000000", "matrix[1][2]": "0.000000",
                "matrix[1][3]": "0.000000", "matrix[2][1]": "0.000000",
                "matrix[2][2]": "1.000000", "matrix[2][3]": "0.000000",
                "matrix[3][1]": "0.000000", "matrix[3][2]": "0.000000", 
                "matrix[3][3]": "-1.000000", "vector[1]": "0.00000",
                "vector[2]": "0.00000", "vector[3]": "0.00000", 
            }]
        })
    

    def test_can_parse_multiple_mtrix(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        filestring = (
            "MTRIX1   1 -1.000000  0.000000  0.000000        0.00000    1        \n"  
            "MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1        \n"  
            "MTRIX3   1  0.000000  0.000000 -1.000000        0.00000    1        \n"
            "MTRIX1   2  0.963457  0.136613  0.230424       16.61000             \n"     
            "MTRIX2   2 -0.158977  0.983924  0.081383       13.72000             \n"     
            "MTRIX3   2 -0.215598 -0.115048  0.969683       37.65000             \n"
        )
        parse_mtrix(filestring, mmcif)
        self.assertEqual(mmcif, {
            "entry": [{"id": "1XXX"}],
            "struct_ncs_oper": [{
                "id": "1", "code": "given", "details": "?",
                "matrix[1][1]": "-1.000000", "matrix[1][2]": "0.000000",
                "matrix[1][3]": "0.000000", "matrix[2][1]": "0.000000",
                "matrix[2][2]": "1.000000", "matrix[2][3]": "0.000000",
                "matrix[3][1]": "0.000000", "matrix[3][2]": "0.000000", 
                "matrix[3][3]": "-1.000000", "vector[1]": "0.00000",
                "vector[2]": "0.00000", "vector[3]": "0.00000", 
            }, {
                "id": "2", "code": "generate", "details": "?",
                "matrix[1][1]": "0.963457", "matrix[1][2]": "0.136613",
                "matrix[1][3]": "0.230424", "matrix[2][1]": "-0.158977",
                "matrix[2][2]": "0.983924", "matrix[2][3]": "0.081383",
                "matrix[3][1]": "-0.215598", "matrix[3][2]": "-0.115048", 
                "matrix[3][3]": "0.969683", "vector[1]": "16.61000",
                "vector[2]": "13.72000", "vector[3]": "37.65000", 
            }]
        })



class CompndSourceParsingTests(TestCase):

    def test_can_handle_no_compnd_or_source(self):
        self.assertEqual(parse_compnd_and_source(""), {})


    @patch("atomium.pdb.parse_entity_string")
    def test_can_parse_compnd_and_source(self, mock_pes):
        filestring = (
            "COMPND    MOL_ID: 1;\n"
            "COMPND   2 MOLECULE: FATTY ACID ACYLATED INSULIN;\n"
            "COMPND   3 MOL_ID: 2;\n"
            "COMPND   4 MOLECULE: DIFFERENT INSULIN;\n"
            "COMPND   5 CHAIN: B, D, F, H;\n"
            "COMPND   6 MOL_ID: 3;\n"
            "SOURCE    MOL_ID: 1;\n"
            "SOURCE   2 ORGANISM_SCIENTIFIC: HOMO SAPIENS;\n"
            "SOURCE   3 ORGANISM_COMMON: HUMAN;\n"
            "SOURCE   4 MOL_ID: 2;\n"
            "SOURCE   5 ORGANISM_SCIENTIFIC: RAT RATTUS;\n"
        )
        mock_pes.side_effect = [
            {"MOL_ID": "1", "PROP1": "A", "CHAIN": ["A"]},
            {"MOL_ID": "2", "PROP1": "B", "CHAIN": ["B"]},
            {"MOL_ID": "3", "PROP1": "X"},
            {"MOL_ID": "1", "PROP2": "C"},
            {"MOL_ID": "2", "PROP2": "D"},
        ]
        entities = parse_compnd_and_source(filestring)
        self.assertEqual([call[0][0] for call in mock_pes.call_args_list], [
            "MOL_ID: 1; MOLECULE: FATTY ACID ACYLATED INSULIN;",
            "MOL_ID: 2; MOLECULE: DIFFERENT INSULIN; CHAIN: B, D, F, H;",
            "MOL_ID: 3;",
            "MOL_ID: 1; ORGANISM_SCIENTIFIC: HOMO SAPIENS; ORGANISM_COMMON: HUMAN;",
            "MOL_ID: 2; ORGANISM_SCIENTIFIC: RAT RATTUS;",
        ])
        self.assertEqual(entities, {
            "1": {"id": "1", "PROP1": "A", "PROP2": "C", "CHAIN": ["A"]},
            "2": {"id": "2", "PROP1": "B", "PROP2": "D", "CHAIN": ["B"]},
        })



class EntityStringParsingTests(TestCase):

    def test_can_parse_entity_string(self):
        string = "MOL_ID: 2; MOLECULE: DIFFERENT INSULIN; CHAIN: B, D, F, H; ENGINEERED: YES; SYNTHETIC: NO;"
        self.assertEqual(parse_entity_string(string), {
            "MOL_ID": "2", "MOLECULE": "DIFFERENT INSULIN", "ENGINEERED": True,
            "CHAIN": ("B", "D", "F", "H"), "SYNTHETIC": False, "molecules": {}
        })



class DbrefParsingTests(TestCase):

    def test_can_handle_no_dbref(self):
        filestring = "HEADER\n"
        self.assertEqual(parse_dbref(filestring), {})
    

    def test_can_parse_dbref(self):
        filestring = (
            "DBREF  4OPJ A   59   196  UNP    Q9KEI9   RNH1_BACHD      58    195\n"
            "DBREF  4OPJ C   59A  196  UNP    Q9KEI9   RNH1_BACHD      58    195\n"
            "DBREF1 5MLU D   29   121  UNP                  A0A1B8Y853_XENTR\n"
            "DBREF2 5MLU D     A0A1B8Y853                         33         125\n"
            "DBREF  4OPJ A    1    12B PDB    4OPJ     4OPJ             2P    13\n"
            "DBREF  4OPJ D    1    12  PDB    4OPJ     4OPJ             2     13C\n"
        )
        self.assertEqual(parse_dbref(filestring), {
            "A": {"dbrefs": [{
                "start": "59", "start_insert": "", "end": "196", "end_insert": "", "database": "UNP",
                "accession": "Q9KEI9", "id": "RNH1_BACHD", "db_start": "58", "db_start_insert": "",
                "order": 0, "db_end": "195", "db_end_insert": ""
            }, {
                "start": "1", "start_insert": "", "end": "12", "end_insert": "B", "database": "PDB",
                "accession": "4OPJ", "id": "4OPJ", "db_start": "2", "db_start_insert": "P", "order": 4,
                "db_end": "13", "db_end_insert": ""
            }]},
            "C": {"dbrefs": [{
                "start": "59", "start_insert": "A", "end": "196", "end_insert": "", "database": "UNP",
                "accession": "Q9KEI9", "id": "RNH1_BACHD", "db_start": "58", "db_start_insert": "",
                "order": 1, "db_end": "195", "db_end_insert": ""
            }]},
            "D": {"dbrefs": [{
                "start": "29", "start_insert": "", "end": "121", "end_insert": "", "database": "UNP", 
                "accession": "A0A1B8Y853", "id": "A0A1B8Y853_XENTR", "db_start": "33", "db_start_insert": "",
                "order": 2, "db_end": "125", "db_end_insert": ""
            }, {
                "start": "1", "start_insert": "", "end": "12", "end_insert": "", "database": "PDB",
                "accession": "4OPJ", "id": "4OPJ", "db_start": "2", "db_start_insert": "", "order": 5,
                "db_end": "13", "db_end_insert": "C"
            }]}
        })



class SeqadvParsingTests(TestCase):

    def test_can_handle_no_seqadv(self):
        filestring = "HEADER\nDBREF1\nDBREF2"
        polymers = {"A": {}}
        parse_seqadv(filestring, polymers)
        self.assertEqual(polymers, {"A": {}})
    

    def test_can_parse_seqadv(self):
        filestring = (
            "SEQADV 4OPJ GLY A   55A UNP  Q9KEI9              EXPRESSION TAG\n"
            "SEQADV 4OPJ ASN A  132  UNP  Q9KEI9    ASP   132 ENGINEERED MUTATION\n"
            "SEQADV 4OPJ MET C   58  UNP  Q9KEI9              EXPRESSION TAG\n"
            "SEQADV 4OPJ ASN C  132  UNP  Q9KEI9    ASP   132 ENGINEERED MUTATION\n"
        )
        polymers = {"A": {"dbrefs": [1, 2]}}
        parse_seqadv(filestring, polymers)
        self.assertEqual(polymers, {
            "A": {"dbrefs": [1, 2], "differences": [{
                "name": "GLY", "number": "55", "insert": "A", "database": "UNP", "accession": "Q9KEI9",
                "db_name": "", "db_number": "", "comment": "EXPRESSION TAG"
            }, {
                "name": "ASN", "number": "132", "insert": "", "database": "UNP", "accession": "Q9KEI9",
                "db_name": "ASP", "db_number": "132", "comment": "ENGINEERED MUTATION"
            }]},
            "C": {"differences": [{
                "name": "MET", "number": "58", "insert": "", "database": "UNP", "accession": "Q9KEI9",
                "db_name": "", "db_number": "", "comment": "EXPRESSION TAG"
            }, {
                "name": "ASN", "number": "132", "insert": "", "database": "UNP", "accession": "Q9KEI9",
                "db_name": "ASP", "db_number": "132", "comment": "ENGINEERED MUTATION"
            }]}
        })



class SeqresParsingTests(TestCase):

    def test_can_handle_no_seqres(self):
        filestring = "HEADER\nDBREF1\nDBREF2"
        polymers = {"A": {}}
        parse_seqres(filestring, polymers)
        self.assertEqual(polymers, {"A": {}})
    

    def test_can_parse_seqres(self):
        filestring = (
            "SEQRES   1 A   10   DC  DC  DT  DC  DT  DA  DG  DA  DG  DG   \n"
            "SEQRES   1 B   10   DC  DC  DT  DC  DT  DA  DG  DA  DG  DG   \n"
            "SEQRES   1 C  229  LEU ARG SER ARG ARG VAL ASP VAL MET ASP VAL MET ASN\n"
            "SEQRES   2 C  229  ARG LEU ILE LEU ALA MET ASP LEU MET ASN ARG ASP ASP\n"
        )
        polymers = {"A": {"dbrefs": [1, 2]}}
        parse_seqres(filestring, polymers)
        self.assertEqual(polymers, {
            "A": {"dbrefs": [1, 2], "residues": ["DC", "DC", "DT", "DC", "DT", "DA", "DG", "DA", "DG", "DG"]},
            "B": {"residues": ["DC", "DC", "DT", "DC", "DT", "DA", "DG", "DA", "DG", "DG"]},
            "C": {"residues": [
                "LEU", "ARG", "SER", "ARG", "ARG", "VAL", "ASP", "VAL", "MET", "ASP", "VAL", "MET", "ASN",
                "ARG", "LEU", "ILE", "LEU", "ALA", "MET", "ASP", "LEU", "MET", "ASN", "ARG", "ASP", "ASP"
            ]}
        })



class ModresParsingTests(TestCase):

    def test_can_handle_no_modres(self):
        filestring = "HEADER\nDBREF1\nDBREF2"
        polymers = {"A": {}}
        parse_modres(filestring, polymers)
        self.assertEqual(polymers, {"A": {}})
    

    def test_can_parse_modres(self):
        filestring = (
            "MODRES 2F1M MSE A  223  MET  SELENOMETHIONINE     \n"                              
            "MODRES 2F1M MSE A  287A MET  HYDROMETH     \n"                              
            "MODRES 2F1M MSE A  288  MET  SELENOMETHIONINE     \n"                              
            "MODRES 2F1M MSE A  291  MET  MEGAMETH     \n"                              
            "MODRES 2F1M MSE B  287  MET  SELENOMETHIONINE     \n"                              
            "MODRES 2F1M MSE B  288  MET  FERROMETH     \n"                                                          
            "MODRES 2F1M MSE C  223  PRO  SELENOMETHIONINE     \n"                              
            "MODRES 2F1M MSE C  224  HIS  SELENOMETHIONINE     \n"
        )
        polymers = {"A": {"dbrefs": [1, 2]}}
        parse_modres(filestring, polymers)
        self.assertEqual(polymers, {
            "A": {"dbrefs": [1, 2], "modified": [
                {"name": "MSE", "number": "223", "insert": "", "standard_name": "MET", "comment": "SELENOMETHIONINE"},
                {"name": "MSE", "number": "287", "insert": "A", "standard_name": "MET", "comment": "HYDROMETH"},
                {"name": "MSE", "number": "288", "insert": "", "standard_name": "MET", "comment": "SELENOMETHIONINE"},
                {"name": "MSE", "number": "291", "insert": "", "standard_name": "MET", "comment": "MEGAMETH"}
            ]},
            "B": {"modified": [
                {"name": "MSE", "number": "287", "insert": "", "standard_name": "MET", "comment": "SELENOMETHIONINE"},
                {"name": "MSE", "number": "288", "insert": "", "standard_name": "MET", "comment": "FERROMETH"}
            ]},
            "C": {"modified": [
                {"name": "MSE", "number": "223", "insert": "", "standard_name": "PRO", "comment": "SELENOMETHIONINE"},
                {"name": "MSE", "number": "224", "insert": "", "standard_name": "HIS", "comment": "SELENOMETHIONINE"}
            ]}
        })



class HetParsingTests(TestCase):

    def test_can_handle_no_het(self):
        entities, molecules = parse_het("")
        self.assertEqual(entities, {})
        self.assertEqual(molecules, {})
    

    def test_can_parse_het(self):
        filestring = (
            "HET    BU2  A5001       6         \n"
            "HET    BU2  B5002       6         \n"
            "HET    XMP  A2001      24         \n"
            "HET    XMP  B2002      24\n"
        )
        entities, molecules = parse_het(filestring)
        self.assertEqual(entities, {"BU2": {}, "XMP": {}})
        self.assertEqual(molecules, {
            ("A", "BU2", "5001", ""): {},
            ("A", "XMP", "2001", ""): {},
            ("B", "BU2", "5002", ""): {},
            ("B", "XMP", "2002", ""): {},
        })



class HetnamParsingTests(TestCase):

    def test_can_handle_no_hetnam(self):
        entities = {"ABC": {}}
        parse_hetnam("", entities)
        self.assertEqual(entities, {"ABC": {}})
    

    def test_can_parse_hetnam(self):
        entities = {"BU2": {}}
        filestring = (
            "HETNAM     BU2 1,3-BUTANEDIOL              \n"                                     
            "HETNAM     XMP XANTHOSINE-5'-MONOPHOSPHATE   "
        )
        parse_hetnam(filestring, entities)
        self.assertEqual(entities, {
            "BU2": {"name": "1,3-BUTANEDIOL"}, "XMP": {"name": "XANTHOSINE-5'-MONOPHOSPHATE"}
        })
    

    def test_can_parse_multiline_hetnam(self):
        entities = {"BU2": {}}
        filestring = (
            "HETNAM     BU2 1,3-BUTANEDIOL              \n"                                     
            "HETNAM     3N6 N-{(1S)-5-AMINO-1-[(4-PYRIDIN-4-YLPIPERAZIN-1-YL)      \n"
            "HETNAM   2 3N6  CARBONYL]PENTYL}-3,5-DIBROMO-NALPHA-{[4-(2-OXO-1,4-   \n"
            "HETNAM   3 3N6  DIHYDROQUINAZOLIN-3 (2H)-YL)PIPERIDIN-1-YL]CARBONYL}- \n"
            "HETNAM   4 3N6  D-TYROSINAMIDE "
        )
        parse_hetnam(filestring, entities)
        self.assertEqual(entities, {
            "BU2": {"name": "1,3-BUTANEDIOL"},
            "3N6": {"name": "N-{(1S)-5-AMINO-1-[(4-PYRIDIN-4-YLPIPERAZIN-1-YL) CARBONYL]PENTYL}-3,5-DIBROMO-NALPHA-{[4-(2-OXO-1,4- DIHYDROQUINAZOLIN-3 (2H)-YL)PIPERIDIN-1-YL]CARBONYL}- D-TYROSINAMIDE"}
        })



class HetsynParsingTests(TestCase):

    def test_can_handle_no_hetsyn(self):
        entities = {"ABC": {}}
        parse_hetsyn("", entities)
        self.assertEqual(entities, {"ABC": {}})
    

    def test_can_parse_hetsyn(self):
        entities = {"BU2": {}}
        filestring = (
            "HETSYN     BU2 BTYLO             \n"                                     
            "HETSYN     XMP XANAX;XANAX2   "
        )
        parse_hetsyn(filestring, entities)
        self.assertEqual(entities, {
            "BU2": {"synonyms": ["BTYLO"]}, "XMP": {"synonyms": ["XANAX", "XANAX2"]}
        })
    

    def test_can_parse_multiline_hetsyn(self):
        entities = {"BU2": {}}
        filestring = (
            "HETSYN     BU2 1,3-BUTANEDIOL              \n"                                     
            "HETSYN     3N6 N-{(1S)-5-AMINO-1-[(4-PYRIDIN-4-YLPIPERAZIN-1-YL)      \n"
            "HETSYN   2 3N6  CARBONYL]PENTYL}-3,5-DIBROMO-NALPHA-{[4-(2-OXO-1,4-   \n"
            "HETSYN   3 3N6  DIHYDROQUINAZOLIN-3 (2H)-YL)PIPERIDIN-1-YL]CARBONYL}- \n"
            "HETSYN   4 3N6  D-TYROSINAMIDE "
        )
        parse_hetsyn(filestring, entities)
        self.assertEqual(entities, {
            "BU2": {"synonyms": ["1,3-BUTANEDIOL"]},
            "3N6": {"synonyms": ["N-{(1S)-5-AMINO-1-[(4-PYRIDIN-4-YLPIPERAZIN-1-YL) CARBONYL]PENTYL}-3,5-DIBROMO-NALPHA-{[4-(2-OXO-1,4- DIHYDROQUINAZOLIN-3 (2H)-YL)PIPERIDIN-1-YL]CARBONYL}- D-TYROSINAMIDE"]}
        })



class FormulParsingTests(TestCase):

    def test_can_handle_no_formul(self):
        entities = {"ABC": {}}
        parse_formul("", entities)
        self.assertEqual(entities, {"ABC": {}})
    

    def test_can_parse_formul(self):
        entities = {"BU2": {}}
        filestring = (
            "FORMUL   3  BU2    2(C4 H10 O2)                 \n"                                 
            "FORMUL   5  XMP    2(C10 H14 N4 O9 P 1+)        \n"                                 
            "FORMUL   7  HOH   *180(H2 O) \n"
        )
        parse_formul(filestring, entities)
        self.assertEqual(entities, {
            "BU2": {"formula": "2(C4 H10 O2)", "is_water": False},
            "XMP": {"formula": "2(C10 H14 N4 O9 P 1+)", "is_water": False},
            "HOH": {"formula": "180(H2 O)", "is_water": True},
        })



class EntitiesFromAtomsTests(TestCase):

    def test_no_atoms(self):
        polymer_entities, polymers, nonpoly_entities, nonpolies = {}, {}, {}, {}
        update_entities_from_atoms("", polymer_entities, polymers, nonpoly_entities, nonpolies)
        self.assertEqual(polymer_entities, {})
        self.assertEqual(polymers, {})
        self.assertEqual(nonpoly_entities, {})
        self.assertEqual(nonpolies, {})
    

    def test_can_generate_from_scratch(self):
        filestring = (
            "ATOM      1  N   VAL A  11    \n"
            "ATOM      2  CA  VAL A  11    \n"
            "ATOM      9  CA  MET A  12    \n"
            "TER    1558      ILE A 222\n"
            "ATOM   1559  N   PRO B1011    \n"
            "ATOM   1567  CA  CYS B1012    \n"
            "TER    3193      GLU B1229\n"
            "HETATM 3194  C1  BU2 A5001A   \n"
            "HETATM 3195  O1  BU2 A5001A   \n"
            "HETATM 3200  P   XMP A2001    \n"  
            "HETATM 3201  O1P XMP A2001    \n"
            "HETATM 3224  C1  BU2 B5002    \n"
            "HETATM 3230  P   XMP B2002    \n"
            "HETATM 3254  O   HOH A3005    \n"
            "HETATM 3350  O   HOH B3001    \n"
        )
        polymer_entities, polymers, nonpoly_entities, nonpolies = {}, {}, {}, {}
        update_entities_from_atoms(filestring, polymer_entities, polymers, nonpoly_entities, nonpolies)
        self.assertEqual(polymer_entities, {"1": {"CHAIN": ("A",)}, "2": {"CHAIN": ("B",)}})
        self.assertEqual(polymers, {
            "A": {"observed_residues": ["VAL", "MET"]},
            "B": {"observed_residues": ["PRO", "CYS"]},
        })
        self.assertEqual(nonpoly_entities, {"BU2": {}, "XMP": {}, "HOH": {}})
        self.assertEqual(nonpolies, {
            ("A", "BU2", "5001", "A"): {}, ("A", "XMP", "2001", ""): {},
            ("B", "BU2", "5002", ""): {}, ("B", "XMP", "2002", ""): {},
            ("A", "HOH", "3005", ""): {}, ("B", "HOH", "3001", ""): {},
        })
    

    def test_can_update_existing(self):
        filestring = (
            "ATOM      1  N   VAL A  11    \n"
            "ATOM      2  CA  VAL A  11    \n"
            "ATOM      9  CA  MET A  12    \n"
            "TER    1558      ILE A 222\n"
            "ATOM   1559  N   PRO B1011    \n"
            "ATOM   1567  CA  CYS B1012    \n"
            "TER    3193      GLU B1229\n"
            "ATOM   1559  N   PRO C1511    \n"
            "ATOM   1567  CA  CYS C1512    \n"
            "TER    3193      GLU C1529\n"
            "HETATM 3194  C1  BU2 A5001A   \n"
            "HETATM 3195  O1  BU2 A5001A   \n"
            "HETATM 3200  P   XMP A2001    \n"  
            "HETATM 3201  O1P XMP A2001    \n"
            "HETATM 3224  C1  BU2 B5002    \n"
            "HETATM 3230  P   XMP B2002    \n"
            "HETATM 3254  O   HOH A3005    \n"
            "HETATM 3350  O   HOH B3001    \n"
        )
        polymer_entities, polymers, nonpoly_entities, nonpolies = {}, {}, {}, {}
        polymer_entities = {"1": {"CHAIN": ("A", "B")}}
        polymers = {"A": {"name": "AA"}, "B": {"name": "BB"}}
        nonpoly_entities = {"BU2": {"name": "BB"}}
        nonpolies = {("A", "BU2", "5001", "A"): {}}
        update_entities_from_atoms(filestring, polymer_entities, polymers, nonpoly_entities, nonpolies)
        self.assertEqual(polymer_entities, {"1": {"CHAIN": ("A", "B")}, "2": {"CHAIN": ("C",)}})
        self.assertEqual(polymers, {
            "A": {"observed_residues": ["VAL", "MET"], "name": "AA"},
            "B": {"observed_residues": ["PRO", "CYS"], "name": "BB"},
            "C": {"observed_residues": ["PRO", "CYS"]},
        })
        self.assertEqual(nonpoly_entities, {"BU2": {"name": "BB"}, "XMP": {}, "HOH": {}})
        self.assertEqual(nonpolies, {
            ("A", "BU2", "5001", "A"): {}, ("A", "XMP", "2001", ""): {},
            ("B", "BU2", "5002", ""): {}, ("B", "XMP", "2002", ""): {},
            ("A", "HOH", "3005", ""): {}, ("B", "HOH", "3001", ""): {},
        })



class EntityFinalizationTests(TestCase):

    def test_can_fully_populate_entities(self):
        polymers = {"1": {}, "2": {}}
        non_polymers = {"XMP": {}, "BU2": {}}
        finalize_entities(polymers, non_polymers)
        self.assertEqual(polymers, {
            "1": {"id": "1", "molecules": {}}, "2": {"id": "2", "molecules": {}}
        })
        self.assertEqual(non_polymers, {
            "XMP": {"id": "3", "name": "", "synonyms": [], "formula": "", "is_water": False, "molecules": {}},
            "BU2": {"id": "4", "name": "", "synonyms": [], "formula": "", "is_water": False, "molecules": {}},
        })
    

    def test_can_extend_populated_entities(self):
        polymers = {"1": {"id": "2"}, "3": {"id": "4", "molecules": {1: 2}}}
        non_polymers = {
            "XMP": {"is_water": True},
            "BU2": {"name": "B", "synonyms": ["S"], "formula": "B2"},
            "HOH": {"is_water": False}
        }
        finalize_entities(polymers, non_polymers)
        self.assertEqual(polymers, {
            "1": {"id": "1", "molecules": {}}, "2": {"id": "2", "molecules": {1: 2}}
        })
        self.assertEqual(non_polymers, {
            "XMP": {"id": "3", "name": "", "synonyms": [], "formula": "", "is_water": True, "molecules": {}},
            "BU2": {"id": "4", "name": "B", "synonyms": ["S"], "formula": "B2", "is_water": False, "molecules": {}},
            "HOH": {"id": "5", "name": "", "synonyms": [], "formula": "", "is_water": False, "molecules": {}},
        })
    

    def test_can_determine_water(self):
        polymers = {"1": {}, "2": {}}
        non_polymers = {"HOH": {}, "WAT": {}}
        finalize_entities(polymers, non_polymers)
        self.assertEqual(polymers, {
            "1": {"id": "1", "molecules": {}}, "2": {"id": "2", "molecules": {}}
        })
        self.assertEqual(non_polymers, {
            "HOH": {"id": "3", "name": "", "synonyms": [], "formula": "", "is_water": True, "molecules": {}},
            "WAT": {"id": "4", "name": "", "synonyms": [], "formula": "", "is_water": True, "molecules": {}},
        })



class PolymerFinalizationTests(TestCase):

    @patch("atomium.pdb.get_alignment_indices")
    def test_can_generate_all_fields(self, mock_align):
        mock_align.side_effect = ["align1", "align2"]
        polymers = {"A": {}, "B": {}}
        finalize_polymers(polymers)
        self.assertEqual(polymers, {
            "A": {"dbrefs": [], "differences": [], "modified": [], "residues": [], "observed_residues": [], "alignment": "align1"},
            "B": {"dbrefs": [], "differences": [], "modified": [], "residues": [], "observed_residues": [], "alignment": "align2"},
        })
        self.assertEqual(
            [c[0] for c in mock_align.call_args_list],
            [([], []), ([], [])]
        )
    

    @patch("atomium.pdb.get_alignment_indices")
    def test_can_update_existing_fields(self, mock_align):
        mock_align.side_effect = ["align1", "align2"]
        polymers = {
            "A": {"dbrefs": [1], "observed_residues": ["P", "C"]},
            "B": {"differences": [3], "modified": [5], "residues": ["A", "T"]}
        }
        finalize_polymers(polymers)
        self.assertEqual(polymers, {
            "A": {"dbrefs": [1], "differences": [], "modified": [], "residues": ["P", "C"], "observed_residues": ["P", "C"], "alignment": "align1"},
            "B": {"dbrefs": [], "differences": [3], "modified": [5], "residues": ["A", "T"], "observed_residues": [], "alignment": "align2"},
        })
        self.assertEqual(
            [c[0] for c in mock_align.call_args_list],
            [(["P", "C"], ["P", "C"]), (["A", "T"], [])]
        )



class MoleculeToEntitiesTests(TestCase):

    def test_can_add_molecules_to_entities(self):
        polymer_entities, polymers, nonpoly_entities, nonpolies = {}, {}, {}, {}
        polymer_entities = {
            "1": {"CHAIN": ("A", "B"), "molecules": {}},
            "2": {"CHAIN": ("C"), "molecules": {}},
        }
        polymers = {"A": 1, "B": 2, "C": 3}
        nonpoly_entities = {"XMP": {"molecules": {}}, "BU2": {"molecules": {}}}
        nonpolies = {
            ("A", "XMP", "100", "A"): {}, ("B", "XMP", "101", ""): {},
            ("A", "BU2", "1000", ""): {}
        }
        add_molecules_to_entities(polymer_entities, polymers, nonpoly_entities, nonpolies)
        self.assertEqual(polymer_entities, {
            "1": {"CHAIN": ("A", "B"), "molecules": {"A": 1, "B": 2}},
            "2": {"CHAIN": ("C"), "molecules": {"C": 3}},
        })
        self.assertEqual(nonpoly_entities, {
            "XMP": {"molecules": {("A", "XMP", "100", "A"): {}, ("B", "XMP", "101", ""): {}}},
            "BU2": {"molecules": {("A", "BU2", "1000", ""): {}}}
        })



class EntityBuildingTests(TestCase):

    def test_can_handle_no_entities(self):
        polymers = {}
        non_polymers = {}
        mmcif = {1: 2}
        build_entity_category(polymers, non_polymers, mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    @patch("atomium.pdb.update_polymer_entity")
    @patch("atomium.pdb.update_non_polymer_entity")
    def test_can_parse_entities(self, mock_nonpoly, mock_poly):
        polymers = {"1": {"id": "1"}, "2": {"id": "3"}}
        non_polymers = {
            "ABC": {"A": "A", "molecules": [1, 2]},
            "XYZ": {"X": "X", "molecules": [1]},
            "ZZZ": {"Z": "Z", "molecules": []},
        }
        mmcif = {1: 2}
        build_entity_category(polymers, non_polymers, mmcif)
        self.assertEqual(mmcif, {1: 2, "entity": [{
            "id": "1", "type": "?", "src_method": "?", "pdbx_description": "?",
            "formula_weight": "?", "pdbx_number_of_molecules": "1", "pdbx_ec": "?",
            "pdbx_mutation": "?", "pdbx_fragment": "?", "details": "?",
        }, {
            "id": "2", "type": "?", "src_method": "?", "pdbx_description": "?",
            "formula_weight": "?", "pdbx_number_of_molecules": "1", "pdbx_ec": "?",
            "pdbx_mutation": "?", "pdbx_fragment": "?", "details": "?",
        }, {
            "id": "3", "type": "?", "src_method": "?", "pdbx_description": "?",
            "formula_weight": "?", "pdbx_number_of_molecules": "1", "pdbx_ec": "?",
            "pdbx_mutation": "?", "pdbx_fragment": "?", "details": "?",
        }, {
            "id": "4", "type": "?", "src_method": "?", "pdbx_description": "?",
            "formula_weight": "?", "pdbx_number_of_molecules": "1", "pdbx_ec": "?",
            "pdbx_mutation": "?", "pdbx_fragment": "?", "details": "?",
        }]})
        mock_poly.assert_any_call(mmcif["entity"][0], {"id": "1"})
        mock_poly.assert_any_call(mmcif["entity"][1], {"id": "2"})
        mock_nonpoly.assert_any_call(mmcif["entity"][2], {"A": "A", "molecules": [1, 2], "id": "3"})
        mock_nonpoly.assert_any_call(mmcif["entity"][3], {"X": "X", "molecules": [1], "id": "4"})



class PolymerEntityUpdatingTests(TestCase):

    def test_can_update_minimal_info(self):
        entity = {
            "id": "2", "type": "?", "src_method": "?", "pdbx_description": "?",
            "formula_weight": "?", "pdbx_number_of_molecules": "1", "pdbx_ec": "?",
            "pdbx_mutation": "?", "pdbx_fragment": "?", "details": "?",
        }
        info = {"CHAIN": ("A", "B", "C")}
        update_polymer_entity(entity, info)
        self.assertEqual(entity, {
            "id": "2", "type": "polymer", "src_method": "nat", "pdbx_description": "?",
            "formula_weight": "?", "pdbx_number_of_molecules": "3", "pdbx_ec": "?",
            "pdbx_mutation": "?", "pdbx_fragment": "?", "details": "?",
        })
    

    def test_can_update_full_info(self):
        entity = {
            "id": "2", "type": "?", "src_method": "?", "pdbx_description": "?",
            "formula_weight": "?", "pdbx_number_of_molecules": "1", "pdbx_ec": "?",
            "pdbx_mutation": "?", "pdbx_fragment": "?", "details": "?",
        }
        info = {
            "CHAIN": ("A", "B", "C"), "MOLECULE": "fullname", "EC": "1.2",
            "FRAGMENT": "frag", "OTHER_DETAILS": "extra", "SYNTHETIC": "yes"
        }
        update_polymer_entity(entity, info)
        self.assertEqual(entity, {
            "id": "2", "type": "polymer", "src_method": "syn", "pdbx_description": "fullname",
            "formula_weight": "?", "pdbx_number_of_molecules": "3", "pdbx_ec": "1.2",
            "pdbx_mutation": "?", "pdbx_fragment": "frag", "details": "extra",
        })
    

    def test_can_update_full_info_with_engineered(self):
        entity = {
            "id": "2", "type": "?", "src_method": "?", "pdbx_description": "?",
            "formula_weight": "?", "pdbx_number_of_molecules": "1", "pdbx_ec": "?",
            "pdbx_mutation": "?", "pdbx_fragment": "?", "details": "?",
        }
        info = {
            "CHAIN": ("A", "B", "C"), "MOLECULE": "fullname", "EC": "1.2",
            "FRAGMENT": "frag", "OTHER_DETAILS": "extra", "ENGINEERED": "yes"
        }
        update_polymer_entity(entity, info)
        self.assertEqual(entity, {
            "id": "2", "type": "polymer", "src_method": "man", "pdbx_description": "fullname",
            "formula_weight": "?", "pdbx_number_of_molecules": "3", "pdbx_ec": "1.2",
            "pdbx_mutation": "?", "pdbx_fragment": "frag", "details": "extra",
        })



class NonPolymerEntityUpdatingTests(TestCase):

    def test_can_update_minimal_info(self):
        entity = {
            "id": "2", "type": "?", "src_method": "?", "pdbx_description": "?",
            "formula_weight": "?", "pdbx_number_of_molecules": "1", "pdbx_ec": "?",
            "pdbx_mutation": "?", "pdbx_fragment": "?", "details": "?",
        }
        info = {"molecules": {1: 2, 3: 4}, "is_water": False}
        update_non_polymer_entity(entity, info)
        self.assertEqual(entity, {
            "id": "2", "type": "non-polymer", "src_method": "syn", "pdbx_description": "?",
            "formula_weight": "?", "pdbx_number_of_molecules": "2", "pdbx_ec": "?",
            "pdbx_mutation": "?", "pdbx_fragment": "?", "details": "?",
        })
    

    def test_can_update_full_info(self):
        entity = {
            "id": "2", "type": "?", "src_method": "?", "pdbx_description": "?",
            "formula_weight": "?", "pdbx_number_of_molecules": "1", "pdbx_ec": "?",
            "pdbx_mutation": "?", "pdbx_fragment": "?", "details": "?",
        }
        info = {"molecules": {1: 2, 3: 4}, "is_water": False, "name": "fullium"}
        update_non_polymer_entity(entity, info)
        self.assertEqual(entity, {
            "id": "2", "type": "non-polymer", "src_method": "syn", "pdbx_description": "fullium",
            "formula_weight": "?", "pdbx_number_of_molecules": "2", "pdbx_ec": "?",
            "pdbx_mutation": "?", "pdbx_fragment": "?", "details": "?",
        })
    

    def test_can_update_full_info_as_water(self):
        entity = {
            "id": "2", "type": "?", "src_method": "?", "pdbx_description": "?",
            "formula_weight": "?", "pdbx_number_of_molecules": "1", "pdbx_ec": "?",
            "pdbx_mutation": "?", "pdbx_fragment": "?", "details": "?",
        }
        info = {"molecules": {1: 2, 3: 4}, "is_water": True, "name": ""}
        update_non_polymer_entity(entity, info)
        self.assertEqual(entity, {
            "id": "2", "type": "water", "src_method": "nat", "pdbx_description": "water",
            "formula_weight": "?", "pdbx_number_of_molecules": "2", "pdbx_ec": "?",
            "pdbx_mutation": "?", "pdbx_fragment": "?", "details": "?",
        })



class EntityNameComTests(TestCase):

    def test_can_not_build_category(self):
        polymers = {"1": {}, "2": {"SYNONYM": ""}}
        mmcif = {1: 2}
        build_entity_name_com(polymers, mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_not_build_category(self):
        polymers = {
            "1": {"id": "1"}, "2": {"SYNONYM": "", "id": "2"},
            "3": {"SYNONYM": "name1", "id": "3"}, "4": {"SYNONYM": "name2", "id": "4"}
        }
        mmcif = {1: 2}
        build_entity_name_com(polymers, mmcif)
        self.assertEqual(mmcif, {1: 2, "entity_name_com": [
            {"entity_id": "3", "name": "name1"}, {"entity_id": "4", "name": "name2"}
        ]})



class EntityPolyTests(TestCase):

    def test_can_handle_no_polymers(self):
        mmcif = {1: 2}
        build_entity_poly({}, mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_handle_no_valid_polymers(self):
        mmcif = {1: 2}
        build_entity_poly({"1": {"molecules": {}}, "2": {"molecules": {}}}, mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    @patch("atomium.pdb.get_sequence_strings")
    def test_can_build_category(self, mock_strings):
        mmcif = {1: 2}
        mock_strings.side_effect = [["seq1", "can1"], ["seq2", "can2"]]
        entities = {
            "1": {"id": "1", "molecules": {
                "A": {"residues": ["HIS", "GLY"], "modified": []},
                "B": {"residues": ["HIS", "GLY"], "modified": [1, 2]},
            }},
            "2": {"id": "2", "molecules": {
                "C": {"residues": ["A", "T"], "modified": [1, 2]},
            }}
        }
        build_entity_poly(entities, mmcif)
        self.assertEqual(mmcif, {1: 2, "entity_poly": [{
            "entity_id": "1",
            "type": "polypeptide(L)",
            "nstd_linkage": "no",
            "nstd_monomer": "no",
            "pdbx_seq_one_letter_code": "seq1",
            "pdbx_seq_one_letter_code_can": "can1",
            "pdbx_strand_id": "A,B",
            "pdbx_target_identifier": "?"
        }, {
            "entity_id": "2",
            "type": "polyribonucleotide",
            "nstd_linkage": "no",
            "nstd_monomer": "yes",
            "pdbx_seq_one_letter_code": "seq2",
            "pdbx_seq_one_letter_code_can": "can2",
            "pdbx_strand_id": "C",
            "pdbx_target_identifier": "?"
        }]})



class SequenceStringTests(TestCase):

    def test_can_get_protein_sequence(self):
        polymer = {"residues": ["HIS", "MET", "CYS"], "modified": []}
        self.assertEqual(get_sequence_strings(polymer), ("HMC", "HMC"))
    

    def test_can_get_nucleotide_sequence(self):
        polymer = {"residues": ["A", "U", "G"], "modified": []}
        self.assertEqual(get_sequence_strings(polymer), ("AUG", "AUG"))
    

    def test_can_get_unknown_sequence(self):
        polymer = {"residues": ["XYZ", "ABC", "UNK"], "modified": []}
        self.assertEqual(get_sequence_strings(polymer), ("(XYZ)(ABC)(UNK)", "XXX"))
    

    def test_can_get_unknown_sequence_with_modified_lookup(self):
        polymer = {"residues": ["XYZ", "ABC", "UNK"], "modified": [{"standard_name": "HIS", "name": "ABC"}]}
        self.assertEqual(get_sequence_strings(polymer), ("(XYZ)(ABC)(UNK)", "XHX"))



class EntityPolySeqTests(TestCase):

    def test_can_handle_no_polymers(self):
        mmcif = {1: 2}
        build_entity_poly_seq({}, mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_handle_no_valid_polymers(self):
        mmcif = {1: 2}
        build_entity_poly_seq({"1": {"molecules": {}}, "2": {"molecules": {}}}, mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_not_build_category(self):
        mmcif = {1: 2}
        entities = {
            "1": {"id": "1", "molecules": {
                "A": {"residues": ["HIS", "GLY"], "modified": []},
                "B": {"residues": ["HIS", "GLY"], "modified": [1, 2]},
            }},
            "2": {"id": "2", "molecules": {
                "C": {"residues": ["A", "T"], "modified": [1, 2]},
            }}
        }
        build_entity_poly_seq(entities, mmcif)
        self.assertEqual(mmcif, {1: 2, "entity_poly_seq": [
            {"entity_id": "1", "num": "1", "mon_id": "HIS", "hetero": "n"},
            {"entity_id": "1", "num": "2", "mon_id": "GLY", "hetero": "n"},
            {"entity_id": "2", "num": "1", "mon_id": "A", "hetero": "n"},
            {"entity_id": "2", "num": "2", "mon_id": "T", "hetero": "n"},
        ]})



class StructRefTests(TestCase):

    def test_can_handle_no_entities(self):
        mmcif = {1: 2}
        build_struct_ref({}, mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_handle_no_valid_polymers(self):
        mmcif = {1: 2}
        build_struct_ref({"1": {"molecules": {}}, "2": {"molecules": {}}}, mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_handle_no_dbrefs(self):
        mmcif = {1: 2}
        build_struct_ref({
            "1": {"molecules": {"A": {"dbrefs": []}}},
            "2": {"molecules": {"B": {"dbrefs": []}}}
        }, mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_build_table(self):
        mmcif = {1: 2}
        build_struct_ref({
            "1": {"id": "1", "molecules": {
                "A": {"dbrefs": [
                    {"id": "1", "database": "UNP", "start": "2", "accession": "ABC"}
                ]},
                "B": {"dbrefs": [
                    {"id": "1", "database": "UNP", "start": "3", "accession": "ABC"}
                ]}
            }},
            "2": {"id": "2", "molecules": {
                "C": {"dbrefs": [
                    {"id": "1", "database": "UNP", "start": "4", "accession": "YYY"}
                ]},
                "D": {"dbrefs": [
                    {"id": "1", "database": "RFS", "start": "5", "accession": "ZZZ"},
                    {"id": "2", "database": "RFS", "start": "6", "accession": "ZZZ"},
                ]}
            }},
        }, mmcif)
        self.assertEqual(mmcif, {1: 2, "struct_ref": [{
            "id": "1", "db_name": "UNP", "db_code": "1", "entity_id": "1", "pdbx_seq_one_letter_code": "?",
            "pdbx_align_begin": "?", "pdbx_db_accession": "ABC", "pdbx_db_isoform": "?"
        }, {
            "id": "2", "db_name": "UNP", "db_code": "1", "entity_id": "2", "pdbx_seq_one_letter_code": "?",
            "pdbx_align_begin": "?", "pdbx_db_accession": "YYY", "pdbx_db_isoform": "?"
        }, {
            "id": "3", "db_name": "RFS", "db_code": "2", "entity_id": "2", "pdbx_seq_one_letter_code": "?",
            "pdbx_align_begin": "?", "pdbx_db_accession": "ZZZ", "pdbx_db_isoform": "?"
        }]})
    

    def test_can_build_empty_table(self):
        mmcif = {1: 2}
        build_struct_ref({
            "1": {"id": "1", "molecules": {
                "A": {"dbrefs": [
                    {"id": "", "database": "", "start": "", "accession": ""}
                ]}
            }},
        }, mmcif)
        self.assertEqual(mmcif, {1: 2, "struct_ref": [{
            "id": "1", "db_name": "?", "db_code": "?", "entity_id": "1", "pdbx_seq_one_letter_code": "?",
            "pdbx_align_begin": "?", "pdbx_db_accession": "?", "pdbx_db_isoform": "?"
        }]})



class StructRefSeqTests(TestCase):

    def test_can_handle_no_entities(self):
        mmcif = {1: 2}
        build_struct_ref_seq({}, mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_handle_no_valid_polymers(self):
        mmcif = {1: 2}
        build_struct_ref_seq({"1": {"molecules": {}}, "2": {"molecules": {}}}, mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_handle_no_dbrefs(self):
        mmcif = {1: 2}
        build_struct_ref_seq({
            "1": {"molecules": {"A": {"dbrefs": []}}},
            "2": {"molecules": {"B": {"dbrefs": []}}}
        }, mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_build_table(self):
        mmcif = {
            1: 2, "entry": [{"id": "1XXX"}],
            "struct_ref": [{"id": "1", "pdbx_db_accession": "ABC"}, {"id": "2", "pdbx_db_accession": "XYZ"}]
        }
        build_struct_ref_seq({
            "1": {"id": "1", "molecules": {
                "A": {"dbrefs": [{
                    "id": "1", "start_insert": "", "end_insert": "", "accession": "ABC",
                    "db_start_insert": "", "db_end_insert": "", "start": "1", "end": "2",
                    "db_start": "3", "db_end": "4", "order": 1
                }]},
                "B": {"dbrefs": [{
                    "id": "1", "start_insert": "", "end_insert": "", "accession": "ABC",
                    "db_start_insert": "", "db_end_insert": "", "start": "5", "end": "6",
                    "db_start": "7", "db_end": "8", "order": 2
                }]}
            }},
            "2": {"id": "2", "molecules": {
                "C": {"dbrefs": [{
                    "id": "1", "start_insert": "", "end_insert": "", "accession": "ABC",
                    "db_start_insert": "", "db_end_insert": "", "start": "9", "end": "10",
                    "db_start": "11", "db_end": "12", "order": 3
                }]},
                "D": {"dbrefs": [{
                    "id": "1", "start_insert": "", "end_insert": "", "accession": "XYZ",
                    "db_start_insert": "", "db_end_insert": "", "start": "13", "end": "14",
                    "db_start": "15", "db_end": "16", "order": 4
                }, {
                    "id": "1", "start_insert": "A", "end_insert": "B", "accession": "ABC",
                    "db_start_insert": "C", "db_end_insert": "D", "start": "17", "end": "18",
                    "db_start": "19", "db_end": "20", "order": 5
                }]}
            }},
        }, mmcif)
        self.assertEqual(mmcif, {
            1: 2, "entry": [{"id": "1XXX"}],
            "struct_ref": [{"id": "1", "pdbx_db_accession": "ABC"}, {"id": "2", "pdbx_db_accession": "XYZ"}],
            "struct_ref_seq": [{
                "align_id": "1", "ref_id": "1", "pdbx_PDB_id_code": "1XXX", "pdbx_strand_id": "A",
                "seq_align_beg": "?", "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
                "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "ABC", "db_align_beg": "3",
                "pdbx_db_align_beg_ins_code": "?", "db_align_end": "4", "pdbx_db_align_end_ins_code": "?",
                "pdbx_auth_seq_align_beg": "1", "pdbx_auth_seq_align_end": "2"
            }, {
                "align_id": "2", "ref_id": "1", "pdbx_PDB_id_code": "1XXX", "pdbx_strand_id": "B",
                "seq_align_beg": "?", "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
                "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "ABC", "db_align_beg": "7",
                "pdbx_db_align_beg_ins_code": "?", "db_align_end": "8", "pdbx_db_align_end_ins_code": "?",
                "pdbx_auth_seq_align_beg": "5", "pdbx_auth_seq_align_end": "6"
            }, {
                "align_id": "3", "ref_id": "1", "pdbx_PDB_id_code": "1XXX", "pdbx_strand_id": "C",
                "seq_align_beg": "?", "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
                "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "ABC", "db_align_beg": "11",
                "pdbx_db_align_beg_ins_code": "?", "db_align_end": "12", "pdbx_db_align_end_ins_code": "?",
                "pdbx_auth_seq_align_beg": "9", "pdbx_auth_seq_align_end": "10"
            }, {
                "align_id": "4", "ref_id": "2", "pdbx_PDB_id_code": "1XXX", "pdbx_strand_id": "D", 
                "seq_align_beg": "?", "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
                "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "XYZ", "db_align_beg": "15",
                "pdbx_db_align_beg_ins_code": "?", "db_align_end": "16", "pdbx_db_align_end_ins_code": "?",
                "pdbx_auth_seq_align_beg": "13", "pdbx_auth_seq_align_end": "14"
            }, {
                "align_id": "5", "ref_id": "1", "pdbx_PDB_id_code": "1XXX", "pdbx_strand_id": "D",
                "seq_align_beg": "?", "pdbx_seq_align_beg_ins_code": "A", "seq_align_end": "?",
                "pdbx_seq_align_end_ins_code": "B", "pdbx_db_accession": "ABC", "db_align_beg": "19",
                "pdbx_db_align_beg_ins_code": "C", "db_align_end": "20", "pdbx_db_align_end_ins_code": "D",
                "pdbx_auth_seq_align_beg": "17", "pdbx_auth_seq_align_end": "18"
            }]
        })
    

    def test_can_build_empty_table(self):
        mmcif = {
            1: 2, "struct_ref": [{"id": "1", "pdbx_db_accession": ""}]
        }
        build_struct_ref_seq({
            "1": {"id": "1", "molecules": {
                "A": {"dbrefs": [{
                    "id": "", "start_insert": "", "end_insert": "", "accession": "",
                    "db_start_insert": "", "db_end_insert": "", "start": "", "end": "",
                    "db_start": "", "db_end": "", "order": 0
                }]}
            }}
        }, mmcif)
        self.assertEqual(mmcif, {1: 2,"struct_ref": [{"id": "1", "pdbx_db_accession": ""}], "struct_ref_seq": [{
            "align_id": "1", "ref_id": "1", "pdbx_PDB_id_code": "?", "pdbx_strand_id": "A",
            "seq_align_beg": "?", "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
            "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "?", "db_align_beg": "?",
            "pdbx_db_align_beg_ins_code": "?", "db_align_end": "?", "pdbx_db_align_end_ins_code": "?",
            "pdbx_auth_seq_align_beg": "?", "pdbx_auth_seq_align_end": "?"
        }]})



class StructRefSeqDiffTests(TestCase):

    def test_can_handle_no_entities(self):
        mmcif = {1: 2}
        build_struct_ref_seq_dif({}, mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_handle_no_valid_polymers(self):
        mmcif = {1: 2}
        build_struct_ref_seq_dif({"1": {"molecules": {}}, "2": {"molecules": {}}}, mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_handle_no_differences(self):
        mmcif = {1: 2}
        build_struct_ref_seq_dif({
            "1": {"molecules": {"A": {"differences": []}}},
            "2": {"molecules": {"B": {"differences": []}}}
        }, mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_build_table(self):
        mmcif = {1: 2, "entry": [{"id": "1XXX"}], "struct_ref_seq": [
            {"align_id": "1", "pdbx_db_accession": "ABC", "pdbx_strand_id": "A"},
            {"align_id": "2", "pdbx_db_accession": "XYZ", "pdbx_strand_id": "A"},
            {"align_id": "3", "pdbx_db_accession": "ABC", "pdbx_strand_id": "B"},
            {"align_id": "4", "pdbx_db_accession": "ZZZ", "pdbx_strand_id": "C"},
            {"align_id": "5", "pdbx_db_accession": "YYY", "pdbx_strand_id": "D"},
        ]}
        build_struct_ref_seq_dif({
            "1": {"id": "1", "molecules": {
                "A": {"differences": [{
                    "accession": "XYZ", "name": "HAS", "insert": "A", "database": "UNP",
                    "db_number": "12", "comment": "COM1", "number": "13", "db_name": "HIS"
                }]},
                "B": {"differences": [{
                    "accession": "ABC", "name": "NET", "insert": "", "database": "UNP",
                    "db_number": "22", "comment": "COM2", "number": "23", "db_name": "MET"
                }]}
            }},
            "2": {"id": "2", "molecules": {
                "C": {"differences": [{
                    "accession": "AAA", "name": "SIS", "insert": "", "database": "UNP",
                    "db_number": "32", "comment": "COM3", "number": "33", "db_name": "CYS"
                }]},
                "D": {"differences": [{
                    "accession": "ZZZ", "name": "CRO", "insert": "", "database": "UNP",
                    "db_number": "42", "comment": "COM4", "number": "43", "db_name": "PRO"
                }, {
                    "accession": "YYY", "name": "FAL", "insert": "", "database": "REF",
                    "db_number": "42", "comment": "COM5", "number": "43", "db_name": "VAL"
                }]}
            }},
        }, mmcif)
        self.assertEqual(mmcif, {1: 2, "entry": [{"id": "1XXX"}], "struct_ref_seq": [
            {"align_id": "1", "pdbx_db_accession": "ABC", "pdbx_strand_id": "A"},
            {"align_id": "2", "pdbx_db_accession": "XYZ", "pdbx_strand_id": "A"},
            {"align_id": "3", "pdbx_db_accession": "ABC", "pdbx_strand_id": "B"},
            {"align_id": "4", "pdbx_db_accession": "ZZZ", "pdbx_strand_id": "C"},
            {"align_id": "5", "pdbx_db_accession": "YYY", "pdbx_strand_id": "D"}
        ],
        "struct_ref_seq_dif": [{
            "align_id": "2", "pdbx_pdb_id_code": "1XXX", "mon_id": "HAS", "pdbx_pdb_strand_id": "A",
            "seq_num": "?", "pdbx_pdb_ins_code": "A", "pdbx_seq_db_name": "UNP",
            "pdbx_seq_db_accession_code": "XYZ", "db_mon_id": "HIS", "pdbx_seq_db_seq_num": "12",
            "details": "COM1", "pdbx_auth_seq_num": "13", "pdbx_ordinal": "1"
        }, {
            "align_id": "3", "pdbx_pdb_id_code": "1XXX", "mon_id": "NET", "pdbx_pdb_strand_id": "B",
            "seq_num": "?", "pdbx_pdb_ins_code": "?", "pdbx_seq_db_name": "UNP",
            "pdbx_seq_db_accession_code": "ABC", "db_mon_id": "MET", "pdbx_seq_db_seq_num": "22",
            "details": "COM2", "pdbx_auth_seq_num": "23", "pdbx_ordinal": "2"
        }, {
            "align_id": "?", "pdbx_pdb_id_code": "1XXX", "mon_id": "SIS", "pdbx_pdb_strand_id": "C",
            "seq_num": "?", "pdbx_pdb_ins_code": "?", "pdbx_seq_db_name": "UNP",
            "pdbx_seq_db_accession_code": "AAA", "db_mon_id": "CYS", "pdbx_seq_db_seq_num": "32",
            "details": "COM3", "pdbx_auth_seq_num": "33", "pdbx_ordinal": "3"
        }, {
            "align_id": "?", "pdbx_pdb_id_code": "1XXX", "mon_id": "CRO", "pdbx_pdb_strand_id": "D",
            "seq_num": "?", "pdbx_pdb_ins_code": "?", "pdbx_seq_db_name": "UNP",
            "pdbx_seq_db_accession_code": "ZZZ", "db_mon_id": "PRO", "pdbx_seq_db_seq_num": "42",
            "details": "COM4", "pdbx_auth_seq_num": "43", "pdbx_ordinal": "4"
        }, {
            "align_id": "5", "pdbx_pdb_id_code": "1XXX", "mon_id": "FAL", "pdbx_pdb_strand_id": "D",
            "seq_num": "?", "pdbx_pdb_ins_code": "?", "pdbx_seq_db_name": "REF",
            "pdbx_seq_db_accession_code": "YYY", "db_mon_id": "VAL", "pdbx_seq_db_seq_num": "42",
            "details": "COM5", "pdbx_auth_seq_num": "43", "pdbx_ordinal": "5"
        }]})
    

    def test_can_build_empty_table(self):
        mmcif = {1: 2, "struct_ref_seq": []}
        build_struct_ref_seq_dif({
            "1": {"id": "1", "molecules": {
                "A": {"differences": [{
                    "accession": "", "name": "", "insert": "", "database": "",
                    "db_number": "", "comment": "", "number": "", "db_name": ""
                }]}
            }}
        }, mmcif)
        self.assertEqual(mmcif, {1: 2, "struct_ref_seq": [], "struct_ref_seq_dif": [{
            "align_id": "?", "pdbx_pdb_id_code": "?", "mon_id": "?", "pdbx_pdb_strand_id": "A",
            "seq_num": "?", "pdbx_pdb_ins_code": "?", "pdbx_seq_db_name": "?",
            "pdbx_seq_db_accession_code": "?", "db_mon_id": "?", "pdbx_seq_db_seq_num": "?",
            "details": "?", "pdbx_auth_seq_num": "?", "pdbx_ordinal": "1"
        }]})




class PdbxStructModResidueTests(TestCase):

    def test_can_handle_no_entities(self):
        mmcif = {1: 2}
        build_pdbx_struct_mod_residue({}, mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_handle_no_valid_polymers(self):
        mmcif = {1: 2}
        build_pdbx_struct_mod_residue({"1": {"molecules": {}}, "2": {"molecules": {}}}, mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_handle_no_differences(self):
        mmcif = {1: 2}
        build_pdbx_struct_mod_residue({
            "1": {"molecules": {"A": {"modified": []}}},
            "2": {"molecules": {"B": {"modified": []}}}
        }, mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_build_table(self):
        mmcif = {1: 2}
        build_pdbx_struct_mod_residue({
            "1": {"id": "1", "molecules": {
                "A": {"modified": [
                    {"name": "HYS", "number": "10", "insert": "", "standard_name": "HIS", "comment": "C1"}
                ]},
                "B": {"modified": [
                    {"name": "CRO", "number": "20", "insert": "", "standard_name": "PRO", "comment": "C2"}
                ]}
            }},
            "2": {"id": "2", "molecules": {
                "C": {"modified": [
                    {"name": "MAL", "number": "30", "insert": "A", "standard_name": "VAL", "comment": "C3"}
                ]},
                "D": {"modified": [
                    {"name": "FEE", "number": "40", "insert": "", "standard_name": "PHE", "comment": "C4"},
                    {"name": "SIS", "number": "50", "insert": "", "standard_name": "CYS", "comment": "C5"}
                ]}
            }},
        }, mmcif)
        self.assertEqual(mmcif, {1: 2, "pdbx_struct_mod_residue": [{
            "id": "1", "label_asym_id": "?", "label_comp_id": "HYS", "label_seq_id": "?",
            "auth_asym_id": "A", "auth_comp_id": "HYS", "auth_seq_id": "10", "PDB_ins_code": "?",
            "parent_comp_id": "HIS", "details": "C1"
        }, {
            "id": "2", "label_asym_id": "?", "label_comp_id": "CRO", "label_seq_id": "?", 
            "auth_asym_id": "B", "auth_comp_id": "CRO", "auth_seq_id": "20", "PDB_ins_code": "?",
            "parent_comp_id": "PRO", "details": "C2"
        }, {
            "id": "3", "label_asym_id": "?", "label_comp_id": "MAL", "label_seq_id": "?",
            "auth_asym_id": "C", "auth_comp_id": "MAL", "auth_seq_id": "30", "PDB_ins_code": "A",
            "parent_comp_id": "VAL", "details": "C3"
        }, {
            "id": "4", "label_asym_id": "?", "label_comp_id": "FEE", "label_seq_id": "?",
            "auth_asym_id": "D", "auth_comp_id": "FEE", "auth_seq_id": "40", "PDB_ins_code": "?",
            "parent_comp_id": "PHE", "details": "C4"
        }, {
            "id": "5", "label_asym_id": "?", "label_comp_id": "SIS", "label_seq_id": "?",
            "auth_asym_id": "D", "auth_comp_id": "SIS", "auth_seq_id": "50", "PDB_ins_code": "?",
            "parent_comp_id": "CYS", "details": "C5"
        }]})
    

    def test_can_build_empty_table(self):
        mmcif = {1: 2}
        build_pdbx_struct_mod_residue({
            "1": {"id": "1", "molecules": {
                "A": {"modified": [
                    {"name": "", "number": "", "insert": "", "standard_name": "", "comment": ""}
                ]}
            }}
        }, mmcif)
        self.assertEqual(mmcif, {1: 2, "pdbx_struct_mod_residue": [{
            "id": "1", "label_asym_id": "?", "label_comp_id": "?", "label_seq_id": "?",
            "auth_asym_id": "A", "auth_comp_id": "?", "auth_seq_id": "?", "PDB_ins_code": "?",
            "parent_comp_id": "?", "details": "?"
        }]})



class PdbxEntityNonpolyTests(TestCase):

    def test_can_handle_no_entities(self):
        mmcif = {1: 2}
        build_pdbx_entity_nonpoly({}, mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_handle_no_valid_molecules(self):
        mmcif = {1: 2}
        build_pdbx_entity_nonpoly({"1": {"molecules": {}}, "2": {"molecules": {}}}, mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_build_table(self):
        mmcif = {"entity": [
            {"id": "1", "pdbx_description": "mol1"},
            {"id": "2", "pdbx_description": "mol2"},
        ]}
        entities = {"XMP": {"id": "1", "molecules": {1: 2}}, "BU2": {"id": "2", "molecules": {1: 2}}}
        build_pdbx_entity_nonpoly(entities, mmcif)
        self.assertEqual(mmcif, {
            "entity": [
                {"id": "1", "pdbx_description": "mol1"},
                {"id": "2", "pdbx_description": "mol2"},
            ],
            "pdbx_entity_nonpoly": [
                {"entity_id": "1", "name": "mol1", "comp_id": "XMP"},
                {"entity_id": "2", "name": "mol2", "comp_id": "BU2"},
            ]
        })



class ChemCompTests(TestCase):

    @patch("atomium.pdb.build_ligand_chem_comp")
    @patch("atomium.pdb.build_residue_chem_comp")
    def test_can_handle_no_chem_comp(self, mock_res, mock_lig):
        entities = {1: 2}
        mmcif = {3: 4}
        build_chem_comp(entities, mmcif)
        mock_lig.assert_called_with(entities, mmcif)
        mock_res.assert_called_with(mock_lig.return_value, mmcif)
        self.assertEqual(mmcif, {3: 4})
    

    @patch("atomium.pdb.build_ligand_chem_comp")
    @patch("atomium.pdb.build_residue_chem_comp")
    def test_can_sort_chem_comp(self, mock_res, mock_lig):
        entities = {1: 2}
        mmcif = {3: 4}
        mock_res.side_effect = lambda _, m: [m["chem_comp"].append({"id": n}) for n in ["X", "A", "C"]]
        build_chem_comp(entities, mmcif)
        mock_lig.assert_called_with(entities, mmcif)
        mock_res.assert_called_with(mock_lig.return_value, mmcif)
        self.assertEqual(mmcif, {3: 4, "chem_comp": [{"id": "A"}, {"id": "C"}, {"id": "X"}]})



class LigandChemCompTests(TestCase):

    def test_can_handle_no_entities(self):
        mmcif = {"chem_comp": []}
        entities = {}
        lookup = build_ligand_chem_comp(entities, mmcif)
        self.assertEqual(mmcif, {"chem_comp": []})
        self.assertEqual(lookup, {})
    

    @patch("atomium.pdb.formula_to_weight")
    def test_can_parse_entities(self, mock_weight):
        mock_weight.side_effect = [1.23456, 9.876543]
        mmcif = {"chem_comp": [], "entity": [
            {"id": "3", "pdbx_description": "x.m.p"},
            {"id": "4", "pdbx_description": "b.u.2"},
        ]}
        entities = {
            "XMP": {"id": "3", "molecules": {1: 2}, "formula": "12(C2 N)", "synonyms": ["syn1", "syn2"]},
            "BU2": {"id": "4", "molecules": {1: 2}, "formula": "C5 N"},
            "HOH": {"id": "5", "molecules": {}, "formula": "8(C9 CL)"},
        }
        lookup = build_ligand_chem_comp(entities, mmcif)
        self.assertEqual(lookup, {"HOH": "C9 CL"})
        self.assertEqual(mmcif, {
            "chem_comp": [{
                "id": "XMP", "type": "non-polymer", "mon_nstd_flag": ".",
                "name": "X.M.P",
                "pdbx_synonyms": "syn1,syn2",
                "formula": "C2 N",
                "formula_weight": f"1.235"
            }, {
                "id": "BU2", "type": "non-polymer", "mon_nstd_flag": ".",
                "name": "B.U.2",
                "pdbx_synonyms": "?",
                "formula": "C5 N",
                "formula_weight": f"9.877"
            }],
            "entity": [
                {"id": "3", "pdbx_description": "x.m.p"},
                {"id": "4", "pdbx_description": "b.u.2"},
            ]
        })



class ResidueChemCompTests(TestCase):

    @patch("atomium.pdb.formula_to_weight")
    def test_can_build_residue_chem_comp(self, mock_weight):
        mmcif = {"chem_comp": [], "entity_poly_seq": [
            {"mon_id": "HIS"}, {"mon_id": "HIS"}, {"mon_id": "VAL"},
            {"mon_id": "MOD"}, {"mon_id": "XYZ"}
        ]}
        mock_weight.return_value = 0.123456
        formulae_lookup = {"MOD": "C2 P", "ABC": "AA"}
        build_residue_chem_comp(formulae_lookup, mmcif)
        self.assertEqual(mmcif, {
            "chem_comp": [{
                "id": "HIS", "type": "L-peptide linking", "mon_nstd_flag": "y", "name": "HISTIDINE",
                "pdbx_synonyms": "?", "formula": "C6 H10 N3 O2 1", "formula_weight": "156.162"
            }, {
                "id": "MOD", "type": "L-peptide linking", "mon_nstd_flag": "n", "name": "?",
                "pdbx_synonyms": "?", "formula": "C2 P", "formula_weight": "0.123"
            }, {
                "id": "VAL", "type": "L-peptide linking", "mon_nstd_flag": "y", "name": "VALINE",
                "pdbx_synonyms": "?", "formula": "C5 H11 N O2", "formula_weight": "117.146"
            }, {
                "id": "XYZ", "type": "L-peptide linking", "mon_nstd_flag": "n", "name": "?",
                "pdbx_synonyms": "?", "formula": "?", "formula_weight": "?"
            }],
            "entity_poly_seq": [
                {"mon_id": "HIS"}, {"mon_id": "HIS"}, {"mon_id": "VAL"},
                {"mon_id": "MOD"}, {"mon_id": "XYZ"}
            ]
        })



class FormulaToWeightTests(TestCase):

    def test_can_weight_of_atom(self):
        self.assertEqual(formula_to_weight("C"), 12.0107)
        self.assertEqual(formula_to_weight("N"), 14.0067)
        self.assertEqual(formula_to_weight("X"), 0)
    

    def test_can_get_weight_of_molecule(self):
        self.assertAlmostEqual(formula_to_weight("C6 H15 N4 O2 1"), 175.2083, delta=0.00001)
        self.assertAlmostEqual(formula_to_weight("ZN2 2+"), 130.78, delta=0.00001)
        self.assertAlmostEqual(formula_to_weight("ZN2 2-"), 130.78, delta=0.00001)



class AtomTypeTests(TestCase):

    def test_can_handle_no_atoms(self):
        mmcif = {1: 2}
        build_atom_type("", mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_parse_atoms(self):
        filestring = (
            "ATOM    269  CA  GLY A  44      -8.089  40.385  49.069  1.00 13.57           C  \n"
            "ATOM    270  N   TYR A  44      -9.857  41.667  48.097  1.00 14.09          PN  \n"
            "ATOM    270  C   GLY A  44      -9.512  40.482  48.572  1.00 14.60           C  \n"
            "HETATM  272  O   GLY A  45     -10.281  39.509  48.629  1.00 13.05           O  \n"
        )
        mmcif = {1: 2}
        build_atom_type(filestring, mmcif)
        self.assertEqual(mmcif, {1: 2, "atom_type": [
            {"symbol": "C"}, {"symbol": "O"}, {"symbol": "PN"}
        ]})



class AtomSiteTests(TestCase):

    def test_can_handle_no_atoms(self):
        mmcif = {1: 2}
        build_atom_site("", {}, {}, mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    @patch("atomium.pdb.get_atom_entity")
    @patch("atomium.pdb.get_atom_label")
    @patch("atomium.pdb.add_atom")
    def test_can_build_single_model(self, mock_add, mock_label, mock_entity):
        filestring = (
            "ATOM      1  N   VAL A  11       3.696  33.898  63.219  1.00 21.50           N  \n"
            "ATOM      2  CA  VAL A  11       3.198  33.218  61.983  1.00 19.76           C  \n"
            "ATOM      8  N   MET A  12       3.155  30.797  61.557  1.00 17.03           N  \n"
            "HETATM 3194  C1  BU2 A5001       2.646  45.112  48.995  1.00 43.24           C  \n"
            "HETATM 3195  O1  BU2 A5001       1.781  45.484  47.929  1.00 42.82           O \n"
        )
        polymers = {"1": {}}
        non_polymers = {"XMP": {}}
        mmcif = {1: 2}
        mock_entity.side_effect = [
            [{"id": "1"}, {"alignment": [1, 2, 3]}, None],
            [{"id": "1"}, {"alignment": [4, 5, 6]}, None],
            [{"id": "2"}, {"alignment": [7, 8, 9]}, None],
            [None, None, {"id": "3"}],
            [None, None, {"id": "4"}],
        ]
        mock_label.side_effect = ["A", "A", "B", "C", "D"]
        build_atom_site(filestring, polymers, non_polymers, mmcif)
        self.assertEqual([c[0] for c in mock_add.call_args_list], [
            (mmcif, filestring.splitlines()[0], "A", "1", "2", "1"),
            (mmcif, filestring.splitlines()[1], "A", "1", "5", "1"),
            (mmcif, filestring.splitlines()[2], "B", "2", "9", "1"),
            (mmcif, filestring.splitlines()[3], "C", "3", ".", "1"),
            (mmcif, filestring.splitlines()[4], "D", "4", ".", "1"),
        ])
    

    @patch("atomium.pdb.get_atom_entity")
    @patch("atomium.pdb.get_atom_label")
    @patch("atomium.pdb.add_atom")
    def test_can_build_multi_model(self, mock_add, mock_label, mock_entity):
        filestring = (
            "MODEL     1\n"
            "ATOM      1  N   VAL A  11       3.696  33.898  63.219  1.00 21.50           N  \n"
            "ATOM      2  CA  VAL A  11       3.198  33.218  61.983  1.00 19.76           C  \n"
            "ATOM      8  N   MET A  12       3.155  30.797  61.557  1.00 17.03           N  \n"
            "HETATM 3194  C1  BU2 A5001       2.646  45.112  48.995  1.00 43.24           C  \n"
            "HETATM 3195  O1  BU2 A5001       1.781  45.484  47.929  1.00 42.82           O \n"
            "ENDMDL  \n"
            "MODEL     2\n"
            "ATOM      1  N   VAL A  11       3.696  33.898  63.219  1.00 21.50           N  \n"
            "ATOM      8  N   MET A  12       3.155  30.797  61.557  1.00 17.03           N  \n"
            "HETATM 3194  C1  BU2 A5001       2.646  45.112  48.995  1.00 43.24           C  \n"
            "ENDMDL  \n"
        )
        polymers = {"1": {}}
        non_polymers = {"XMP": {}}
        mmcif = {1: 2}
        mock_entity.side_effect = [
            [{"id": "1"}, {"alignment": [1, 2, 3]}, None],
            [{"id": "1"}, {"alignment": [4, 5, 6]}, None],
            [{"id": "2"}, {"alignment": [7, 8, 9]}, None],
            [None, None, {"id": "3"}],
            [None, None, {"id": "4"}],
            [{"id": "1"}, {"alignment": [1, 2, 3]}, None],
            [{"id": "2"}, {"alignment": [7, 8, 9]}, None],
            [None, None, {"id": "3"}],
        ]
        mock_label.side_effect = ["A", "A", "B", "C", "D", "A", "B", "C"]
        build_atom_site(filestring, polymers, non_polymers, mmcif)
        self.assertEqual([c[0] for c in mock_add.call_args_list], [
            (mmcif, filestring.splitlines()[1], "A", "1", "2", "1"),
            (mmcif, filestring.splitlines()[2], "A", "1", "5", "1"),
            (mmcif, filestring.splitlines()[3], "B", "2", "9", "1"),
            (mmcif, filestring.splitlines()[4], "C", "3", ".", "1"),
            (mmcif, filestring.splitlines()[5], "D", "4", ".", "1"),
            (mmcif, filestring.splitlines()[8], "A", "1", "2", "2"),
            (mmcif, filestring.splitlines()[9], "B", "2", "9", "2"),
            (mmcif, filestring.splitlines()[10], "C", "3", ".", "2"),
        ])



class AtomEntityTests(TestCase):

    def test_can_get_non_polymer(self):
        sig = ("D", "HIS", "101", "A")
        polymers = {}
        non_polymers = {
            "XMP": {"molecules": {1: 2, 3: 4}},
            "BU2": {"molecules": {1: 2, ("D", "HIS", "101", "A"): 4}},
        }
        polymers = {
            "1": {"molecules": {1: 2, 3: 4}},
            "2": {"molecules": {1: 2, "D": 4}}
        }
        entity = get_atom_entity(sig, polymers, non_polymers)
        self.assertEqual(entity, (None, None, {"molecules": {1: 2, ("D", "HIS", "101", "A"): 4}}))
    

    def test_can_get_polymer(self):
        sig = ("D", "HIS", "101", "A")
        polymers = {}
        non_polymers = {
            "XMP": {"molecules": {1: 2, 3: 4}},
            "BU2": {"molecules": {1: 2, 3: 4}},
        }
        polymers = {
            "1": {"molecules": {1: 2, 3: 4}},
            "2": {"molecules": {1: 2, "D": 4}}
        }
        entity = get_atom_entity(sig, polymers, non_polymers)
        self.assertEqual(entity, ({"molecules": {1: 2, "D": 4}}, 4, None))



class AtomLabelTests(TestCase):

    def setUp(self):
        self.patch1 = patch("atomium.pdb.next_id")
        self.mock_next = self.patch1.start()
        self.mock_next.return_value = "NEXT"
    

    def tearDown(self):
        self.patch1.stop()


    def test_can_get_first_polymer_label(self):
        sig = ("D", "HIS", "101", "A")
        polymer = {1: 2}
        non_polymer = {}
        labels = {}
        label = get_atom_label(sig, polymer, non_polymer, labels)
        self.assertEqual(label, "NEXT")
        self.assertEqual(labels, {"D": "NEXT"})
        self.mock_next.assert_called_with("@")
    

    def test_can_get_existing_polymer_label(self):
        sig = ("D", "HIS", "101", "A")
        polymer = {1: 2}
        non_polymer = {}
        labels = {"D": "X"}
        label = get_atom_label(sig, polymer, non_polymer, labels)
        self.assertEqual(label, "X")
        self.assertEqual(labels, {"D": "X"})
        self.mock_next.assert_called_with("X")
    

    def test_can_get_new_polymer_label(self):
        sig = ("D", "HIS", "101", "A")
        polymer = {1: 2}
        non_polymer = {}
        labels = {"A": "G", "B": "AB", "C": "BAA"}
        label = get_atom_label(sig, polymer, non_polymer, labels)
        self.assertEqual(label, "NEXT")
        self.assertEqual(labels, {"A": "G", "B": "AB", "C": "BAA", "D": "NEXT"})
        self.mock_next.assert_called_with("BAA")
    

    def test_can_get_first_non_polymer_label(self):
        sig = ("D", "XMP", "101", "A")
        polymer = None
        non_polymer = {"is_water": False}
        labels = {}
        label = get_atom_label(sig, polymer, non_polymer, labels)
        self.assertEqual(label, "NEXT")
        self.assertEqual(labels, {("D", "XMP", "101", "A"): "NEXT"})
        self.mock_next.assert_called_with("@")
    

    def test_can_get_existing_non_polymer_label(self):
        sig = ("D", "HIS", "101", "A")
        polymer = None
        non_polymer = {"is_water": False}
        labels = {("D", "HIS", "101", "A"): "X"}
        label = get_atom_label(sig, polymer, non_polymer, labels)
        self.assertEqual(label, "X")
        self.assertEqual(labels, {("D", "HIS", "101", "A"): "X"})
        self.mock_next.assert_called_with("X")
    

    def test_can_get_new_non_polymer_label(self):
        sig = ("D", "HIS", "102", "A")
        polymer = None
        non_polymer = {"is_water": False}
        labels = {"A": "G", "B": "AB", "C": "BAA"}
        label = get_atom_label(sig, polymer, non_polymer, labels)
        self.assertEqual(label, "NEXT")
        self.assertEqual(labels, {"A": "G", "B": "AB", "C": "BAA", ("D", "HIS", "102", "A"): "NEXT"})
        self.mock_next.assert_called_with("BAA")
    

    def test_can_get_first_water_label(self):
        sig = ("D", "XMP", "101", "A")
        polymer = None
        non_polymer = {"is_water": True}
        labels = {}
        label = get_atom_label(sig, polymer, non_polymer, labels)
        self.assertEqual(label, "NEXT")
        self.assertEqual(labels, {("D", "W"): "NEXT"})
        self.mock_next.assert_called_with("@")
    

    def test_can_get_existing_water_label(self):
        sig = ("D", "HIS", "101", "A")
        polymer = None
        non_polymer = {"is_water": True}
        labels = {"D": "X", ("D", "W"): "Y"}
        label = get_atom_label(sig, polymer, non_polymer, labels)
        self.assertEqual(label, "Y")
        self.assertEqual(labels, {"D": "X", ("D", "W"): "Y"})
        self.mock_next.assert_called_with("Y")
    

    def test_can_get_new_water_label(self):
        sig = ("D", "HIS", "102", "A")
        polymer = None
        non_polymer = {"is_water": True}
        labels = {"D": "X"}
        label = get_atom_label(sig, polymer, non_polymer, labels)
        self.assertEqual(label, "NEXT")
        self.assertEqual(labels, {"D": "X", ("D", "W"): "NEXT"})
        self.mock_next.assert_called_with("X")



class NextIdTests(TestCase):

    def test_can_get_first_id(self):
        self.assertEqual(next_id("@"), "A")
    

    def test_can_get_single_letter_id(self):
        self.assertEqual(next_id("A"), "B")
        self.assertEqual(next_id("B"), "C")
        self.assertEqual(next_id("M"), "N")
        self.assertEqual(next_id("Y"), "Z")
    

    def test_can_get_two_letter_id(self):
        self.assertEqual(next_id("Z"), "AA")
        self.assertEqual(next_id("AA"), "BA")
        self.assertEqual(next_id("BA"), "CA")
        self.assertEqual(next_id("ZA"), "AB")
        self.assertEqual(next_id("AB"), "BB")
        self.assertEqual(next_id("BB"), "CB")
    

    def test_can_get_three_letter_id(self):
        self.assertEqual(next_id("ZZ"), "AAA")
        self.assertEqual(next_id("AAA"), "BAA")
        self.assertEqual(next_id("BAA"), "CAA")
        self.assertEqual(next_id("ZAA"), "ABA")
        self.assertEqual(next_id("ABA"), "BBA")
        self.assertEqual(next_id("ZBA"), "ACA")
        self.assertEqual(next_id("ZZA"), "AAB")



class AtomAddingTests(TestCase):

    def test_can_add_full_atom(self):
        mmcif = {"atom_site": []}
        line = "ATOM      1  N  ZVAL A  11F      3.696  33.898  63.219  1.00 21.50           N+2"
        add_atom(mmcif, line, "P", "4", "123", "4")   
        self.assertEqual(mmcif, {"atom_site": [{
            "group_PDB": "ATOM", "id": "1", "type_symbol": "N", "label_atom_id": "N",
            "label_alt_id": "Z", "label_comp_id": "VAL", "label_asym_id": "P", "label_entity_id": "4",
            "label_seq_id": "123", "pdbx_PDB_ins_code": "F", "Cartn_x": "3.696", "Cartn_y": "33.898",
            "Cartn_z": "63.219", "occupancy": "1.00", "B_iso_or_equiv": "21.50",
            "pdbx_formal_charge": "+2", "auth_seq_id": "11", "auth_comp_id": "VAL",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "4"
        }]})
    

    def test_can_add_minimal_hetatm(self):
        mmcif = {"atom_site": []}
        line = "HETATM    1  N   VAL A  11       3.696  33.898  63.219  1.00 21.50           N"
        add_atom(mmcif, line, "P", "4", "123", "4")   
        self.assertEqual(mmcif, {"atom_site": [{
            "group_PDB": "HETATM", "id": "1", "type_symbol": "N", "label_atom_id": "N",
            "label_alt_id": ".", "label_comp_id": "VAL", "label_asym_id": "P", "label_entity_id": "4",
            "label_seq_id": "123", "pdbx_PDB_ins_code": "?", "Cartn_x": "3.696", "Cartn_y": "33.898",
            "Cartn_z": "63.219", "occupancy": "1.00", "B_iso_or_equiv": "21.50",
            "pdbx_formal_charge": "?", "auth_seq_id": "11", "auth_comp_id": "VAL",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "4"
        }]})



class AtomSiteAnisotropTests(TestCase):

    def test_can_handle_no_anisou(self):
        mmcif = {1: 2}
        build_atom_site_anisotrop("", mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test(self):
        filestring = (
            "ANISOU   14  N   TYR A  74      657    661    657    178    -47     84       N  \n"
            "ANISOU   15  CA  TYR A  74      477    497    479     92    -20    105       C  \n"
            "ANISOU   16  C   TYR A  74      483    432    493     77     25     95       C  \n"
        )
        mmcif = {"atom_site": [{
            "id": "14", "type_symbol": "N", "label_atom_id": "NN", "label_alt_id": ".",
            "label_comp_id": "XYZ", "label_asym_id": "B", "label_seq_id": "10", "pdbx_PDB_ins_code": "?",
            "auth_seq_id": "101", "auth_comp_id": "AAA", "auth_asym_id": "A", "auth_atom_id": "X"
        }, {
            "id": "15", "type_symbol": "C", "label_atom_id": "CC", "label_alt_id": "R",
            "label_comp_id": "HIS", "label_asym_id": "C", "label_seq_id": "20", "pdbx_PDB_ins_code": "S",
            "auth_seq_id": "200", "auth_comp_id": "PRO", "auth_asym_id": "D", "auth_atom_id": "L"
        }]}
        build_atom_site_anisotrop(filestring, mmcif)
        self.assertEqual(mmcif, {
            "atom_site": [{
                "id": "14", "type_symbol": "N", "label_atom_id": "NN", "label_alt_id": ".",
                "label_comp_id": "XYZ", "label_asym_id": "B", "label_seq_id": "10", "pdbx_PDB_ins_code": "?",
                "auth_seq_id": "101", "auth_comp_id": "AAA", "auth_asym_id": "A", "auth_atom_id": "X"
            }, {
                "id": "15", "type_symbol": "C", "label_atom_id": "CC", "label_alt_id": "R",
                "label_comp_id": "HIS", "label_asym_id": "C", "label_seq_id": "20", "pdbx_PDB_ins_code": "S",
                "auth_seq_id": "200", "auth_comp_id": "PRO", "auth_asym_id": "D", "auth_atom_id": "L"
            }],
            "atom_site_anisotrop": [{
                "id": "14", "type_symbol": "N", "pdbx_label_atom_id": "NN", "pdbx_label_alt_id": ".",
                "pdbx_label_comp_id": "XYZ", "pdbx_label_asym_id": "B", "pdbx_label_seq_id": "10",
                "pdbx_PDB_ins_code": "?", "U[1][1]": "0.0657", "U[2][2]": "0.0661", "U[3][3]": "0.0657",
                "U[1][2]": "0.0178", "U[1][3]": "-0.0047", "U[2][3]": "0.0084", "pdbx_auth_seq_id": "101",
                "pdbx_auth_comp_id": "AAA", "pdbx_auth_asym_id": "A", "pdbx_auth_atom_id ": "X"
            }, {
                "id": "15", "type_symbol": "C", "pdbx_label_atom_id": "CC", "pdbx_label_alt_id": "R",
                "pdbx_label_comp_id": "HIS", "pdbx_label_asym_id": "C", "pdbx_label_seq_id": "20",
                "pdbx_PDB_ins_code": "S", "U[1][1]": "0.0477", "U[2][2]": "0.0497", "U[3][3]": "0.0479",
                "U[1][2]": "0.0092", "U[1][3]": "-0.002", "U[2][3]": "0.0105", "pdbx_auth_seq_id": "200",
                "pdbx_auth_comp_id": "PRO", "pdbx_auth_asym_id": "D", "pdbx_auth_atom_id ": "L"
            }]
        })



class AtomIdUpdatingTests(TestCase):

    def test_can_update_atom_ids(self):
        mmcif = {
            "atom_site": [{"id": "2"}, {"id": "A"}, {"id": "9"}],
            "atom_site_anisotrop": [{"id": "A"}, {"id": "2"}]
        }
        update_atom_ids(mmcif)
        self.assertEqual(mmcif, {
            "atom_site": [{"id": "1"}, {"id": "2"}, {"id": "3"}],
            "atom_site_anisotrop": [{"id": "2"}, {"id": "1"}]
        })
    

    def test_can_handle_no_anisotrop(self):
        mmcif = {"atom_site": [{"id": "2"}, {"id": "A"}, {"id": "9"}]}
        update_atom_ids(mmcif)
        self.assertEqual(mmcif, {"atom_site": [{"id": "1"}, {"id": "2"}, {"id": "3"}]})
    

    def test_can_handle_no_atoms(self):
        mmcif = {1: 2}
        update_atom_ids(mmcif)
        self.assertEqual(mmcif, {1: 2})



class StructAsymTests(TestCase):

    def test_can_handle_no_atoms(self):
        mmcif = {1: 2}
        build_struct_asym(mmcif)
        self.assertEqual(mmcif, {1: 2})


    def test_can_make_struct_asym(self):
        mmcif = {"atom_site": [
            {"label_asym_id": "A", "label_entity_id": "1"},
            {"label_asym_id": "B", "label_entity_id": "2"},
            {"label_asym_id": "C", "label_entity_id": "1"},
            {"label_asym_id": "D", "label_entity_id": "2"},
            {"label_asym_id": "E", "label_entity_id": "3"},
            {"label_asym_id": "F", "label_entity_id": "4"},
            {"label_asym_id": "G", "label_entity_id": "5"},
            {"label_asym_id": "H", "label_entity_id": "5"},
        ]}
        build_struct_asym(mmcif)
        self.assertEqual(mmcif, {"atom_site": [
            {"label_asym_id": "A", "label_entity_id": "1"},
            {"label_asym_id": "B", "label_entity_id": "2"},
            {"label_asym_id": "C", "label_entity_id": "1"},
            {"label_asym_id": "D", "label_entity_id": "2"},
            {"label_asym_id": "E", "label_entity_id": "3"},
            {"label_asym_id": "F", "label_entity_id": "4"},
            {"label_asym_id": "G", "label_entity_id": "5"},
            {"label_asym_id": "H", "label_entity_id": "5"},
        ], "struct_asym": [{
            "id": "A", "pdbx_blank_PDB_chainid_flag": "N",
            "pdbx_modified": "N", "entity_id": "1", "details": "?"
        }, {
            "id": "B", "pdbx_blank_PDB_chainid_flag": "N",
            "pdbx_modified": "N", "entity_id": "2", "details": "?"
        }, {
            "id": "C", "pdbx_blank_PDB_chainid_flag": "N",
            "pdbx_modified": "N", "entity_id": "1", "details": "?"
        }, {
            "id": "D", "pdbx_blank_PDB_chainid_flag": "N",
            "pdbx_modified": "N", "entity_id": "2", "details": "?"
        }, {
            "id": "E", "pdbx_blank_PDB_chainid_flag": "N",
            "pdbx_modified": "N", "entity_id": "3", "details": "?"
        }, {
            "id": "F", "pdbx_blank_PDB_chainid_flag": "N",
            "pdbx_modified": "N", "entity_id": "4", "details": "?"
        }, {
            "id": "G", "pdbx_blank_PDB_chainid_flag": "N",
            "pdbx_modified": "N", "entity_id": "5", "details": "?"
        }, {
            "id": "H", "pdbx_blank_PDB_chainid_flag": "N",
            "pdbx_modified": "N", "entity_id": "5", "details": "?"
        }]})



class AnnotationParsingTests(TestCase):

    @patch("atomium.pdb.parse_helix")
    @patch("atomium.pdb.parse_sheet")
    @patch("atomium.pdb.parse_ssbond")
    @patch("atomium.pdb.parse_link")
    @patch("atomium.pdb.parse_cispep")
    @patch("atomium.pdb.parse_site")
    def test_can_parse_annotation(self, *mocks):
        parse_annotation("filestring", {"mmcif": 1})
        for mock in mocks:
            mock.assert_called_with("filestring", {"mmcif": 1})



class HelixParsingTests(TestCase):

    def test_can_handle_no_helix(self):
        mmcif = {1: 2}
        parse_helix("", mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_parse_helix(self):
        filestring = (
            "HELIX    6 AE1 LYS A  130S GLY A  137Q 1CONTIGUOUS WITH HELIX AE2          8    \n"
            "HELIX    7 AE2 ILE A  138T VAL A  154R 1CONTIGUOUS WITH HELIX AE1         17    \n"
        )
        mmcif = {1: 2}
        parse_helix(filestring, mmcif)
        self.assertEqual(mmcif, {1: 2, "struct_conf": [{
            "conf_type_id": "HELX_P", "id": "HELX_P1", "pdbx_PDB_helix_id": "1", "beg_label_comp_id": "LYS",
            "beg_label_asym_id": "A", "beg_label_seq_id": "130", "pdbx_beg_PDB_ins_code": "S",
            "end_label_comp_id": "GLY", "end_label_asym_id": "A", "end_label_seq_id": "137",
            "pdbx_end_PDB_ins_code": "Q", "beg_auth_comp_id": "LYS", "beg_auth_asym_id": "A",
            "beg_auth_seq_id": "130", "end_auth_comp_id": "GLY", "end_auth_asym_id": "A",
            "end_auth_seq_id": "137", "pdbx_PDB_helix_class": "1", "details": "CONTIGUOUS WITH HELIX AE2", "pdbx_PDB_helix_length": "8"
        }, {
            "conf_type_id": "HELX_P", "id": "HELX_P2", "pdbx_PDB_helix_id": "2", "beg_label_comp_id": "ILE",
            "beg_label_asym_id": "A", "beg_label_seq_id": "138", "pdbx_beg_PDB_ins_code": "T",
            "end_label_comp_id": "VAL", "end_label_asym_id": "A", "end_label_seq_id": "154",
            "pdbx_end_PDB_ins_code": "R", "beg_auth_comp_id": "ILE", "beg_auth_asym_id": "A",
            "beg_auth_seq_id": "138", "end_auth_comp_id": "VAL", "end_auth_asym_id": "A",
            "end_auth_seq_id": "154", "pdbx_PDB_helix_class": "1", "details": "CONTIGUOUS WITH HELIX AE1", "pdbx_PDB_helix_length": "17"
        }]})
    

    def test_can_parse_empty_helix(self):
        filestring = (
            "HELIX    6 AE1 LYS A  130  GLY A  137                                      8    \n"
            "HELIX    7 AE2 ILE A  138  VAL A  154                                     17    \n"
        )
        mmcif = {1: 2}
        parse_helix(filestring, mmcif)
        self.assertEqual(mmcif, {1: 2, "struct_conf": [{
            "conf_type_id": "HELX_P", "id": "HELX_P1", "pdbx_PDB_helix_id": "1", "beg_label_comp_id": "LYS",
            "beg_label_asym_id": "A", "beg_label_seq_id": "130", "pdbx_beg_PDB_ins_code": "?",
            "end_label_comp_id": "GLY", "end_label_asym_id": "A", "end_label_seq_id": "137",
            "pdbx_end_PDB_ins_code": "?", "beg_auth_comp_id": "LYS", "beg_auth_asym_id": "A",
            "beg_auth_seq_id": "130", "end_auth_comp_id": "GLY", "end_auth_asym_id": "A",
            "end_auth_seq_id": "137", "pdbx_PDB_helix_class": "?", "details": "?", "pdbx_PDB_helix_length": "8"
        }, {
            "conf_type_id": "HELX_P", "id": "HELX_P2", "pdbx_PDB_helix_id": "2", "beg_label_comp_id": "ILE",
            "beg_label_asym_id": "A", "beg_label_seq_id": "138", "pdbx_beg_PDB_ins_code": "?",
            "end_label_comp_id": "VAL", "end_label_asym_id": "A", "end_label_seq_id": "154",
            "pdbx_end_PDB_ins_code": "?", "beg_auth_comp_id": "ILE", "beg_auth_asym_id": "A",
            "beg_auth_seq_id": "138", "end_auth_comp_id": "VAL", "end_auth_asym_id": "A",
            "end_auth_seq_id": "154", "pdbx_PDB_helix_class": "?", "details": "?", "pdbx_PDB_helix_length": "17"
        }]})



class SheetParsingTests(TestCase):

    def test_can_handle_no_sheet(self):
        mmcif = {1: 2}
        parse_sheet("", mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    @patch("atomium.pdb.parse_inter_strand")
    def test_can_parse_sheet(self, mock_strand):
        filestring = (
            "SHEET    1   A 9 LEU A  15A MET A  19B 0                              \n"          
            "SHEET    2   A 9 THR A  40C GLY A  44D 1  O  LYS A  42E  N  LEU A  17F\n"                   
            "SHEET    1   B 9 LEU B1015G ALA B1018H 0                              \n"          
            "SHEET    2   B 9 THR B1040  GLY B1044  1  O  LYS B1042   N  LEU B1017 \n"
        )
        mmcif = {1: 2}
        parse_sheet(filestring, mmcif)
        self.assertEqual(mmcif, {1: 2, "struct_sheet": [
            {"id": "A", "type": "?", "number_strands": "9", "details": "?"},
            {"id": "B", "type": "?", "number_strands": "9", "details": "?"}
        ], "struct_sheet_range": [{
            "sheet_id": "A", "id": "1", "beg_label_comp_id": "LEU", "beg_label_asym_id": "A",
            "beg_label_seq_id": "15", "pdbx_beg_PDB_ins_code": "A", "end_label_comp_id": "MET",
            "end_label_asym_id": "A", "end_label_seq_id": "19", "pdbx_end_PDB_ins_code": "B",
            "beg_auth_comp_id": "LEU", "beg_auth_asym_id": "A", "beg_auth_seq_id": "15",
            "end_auth_comp_id": "MET", "end_auth_asym_id": "A", "end_auth_seq_id": "19"
        }, {
            "sheet_id": "A", "id": "2", "beg_label_comp_id": "THR", "beg_label_asym_id": "A",
            "beg_label_seq_id": "40", "pdbx_beg_PDB_ins_code": "C", "end_label_comp_id": "GLY",
            "end_label_asym_id": "A", "end_label_seq_id": "44", "pdbx_end_PDB_ins_code": "D",
            "beg_auth_comp_id": "THR", "beg_auth_asym_id": "A", "beg_auth_seq_id": "40",
            "end_auth_comp_id": "GLY", "end_auth_asym_id": "A", "end_auth_seq_id": "44"
        }, {
            "sheet_id": "B", "id": "1", "beg_label_comp_id": "LEU", "beg_label_asym_id": "B",
            "beg_label_seq_id": "1015", "pdbx_beg_PDB_ins_code": "G", "end_label_comp_id": "ALA",
            "end_label_asym_id": "B", "end_label_seq_id": "1018", "pdbx_end_PDB_ins_code": "H",
            "beg_auth_comp_id": "LEU", "beg_auth_asym_id": "B", "beg_auth_seq_id": "1015",
            "end_auth_comp_id": "ALA", "end_auth_asym_id": "B", "end_auth_seq_id": "1018"
        }, {
            "sheet_id": "B", "id": "2", "beg_label_comp_id": "THR", "beg_label_asym_id": "B",
            "beg_label_seq_id": "1040", "pdbx_beg_PDB_ins_code": "?", "end_label_comp_id": "GLY",
            "end_label_asym_id": "B", "end_label_seq_id": "1044", "pdbx_end_PDB_ins_code": "?",
            "beg_auth_comp_id": "THR", "beg_auth_asym_id": "B", "beg_auth_seq_id": "1040",
            "end_auth_comp_id": "GLY", "end_auth_asym_id": "B", "end_auth_seq_id": "1044"
        }]})
        mock_strand.assert_called_with(filestring.splitlines(), mmcif)



class StrandParsingTests(TestCase):

    def test_can_parse_strands(self):
        lines = [
            "SHEET    1   A 9 LEU A  15A MET A  19B 0                              ",
            "SHEET    2   A 9 THR A  40C GLY A  44D 1  O  LYS A  42E  N  LEU A  17F",
            "SHEET    3   A 9 ARG A  66  VAL A  73  1  O  ILE A  68   N  VAL A  41 ",
            "SHEET    1   B 9 LEU B1015G ALA B1018H 0                              ",
            "SHEET    2   B 9 THR B1040  GLY B1044  -  O  LYS B1042   N  LEU B1017 ",
            "SHEET    3   B 9 ARG B1066  VAL B1073  -  O  ARG B1066L  N  VAL B1041Q"
        ]
        mmcif = {1: 2}
        parse_inter_strand(lines, mmcif)
        self.assertEqual(mmcif, {1: 2, "struct_sheet_order": [
            {"sheet_id": "A", "range_id_1": "1", "range_id_2": "2", "offset": "?", "sense": "parallel"},
            {"sheet_id": "A", "range_id_1": "2", "range_id_2": "3", "offset": "?", "sense": "parallel"},
            {"sheet_id": "B", "range_id_1": "1", "range_id_2": "2", "offset": "?", "sense": "anti-parallel"},
            {"sheet_id": "B", "range_id_1": "2", "range_id_2": "3", "offset": "?", "sense": "anti-parallel"}
        ],
        "pdbx_struct_sheet_hbond": [{
            "sheet_id": "A", "range_id_1": "1", "range_id_2": "2", "range_1_label_atom_id": "N", "range_1_label_comp_id": "LEU",
            "range_1_label_asym_id": "A", "range_1_label_seq_id": "17", "range_1_PDB_ins_code": "F", "range_1_auth_atom_id": "N",
            "range_1_auth_comp_id": "LEU", "range_1_auth_asym_id": "A", "range_1_auth_seq_id": "17", "range_2_label_atom_id": "O",
            "range_2_label_comp_id": "LYS", "range_2_label_asym_id": "A", "range_2_label_seq_id": "44", "range_2_PDB_ins_code": "E",
            "range_2_auth_atom_id": "O", "range_2_auth_comp_id": "LYS", "range_2_auth_asym_id": "A", "range_2_auth_seq_id": "44"
        }, {
            "sheet_id": "A", "range_id_1": "2", "range_id_2": "3", "range_1_label_atom_id": "N", "range_1_label_comp_id": "VAL",
            "range_1_label_asym_id": "A", "range_1_label_seq_id": "41", "range_1_PDB_ins_code": "?", "range_1_auth_atom_id": "N",
            "range_1_auth_comp_id": "VAL", "range_1_auth_asym_id": "A", "range_1_auth_seq_id": "41", "range_2_label_atom_id": "O",
            "range_2_label_comp_id": "ILE", "range_2_label_asym_id": "A", "range_2_label_seq_id": "73", "range_2_PDB_ins_code": "?",
            "range_2_auth_atom_id": "O", "range_2_auth_comp_id": "ILE", "range_2_auth_asym_id": "A", "range_2_auth_seq_id": "73"
        }, {
            "sheet_id": "B", "range_id_1": "1", "range_id_2": "2", "range_1_label_atom_id": "N", "range_1_label_comp_id": "LEU",
            "range_1_label_asym_id": "B", "range_1_label_seq_id": "1017", "range_1_PDB_ins_code": "?", "range_1_auth_atom_id": "N",
            "range_1_auth_comp_id": "LEU", "range_1_auth_asym_id": "B", "range_1_auth_seq_id": "1017", "range_2_label_atom_id": "O",
            "range_2_label_comp_id": "LYS", "range_2_label_asym_id": "B", "range_2_label_seq_id": "1044", "range_2_PDB_ins_code": "?",
            "range_2_auth_atom_id": "O", "range_2_auth_comp_id": "LYS", "range_2_auth_asym_id": "B", "range_2_auth_seq_id": "1044"
        }, {
            "sheet_id": "B", "range_id_1": "2", "range_id_2": "3", "range_1_label_atom_id": "N", "range_1_label_comp_id": "VAL",
            "range_1_label_asym_id": "B", "range_1_label_seq_id": "1041", "range_1_PDB_ins_code": "Q", "range_1_auth_atom_id": "N",
            "range_1_auth_comp_id": "VAL", "range_1_auth_asym_id": "B", "range_1_auth_seq_id": "1041", "range_2_label_atom_id": "O",
            "range_2_label_comp_id": "ARG", "range_2_label_asym_id": "B", "range_2_label_seq_id": "1073", "range_2_PDB_ins_code": "L",
            "range_2_auth_atom_id": "O", "range_2_auth_comp_id": "ARG", "range_2_auth_asym_id": "B", "range_2_auth_seq_id": "1073"
        }]})



class SsbondParsingTests(TestCase):

    def test_can_handle_no_ssbond(self):
        mmcif = {1: 2}
        parse_ssbond("", mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_parse_ssbond(self):
        filestring = (
            "SSBOND   1 CYS A    6Z   CYS A  127                          1555   1555  2.03  \n"
            "SSBOND   2 CYS A   30    CYS A  115P                         1555   1555  2.04  \n"
        )
        mmcif = {1: 2}
        parse_ssbond(filestring, mmcif)
        self.assertEqual(mmcif, {1: 2, "struct_conn": [{
            "id": "disulf1", "conn_type_id": "disulf", "pdbx_leaving_atom_flag": "?",
            "pdbx_PDB_id": "?", "ptnr1_label_asym_id": "A", "ptnr1_label_comp_id": "CYS",
            "ptnr1_label_seq_id": "6", "ptnr1_label_atom_id": "?", "pdbx_ptnr1_label_alt_id": "?",
            "pdbx_ptnr1_PDB_ins_code": "Z", "pdbx_ptnr1_standard_comp_id": "?", "ptnr1_symmetry": "1555",
            "ptnr2_label_asym_id": "A", "ptnr2_label_comp_id": "CYS", "ptnr2_label_seq_id": "127",
            "ptnr2_label_atom_id": "?", "pdbx_ptnr2_label_alt_id": "?", "pdbx_ptnr2_PDB_ins_code": "?",
            "ptnr1_auth_asym_id": "A", "ptnr1_auth_comp_id": "CYS", "ptnr1_auth_seq_id": "6",
            "ptnr2_auth_asym_id": "A", "ptnr2_auth_comp_id": "CYS", "ptnr2_auth_seq_id": "127",
            "ptnr2_symmetry": "1555", "pdbx_ptnr3_label_atom_id": "?", "pdbx_ptnr3_label_seq_id": "?",
            "pdbx_ptnr3_label_comp_id": "?", "pdbx_ptnr3_label_asym_id": "?", "pdbx_ptnr3_label_alt_id": "?",
            "pdbx_ptnr3_PDB_ins_code": "?", "details": "?", "pdbx_dist_value": "2.03", "pdbx_value_order": "?"
        }, {
            "id": "disulf2", "conn_type_id": "disulf", "pdbx_leaving_atom_flag": "?", "pdbx_PDB_id": "?",
            "ptnr1_label_asym_id": "A", "ptnr1_label_comp_id": "CYS", "ptnr1_label_seq_id": "30",
            "ptnr1_label_atom_id": "?", "pdbx_ptnr1_label_alt_id": "?", "pdbx_ptnr1_PDB_ins_code": "?",
            "pdbx_ptnr1_standard_comp_id": "?", "ptnr1_symmetry": "1555", "ptnr2_label_asym_id": "A",
            "ptnr2_label_comp_id": "CYS", "ptnr2_label_seq_id": "115", "ptnr2_label_atom_id": "?",
            "pdbx_ptnr2_label_alt_id": "?", "pdbx_ptnr2_PDB_ins_code": "P", "ptnr1_auth_asym_id": "A",
            "ptnr1_auth_comp_id": "CYS", "ptnr1_auth_seq_id": "30", "ptnr2_auth_asym_id": "A",
            "ptnr2_auth_comp_id": "CYS", "ptnr2_auth_seq_id": "115", "ptnr2_symmetry": "1555",
            "pdbx_ptnr3_label_atom_id": "?", "pdbx_ptnr3_label_seq_id": "?", "pdbx_ptnr3_label_comp_id": "?",
            "pdbx_ptnr3_label_asym_id": "?", "pdbx_ptnr3_label_alt_id": "?", "pdbx_ptnr3_PDB_ins_code": "?",
            "details": "?", "pdbx_dist_value": "2.04", "pdbx_value_order": "?"
        }]})



class LinkParsingTests(TestCase):

    def test_can_handle_no_link(self):
        mmcif = {1: 2}
        parse_link("", mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_parse_link(self):
        filestring = (
            "LINK         OH MOMZ D   2P                C5  GHP D   4     1555   1555  1.39 \n" 
            "LINK         C   ASN D   3                 N  QGHP D   4Z    1555   1555  1.32 \n"
        )
        mmcif = {1: 2}
        parse_link(filestring, mmcif)
        self.assertEqual(mmcif, {1: 2, "struct_conn": [{
            "id": "covale1", "conn_type_id": "covale", "pdbx_leaving_atom_flag": "?", "pdbx_PDB_id": "?",
            "ptnr1_label_asym_id": "D", "ptnr1_label_comp_id": "OMZ", "ptnr1_label_seq_id": "2",
            "ptnr1_label_atom_id": "OH", "pdbx_ptnr1_label_alt_id": "M", "pdbx_ptnr1_PDB_ins_code": "P",
            "pdbx_ptnr1_standard_comp_id": "?", "ptnr1_symmetry": "1555", "ptnr2_label_asym_id": "D",
            "ptnr2_label_comp_id": "GHP", "ptnr2_label_seq_id": "4", "ptnr2_label_atom_id": "C5",
            "pdbx_ptnr2_label_alt_id": "?", "pdbx_ptnr2_PDB_ins_code": "?", "ptnr1_auth_asym_id": "D",
            "ptnr1_auth_comp_id": "OMZ", "ptnr1_auth_seq_id": "2", "ptnr2_auth_asym_id": "D", 
            "ptnr2_auth_comp_id": "GHP", "ptnr2_auth_seq_id": "4", "ptnr2_symmetry": "1555",
            "pdbx_ptnr3_label_atom_id": "?", "pdbx_ptnr3_label_seq_id": "?", "pdbx_ptnr3_label_comp_id": "?",
            "pdbx_ptnr3_label_asym_id": "?", "pdbx_ptnr3_label_alt_id": "?", "pdbx_ptnr3_PDB_ins_code": "?",
            "details": "?", "pdbx_dist_value": "1.39", "pdbx_value_order": "?"
        }, {
            "id": "covale2", "conn_type_id": "covale", "pdbx_leaving_atom_flag": "?", "pdbx_PDB_id": "?",
            "ptnr1_label_asym_id": "D", "ptnr1_label_comp_id": "ASN", "ptnr1_label_seq_id": "3",
            "ptnr1_label_atom_id": "C", "pdbx_ptnr1_label_alt_id": "?", "pdbx_ptnr1_PDB_ins_code": "?",
            "pdbx_ptnr1_standard_comp_id": "?", "ptnr1_symmetry": "1555", "ptnr2_label_asym_id": "D",
            "ptnr2_label_comp_id": "GHP", "ptnr2_label_seq_id": "4", "ptnr2_label_atom_id": "N",
            "pdbx_ptnr2_label_alt_id": "Q", "pdbx_ptnr2_PDB_ins_code": "Z", "ptnr1_auth_asym_id": "D",
            "ptnr1_auth_comp_id": "ASN", "ptnr1_auth_seq_id": "3", "ptnr2_auth_asym_id": "D",
            "ptnr2_auth_comp_id": "GHP", "ptnr2_auth_seq_id": "4", "ptnr2_symmetry": "1555",
            "pdbx_ptnr3_label_atom_id": "?", "pdbx_ptnr3_label_seq_id": "?", "pdbx_ptnr3_label_comp_id": "?",
            "pdbx_ptnr3_label_asym_id": "?", "pdbx_ptnr3_label_alt_id": "?", "pdbx_ptnr3_PDB_ins_code": "?",
            "details": "?", "pdbx_dist_value": "1.32", "pdbx_value_order": "?"
        }]})
    

    def test_can_parse_link_with_existing_category(self):
        filestring = (
            "LINK         OH MOMZ D   2P                C5  GHP D   4     1555   1555  1.39 \n" 
            "LINK         C   ASN D   3                 N  QGHP D   4Z    1555   1555  1.32 \n"
        )
        mmcif = {1: 2, "struct_conn": [3, 4]}
        parse_link(filestring, mmcif)
        self.assertEqual(mmcif, {1: 2, "struct_conn": [3, 4, {
            "id": "covale1", "conn_type_id": "covale", "pdbx_leaving_atom_flag": "?", "pdbx_PDB_id": "?",
            "ptnr1_label_asym_id": "D", "ptnr1_label_comp_id": "OMZ", "ptnr1_label_seq_id": "2",
            "ptnr1_label_atom_id": "OH", "pdbx_ptnr1_label_alt_id": "M", "pdbx_ptnr1_PDB_ins_code": "P",
            "pdbx_ptnr1_standard_comp_id": "?", "ptnr1_symmetry": "1555", "ptnr2_label_asym_id": "D",
            "ptnr2_label_comp_id": "GHP", "ptnr2_label_seq_id": "4", "ptnr2_label_atom_id": "C5",
            "pdbx_ptnr2_label_alt_id": "?", "pdbx_ptnr2_PDB_ins_code": "?", "ptnr1_auth_asym_id": "D",
            "ptnr1_auth_comp_id": "OMZ", "ptnr1_auth_seq_id": "2", "ptnr2_auth_asym_id": "D", 
            "ptnr2_auth_comp_id": "GHP", "ptnr2_auth_seq_id": "4", "ptnr2_symmetry": "1555",
            "pdbx_ptnr3_label_atom_id": "?", "pdbx_ptnr3_label_seq_id": "?", "pdbx_ptnr3_label_comp_id": "?",
            "pdbx_ptnr3_label_asym_id": "?", "pdbx_ptnr3_label_alt_id": "?", "pdbx_ptnr3_PDB_ins_code": "?",
            "details": "?", "pdbx_dist_value": "1.39", "pdbx_value_order": "?"
        }, {
            "id": "covale2", "conn_type_id": "covale", "pdbx_leaving_atom_flag": "?", "pdbx_PDB_id": "?",
            "ptnr1_label_asym_id": "D", "ptnr1_label_comp_id": "ASN", "ptnr1_label_seq_id": "3",
            "ptnr1_label_atom_id": "C", "pdbx_ptnr1_label_alt_id": "?", "pdbx_ptnr1_PDB_ins_code": "?",
            "pdbx_ptnr1_standard_comp_id": "?", "ptnr1_symmetry": "1555", "ptnr2_label_asym_id": "D",
            "ptnr2_label_comp_id": "GHP", "ptnr2_label_seq_id": "4", "ptnr2_label_atom_id": "N",
            "pdbx_ptnr2_label_alt_id": "Q", "pdbx_ptnr2_PDB_ins_code": "Z", "ptnr1_auth_asym_id": "D",
            "ptnr1_auth_comp_id": "ASN", "ptnr1_auth_seq_id": "3", "ptnr2_auth_asym_id": "D",
            "ptnr2_auth_comp_id": "GHP", "ptnr2_auth_seq_id": "4", "ptnr2_symmetry": "1555",
            "pdbx_ptnr3_label_atom_id": "?", "pdbx_ptnr3_label_seq_id": "?", "pdbx_ptnr3_label_comp_id": "?",
            "pdbx_ptnr3_label_asym_id": "?", "pdbx_ptnr3_label_alt_id": "?", "pdbx_ptnr3_PDB_ins_code": "?",
            "details": "?", "pdbx_dist_value": "1.32", "pdbx_value_order": "?"
        }]})



class SiteParsingTests(TestCase):

    def test_can_handle_no_site(self):
        mmcif = {1: 2}
        parse_site("", mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_parse_site(self):
        filestring = (
            "SITE     1 AC1  4 ASP A  70  LYS A  72P LEU A 123  VAL A 155    \n"
            "SITE     1 AC2  5 LYS B1042  ASP B1070  LYS B1072  ILE B1096    \n"
            "SITE     2 AC2  5 VAL B1155                                     \n"
        )
        mmcif = {1: 2}
        parse_site(filestring, mmcif)
        self.assertEqual(mmcif, {1: 2, "struct_site_gen": [{
            "id": "1", "site_id": "AC1", "pdbx_num_res": "4", "label_comp_id": "ASP", "label_asym_id": "A",
            "label_seq_id": "70", "pdbx_auth_ins_code": "?", "auth_comp_id": "ASP", "auth_asym_id": "A",
            "auth_seq_id": "70", "label_atom_id": ".", "label_alt_id": "?", "symmetry": "1_555", "details": "?"
        }, {
            "id": "2", "site_id": "AC1", "pdbx_num_res": "4", "label_comp_id": "LYS", "label_asym_id": "A",
            "label_seq_id": "72", "pdbx_auth_ins_code": "P", "auth_comp_id": "LYS", "auth_asym_id": "A",
            "auth_seq_id": "72", "label_atom_id": ".", "label_alt_id": "?", "symmetry": "1_555", "details": "?"
        }, {
            "id": "3", "site_id": "AC1", "pdbx_num_res": "4", "label_comp_id": "LEU", "label_asym_id": "A",
            "label_seq_id": "123", "pdbx_auth_ins_code": "?", "auth_comp_id": "LEU", "auth_asym_id": "A",
            "auth_seq_id": "123", "label_atom_id": ".", "label_alt_id": "?", "symmetry": "1_555", "details": "?"
        }, {
            "id": "4", "site_id": "AC1", "pdbx_num_res": "4", "label_comp_id": "VAL", "label_asym_id": "A",
            "label_seq_id": "155", "pdbx_auth_ins_code": "?", "auth_comp_id": "VAL", "auth_asym_id": "A",
            "auth_seq_id": "155", "label_atom_id": ".", "label_alt_id": "?", "symmetry": "1_555", "details": "?"
        }, {
            "id": "5", "site_id": "AC2", "pdbx_num_res": "5", "label_comp_id": "LYS", "label_asym_id": "B",
            "label_seq_id": "1042", "pdbx_auth_ins_code": "?", "auth_comp_id": "LYS", "auth_asym_id": "B",
            "auth_seq_id": "1042", "label_atom_id": ".", "label_alt_id": "?", "symmetry": "1_555", "details": "?"
        }, {
            "id": "6", "site_id": "AC2", "pdbx_num_res": "5", "label_comp_id": "ASP", "label_asym_id": "B",
            "label_seq_id": "1070", "pdbx_auth_ins_code": "?", "auth_comp_id": "ASP", "auth_asym_id": "B",
            "auth_seq_id": "1070", "label_atom_id": ".", "label_alt_id": "?", "symmetry": "1_555", "details": "?"
        }, {
            "id": "7", "site_id": "AC2", "pdbx_num_res": "5", "label_comp_id": "LYS", "label_asym_id": "B",
            "label_seq_id": "1072", "pdbx_auth_ins_code": "?", "auth_comp_id": "LYS", "auth_asym_id": "B",
            "auth_seq_id": "1072", "label_atom_id": ".", "label_alt_id": "?", "symmetry": "1_555", "details": "?"
        }, {
            "id": "8", "site_id": "AC2", "pdbx_num_res": "5", "label_comp_id": "ILE", "label_asym_id": "B",
            "label_seq_id": "1096", "pdbx_auth_ins_code": "?", "auth_comp_id": "ILE", "auth_asym_id": "B",
            "auth_seq_id": "1096", "label_atom_id": ".", "label_alt_id": "?", "symmetry": "1_555", "details": "?"
        }, {
            "id": "9", "site_id": "AC2", "pdbx_num_res": "5", "label_comp_id": "VAL", "label_asym_id": "B",
            "label_seq_id": "1155", "pdbx_auth_ins_code": "?", "auth_comp_id": "VAL", "auth_asym_id": "B",
            "auth_seq_id": "1155", "label_atom_id": ".", "label_alt_id": "?", "symmetry": "1_555", "details": "?"
        }]})



class CispepParsingTests(TestCase):

    def test_can_handle_no_cicpep(self):
        mmcif = {1: 2}
        parse_cispep("", mmcif)
        self.assertEqual(mmcif, {1: 2})
    

    def test_can_parse_cispep(self):
        filestring = (
            "CISPEP   1 SER A    7P   PRO A    8          0        -2.43     \n"                
            "CISPEP   2 THR A   94    PRO A   95Q         0        -3.63     \n"
        )
        mmcif = {1: 2}
        parse_cispep(filestring, mmcif)
        self.assertEqual(mmcif, {1: 2, "struct_mon_prot_cis": [{
            "pdbx_id": "1", "label_comp_id": "SER", "label_seq_id": "7", "label_asym_id": "A",
            "label_alt_id": ".", "pdbx_PDB_ins_code": "P", "auth_comp_id": "SER", "auth_seq_id": "7", 
            "auth_asym_id": "A", "pdbx_label_comp_id_2": "PRO", "pdbx_label_seq_id_2": "8",
            "pdbx_label_asym_id_2": "A", "pdbx_PDB_ins_code_2": "?", "pdbx_auth_comp_id_2": "PRO",
            "pdbx_auth_seq_id_2": "8", "pdbx_auth_asym_id_2": "A", "pdbx_PDB_model_num": "1",
            "pdbx_omega_angle": "-2.43"
        }, {
            "pdbx_id": "2", "label_comp_id": "THR", "label_seq_id": "94", "label_asym_id": "A",
            "label_alt_id": ".", "pdbx_PDB_ins_code": "?", "auth_comp_id": "THR", "auth_seq_id": "94",
            "auth_asym_id": "A", "pdbx_label_comp_id_2": "PRO", "pdbx_label_seq_id_2": "95",
            "pdbx_label_asym_id_2": "A", "pdbx_PDB_ins_code_2": "Q", "pdbx_auth_comp_id_2": "PRO",
            "pdbx_auth_seq_id_2": "95", "pdbx_auth_asym_id_2": "A", "pdbx_PDB_model_num": "1",
            "pdbx_omega_angle": "-3.63"
        }]})



class AuthLabelUpdatingTests(TestCase):

    @patch("atomium.pdb.update_auth_ids_in_pdbx_struct_assembly_gen")
    def test_can_update_labels(self, mock_auth):
        mmcif = {
            "atom_site": [
                {"auth_asym_id": "P", "auth_seq_id": "101", "pdbx_PDB_ins_code": "?", "label_asym_id": "A", "label_seq_id": "1"},
                {"auth_asym_id": "P", "auth_seq_id": "102", "pdbx_PDB_ins_code": "?", "label_asym_id": "A", "label_seq_id": "2"},
                {"auth_asym_id": "P", "auth_seq_id": "103", "pdbx_PDB_ins_code": "?", "label_asym_id": "A", "label_seq_id": "3"},
                {"auth_asym_id": "P", "auth_seq_id": "103", "pdbx_PDB_ins_code": "A", "label_asym_id": "A", "label_seq_id": "4"},
                {"auth_asym_id": "P", "auth_seq_id": "103", "pdbx_PDB_ins_code": "B", "label_asym_id": "A", "label_seq_id": "5"},
                {"auth_asym_id": "Q", "auth_seq_id": "201", "pdbx_PDB_ins_code": "?", "label_asym_id": "B", "label_seq_id": "6"},
                {"auth_asym_id": "Q", "auth_seq_id": "202", "pdbx_PDB_ins_code": "?", "label_asym_id": "B", "label_seq_id": "7"},
                {"auth_asym_id": "Q", "auth_seq_id": "203", "pdbx_PDB_ins_code": "?", "label_asym_id": "B", "label_seq_id": "8"},
                {"auth_asym_id": "P", "auth_seq_id": "1000", "pdbx_PDB_ins_code": "?", "label_asym_id": "C", "label_seq_id": "."},
                {"auth_asym_id": "P", "auth_seq_id": "2000", "pdbx_PDB_ins_code": "?", "label_asym_id": "C", "label_seq_id": "."},
                {"auth_asym_id": "P", "auth_seq_id": "3000", "pdbx_PDB_ins_code": "?", "label_asym_id": "C", "label_seq_id": "."},
            ],
            "category1": [
                {"auth_asym_id": "P", "auth_seq_id": "101", "pdbx_PDB_ins_code": "?", "label_asym_id": "P", "label_seq_id": "101"},
                {"auth_asym_id": "P", "auth_seq_id": "103", "pdbx_PDB_ins_code": "A", "label_asym_id": "P", "label_seq_id": "104"},
                {"auth_asym_id": "Z", "auth_seq_id": "103", "pdbx_PDB_ins_code": "A", "label_asym_id": "P", "label_seq_id": "104"},
                {"auth_asym_id": "Q", "auth_seq_id": "203", "pdbx_PDB_ins_code": "?", "label_asym_id": "Q", "label_seq_id": "108"},
            ],
            "category2": [{
                "auth_asym_id_1": "P", "auth_seq_id_1": "101", "pdbx_PDB_ins_code_1": "?", "label_asym_id_1": "P", "label_seq_id_1": "101",
                "auth_asym_id_2": "Q", "auth_seq_id_2": "201", "pdbx_PDB_ins_code_2": "?", "label_asym_id_2": "Q", "label_seq_id_2": "201",
            }],
            "category3": [{
                "pre_auth_asym_id": "P", "pre_auth_seq_id": "101", "pre_pdbx_PDB_ins_code": "?", "pre_label_asym_id": "P", "pre_label_seq_id": "101",
                "post_auth_asym_id": "Q", "post_auth_seq_id": "201", "post_pdbx_PDB_ins_code": "?", "post_label_asym_id": "Q", "post_label_seq_id": "201",
            }]
        }
        update_auth_and_label(mmcif)
        mock_auth.assert_called_with(mmcif)
        self.assertEqual(mmcif, {
            "atom_site": [
                {"auth_asym_id": "P", "auth_seq_id": "101", "pdbx_PDB_ins_code": "?", "label_asym_id": "A", "label_seq_id": "1"},
                {"auth_asym_id": "P", "auth_seq_id": "102", "pdbx_PDB_ins_code": "?", "label_asym_id": "A", "label_seq_id": "2"},
                {"auth_asym_id": "P", "auth_seq_id": "103", "pdbx_PDB_ins_code": "?", "label_asym_id": "A", "label_seq_id": "3"},
                {"auth_asym_id": "P", "auth_seq_id": "103", "pdbx_PDB_ins_code": "A", "label_asym_id": "A", "label_seq_id": "4"},
                {"auth_asym_id": "P", "auth_seq_id": "103", "pdbx_PDB_ins_code": "B", "label_asym_id": "A", "label_seq_id": "5"},
                {"auth_asym_id": "Q", "auth_seq_id": "201", "pdbx_PDB_ins_code": "?", "label_asym_id": "B", "label_seq_id": "6"},
                {"auth_asym_id": "Q", "auth_seq_id": "202", "pdbx_PDB_ins_code": "?", "label_asym_id": "B", "label_seq_id": "7"},
                {"auth_asym_id": "Q", "auth_seq_id": "203", "pdbx_PDB_ins_code": "?", "label_asym_id": "B", "label_seq_id": "8"},
                {"auth_asym_id": "P", "auth_seq_id": "1000", "pdbx_PDB_ins_code": "?", "label_asym_id": "C", "label_seq_id": "."},
                {"auth_asym_id": "P", "auth_seq_id": "2000", "pdbx_PDB_ins_code": "?", "label_asym_id": "C", "label_seq_id": "."},
                {"auth_asym_id": "P", "auth_seq_id": "3000", "pdbx_PDB_ins_code": "?", "label_asym_id": "C", "label_seq_id": "."},
            ],
            "category1": [
                {"auth_asym_id": "P", "auth_seq_id": "101", "pdbx_PDB_ins_code": "?", "label_asym_id": "A", "label_seq_id": "1"},
                {"auth_asym_id": "P", "auth_seq_id": "103", "pdbx_PDB_ins_code": "A", "label_asym_id": "A", "label_seq_id": "4"},
                {"auth_asym_id": "Z", "auth_seq_id": "103", "pdbx_PDB_ins_code": "A", "label_asym_id": "P", "label_seq_id": "104"},
                {"auth_asym_id": "Q", "auth_seq_id": "203", "pdbx_PDB_ins_code": "?", "label_asym_id": "B", "label_seq_id": "8"},
            ],
            "category2": [{
                "auth_asym_id_1": "P", "auth_seq_id_1": "101", "pdbx_PDB_ins_code_1": "?", "label_asym_id_1": "A", "label_seq_id_1": "1",
                "auth_asym_id_2": "Q", "auth_seq_id_2": "201", "pdbx_PDB_ins_code_2": "?", "label_asym_id_2": "B", "label_seq_id_2": "6",
            }],
            "category3": [{
                "pre_auth_asym_id": "P", "pre_auth_seq_id": "101", "pre_pdbx_PDB_ins_code": "?", "pre_label_asym_id": "A", "pre_label_seq_id": "1",
                "post_auth_asym_id": "Q", "post_auth_seq_id": "201", "post_pdbx_PDB_ins_code": "?", "post_label_asym_id": "B", "post_label_seq_id": "6",
            }]
        })



class AuthIdsInStructAssemblyUpdatingTests(TestCase):

    def test_can_handle_no_categories(self):
        update_auth_ids_in_pdbx_struct_assembly_gen({})
    

    def test_can_replace_ids(self):
        mmcif = {
            "pdbx_struct_assembly_gen": [{"asym_id_list": "P,B"}, {"asym_id_list": "C"}],
            "atom_site": [
                {"auth_asym_id": "P", "label_asym_id": "A"},
                {"auth_asym_id": "P", "label_asym_id": "A"},
                {"auth_asym_id": "B", "label_asym_id": "B"},
                {"auth_asym_id": "B", "label_asym_id": "B"},
                {"auth_asym_id": "C", "label_asym_id": "C"},
                {"auth_asym_id": "P", "label_asym_id": "D"},
                {"auth_asym_id": "B", "label_asym_id": "E"},
                {"auth_asym_id": "C", "label_asym_id": "F"},
            ]
        }
        update_auth_ids_in_pdbx_struct_assembly_gen(mmcif)
        self.assertEqual(mmcif["pdbx_struct_assembly_gen"], [
            {"asym_id_list": "A,B,D,E"}, {"asym_id_list": "C,F"}
        ])



class PdbDateToMmcifDateTests(TestCase):

    def test_can_convert_date(self):
        self.assertEqual(pdb_date_to_mmcif_date("06-MAY-02"), "2002-05-06")
        self.assertEqual(pdb_date_to_mmcif_date("19-SEP-76"), "1976-09-19")
    

    def test_can_handle_invalid(self):
        self.assertIsNone(pdb_date_to_mmcif_date("  "))



class PdbNamesToMmcifNamesTests(TestCase):

    def test_single_name(self):
        lines = ["N.WU"]
        self.assertEqual(pdb_names_to_mmcif_names(lines), ["Wu, N"])
    

    def test_multiple_names(self):
        lines = ["J.L.WHITTINGHAM,S.HAVELUND,I.JONASSEN"]
        self.assertEqual(pdb_names_to_mmcif_names(lines), ["Whittingham, J.L", "Havelund, S", "Jonassen, I"])
    

    def test_multiple_lines(self):
        lines = [
            "J.MARKUSSEN,S.HAVELUND,P.KURTZHALS,A.S.ANDERSEN,             ",
            "J.HALSTROM,E.HASSELAGER,U.D.LARSEN,U.RIBEL,                  ",
            "L.SCHAFFER,K.VAD,I.JONASSEN        "
        ]
        self.assertEqual(pdb_names_to_mmcif_names(lines), [
            "Markussen, J", "Havelund, S", "Kurtzhals, P", "Andersen, A.S", "Halstrom, J",
            "Hasselager, E", "Larsen, U.D", "Ribel, U", "Schaffer, L", "Vad, K", "Jonassen, I"
        ])



class MmcifDictSavingTests(TestCase):

    @patch("atomium.pdb.create_header_line")
    @patch("atomium.pdb.create_obslte_lines")
    @patch("atomium.pdb.create_title_lines")
    @patch("atomium.pdb.create_split_lines")
    @patch("atomium.pdb.create_caveat_lines")
    @patch("atomium.pdb.create_compnd_lines")
    @patch("atomium.pdb.create_source_lines")
    @patch("atomium.pdb.create_keywds_lines")
    @patch("atomium.pdb.create_expdta_lines")
    @patch("atomium.pdb.create_nummdl_lines")
    @patch("atomium.pdb.create_mdltyp_lines")
    @patch("atomium.pdb.create_author_lines")
    @patch("atomium.pdb.create_revdat_lines")
    @patch("atomium.pdb.create_sprsde_lines")
    @patch("atomium.pdb.create_jrnl_lines")
    @patch("atomium.pdb.create_seqadv_lines")
    @patch("atomium.pdb.create_seqres_lines")
    @patch("atomium.pdb.create_modres_lines")
    @patch("atomium.pdb.create_het_lines")
    @patch("atomium.pdb.create_hetnam_lines")
    @patch("atomium.pdb.create_hetsyn_lines")
    @patch("atomium.pdb.create_formul_lines")
    @patch("atomium.pdb.create_helix_lines")
    @patch("atomium.pdb.create_sheet_lines")
    @patch("atomium.pdb.create_ssbond_lines")
    @patch("atomium.pdb.create_link_lines")
    @patch("atomium.pdb.create_cispep_lines")
    @patch("atomium.pdb.create_site_lines")
    @patch("atomium.pdb.create_cryst1_line")
    @patch("atomium.pdb.create_origxn_lines")
    @patch("atomium.pdb.create_scalen_lines")
    @patch("atomium.pdb.create_mtrixn_lines")
    @patch("atomium.pdb.create_atom_lines")
    @patch("builtins.open")
    def test_can_save_dict(self, mock_open, *mocks):
        mmcif = {"mmcif": 1}
        for i, mock in enumerate(mocks[::-1]):
            mock.return_value = [str(i)]
        save_mmcif_dict(mmcif, "/path")
        mock_open.assert_called_with("/path", "w")
        mock_open.return_value.__enter__.return_value.write.assert_called_with(
            "\n".join(map(str, range(33))) + "\nEND"
        )
        for mock in mocks:
            mock.assert_called_with(mmcif)



class HeaderLineSavingTests(TestCase):

    @patch("atomium.pdb.create_pdb_date")
    def test_can_handle_no_categories(self, mock_date):
        mock_date.return_value = None
        self.assertEqual(create_header_line({}), [])
        mock_date.assert_called_with("")
    

    @patch("atomium.pdb.create_pdb_date")
    def test_can_handle_no_values(self, mock_date):
        mock_date.return_value = None
        self.assertEqual(create_header_line({
            "entry": [{"id": "?"}],
            "struct_keywords": [{"pdbx_keywords": "?"}],
            "pdbx_database_status": [{"recvd_initial_deposition_date": "?"}],
        }), [])
        self.assertFalse(mock_date.called)
    

    @patch("atomium.pdb.create_pdb_date")
    def test_can_parse_code(self, mock_date):
        mock_date.return_value = None
        self.assertEqual(create_header_line({
            "entry": [{"id": "1XXX"}],
            "struct_keywords": [{"pdbx_keywords": "?"}],
            "pdbx_database_status": [{"recvd_initial_deposition_date": "?"}],
        }), ["HEADER                                                        1XXX"])
        self.assertFalse(mock_date.called)
    

    @patch("atomium.pdb.create_pdb_date")
    def test_can_parse_keywords(self, mock_date):
        mock_date.return_value = None
        self.assertEqual(create_header_line({
            "entry": [{"id": "?"}],
            "struct_keywords": [{"pdbx_keywords": "APOPTOSIS"}],
            "pdbx_database_status": [{"recvd_initial_deposition_date": "?"}],
        }), ["HEADER    APOPTOSIS                                               "])
        self.assertFalse(mock_date.called)
    

    @patch("atomium.pdb.create_pdb_date")
    def test_can_parse_date(self, mock_date):
        mock_date.return_value = "02-FEB-02"
        self.assertEqual(create_header_line({
            "entry": [{"id": "?"}],
            "struct_keywords": [{"pdbx_keywords": "?"}],
            "pdbx_database_status": [{"recvd_initial_deposition_date": "2002-01-01"}],
        }), ["HEADER                                            02-FEB-02       "])
        mock_date.assert_called_with("2002-01-01")
    

    @patch("atomium.pdb.create_pdb_date")
    def test_can_parse_all_values(self, mock_date):
        mock_date.return_value = "02-FEB-02"
        self.assertEqual(create_header_line({
            "entry": [{"id": "1XXX"}],
            "struct_keywords": [{"pdbx_keywords": "-123456789" * 5}],
            "pdbx_database_status": [{"recvd_initial_deposition_date": "2002-01-01"}],
        }), ["HEADER    -123456789-123456789-123456789-12345678902-FEB-02   1XXX"])
        mock_date.assert_called_with("2002-01-01")



class ObslteLinesSavingTests(TestCase):

    def test_can_handle_no_category(self):
        self.assertEqual(create_obslte_lines({}), [])
    

    def test_can_handle_no_obslte(self):
        self.assertEqual(create_obslte_lines({
            "pdbx_database_PDB_obs_spr": [
                {"id": "SPRSDE", "pdb_id": "2XXX", "replace_pdb_id": "1XXX", "date": "2002"},
            ]
        }), [])
    

    @patch("atomium.pdb.create_pdb_date")
    def test_can_create_obslte(self, mock_date):
        mock_date.return_value = "25-JUL-06"
        self.assertEqual(create_obslte_lines({
            "pdbx_database_PDB_obs_spr": [
                {"id": "SPRSDE", "pdb_id": "XXXX", "replace_pdb_id": "YXXX", "date": "2002"},
                {"id": "OBSLTE", "pdb_id": "2XXX", "replace_pdb_id": "1XXX", "date": "2002"},
            ]
        }), ["OBSLTE     25-JUL-06 1XXX      2XXX"])
        mock_date.assert_called_with("2002")



class TitleLinesSavingTests(TestCase):

    def test_can_handle_no_struct(self):
        self.assertEqual(create_title_lines({}), [])
    

    def test_can_handle_no_title(self):
        mmcif = {"struct": [{"title": "?"}]}
        self.assertEqual(create_title_lines(mmcif), [])
    

    @patch("atomium.pdb.split_lines")
    def test_can_save_single_line_title(self, mock_split):
        mmcif = {"struct": [{"title": "TTT"}]}
        mock_split.return_value = ["The title"]
        self.assertEqual(create_title_lines(mmcif), ["TITLE     THE TITLE"])
        mock_split.assert_called_with("TTT", 60)
    

    @patch("atomium.pdb.split_lines")
    def test_can_save_single_multi_line_title(self, mock_split):
        mmcif = {"struct": [{"title": "TTT"}]}
        mock_split.return_value = ["The title", "and the rest", "of the title"]
        self.assertEqual(create_title_lines(mmcif), [
            "TITLE     THE TITLE",
            "TITLE    2 AND THE REST",
            "TITLE    3 OF THE TITLE",
        ])
        mock_split.assert_called_with("TTT", 60)



class SplitLinesSavingTests(TestCase):

    def test_can_handle_no_category(self):
        self.assertEqual(create_split_lines({}), [])
    

    def test_can_get_single_line(self):
        mmcif = {"pdbx_database_related": [
            {"db_id": "1XXX"}, {"db_id": "2XXX"}, {"db_id": "3XXX"}
        ]}
        self.assertEqual(create_split_lines(mmcif), [
            "SPLIT      1XXX 2XXX 3XXX"
        ])
    

    def test_can_get_multiple_lines(self):
        mmcif = {"pdbx_database_related": [
            {"db_id": "1VUU"}, {"db_id": "1VUV"}, {"db_id": "1VUW"},
            {"db_id": "1VUX"}, {"db_id": "1VUY"}, {"db_id": "1VUZ"},
            {"db_id": "1VV0"}, {"db_id": "1VV1"}, {"db_id": "1VV2"},
            {"db_id": "1VV3"}, {"db_id": "1VV4"}, {"db_id": "1VV5"},
            {"db_id": "1VV6"}, {"db_id": "1VV7"}, {"db_id": "1VV8"},
            {"db_id": "1VV9"}, {"db_id": "1VVA"}, {"db_id": "1VVB"},
        ]}
        self.assertEqual(create_split_lines(mmcif), [
            "SPLIT      1VUU 1VUV 1VUW 1VUX 1VUY 1VUZ 1VV0 1VV1 1VV2 1VV3 1VV4 1VV5 1VV6 1VV7",
            "SPLIT    2 1VV8 1VV9 1VVA 1VVB"
        ])



class CaveatLinesSavingTests(TestCase):

    def test_can_handle_no_category(self):
        self.assertEqual(create_caveat_lines({}), [])
    

    @patch("atomium.pdb.split_lines")
    def test_can_save_single_row_to_single_line(self, mock_split):
        mock_split.return_value = ["A CAVEAT"]
        mmcif = {"database_PDB_caveat": [{"text": "A caveat"}], "entry": [{"id": "1XXX"}]}
        self.assertEqual(create_caveat_lines(mmcif), [
            "CAVEAT     1XXX    A CAVEAT",
        ])
        mock_split.assert_called_with("A CAVEAT", 60)
    

    @patch("atomium.pdb.split_lines")
    def test_can_save_multiple_rows_to_single_line(self, mock_split):
        mock_split.return_value = ["A CAVEAT ANOTHER CAVEAT"]
        mmcif = {
            "database_PDB_caveat": [{"text": "A caveat"}, {"text": "Another caveat"}],
            "entry": [{"id": "1XXX"}]
        }
        self.assertEqual(create_caveat_lines(mmcif), [
            "CAVEAT     1XXX    A CAVEAT ANOTHER CAVEAT",
        ])
        mock_split.assert_called_with("A CAVEAT ANOTHER CAVEAT", 60)
    

    @patch("atomium.pdb.split_lines")
    def test_can_save_single_row_to_multiple_lines(self, mock_split):
        mock_split.return_value = ["A CAVEAT", "ANOTHER"]
        mmcif = {"database_PDB_caveat": [{"text": "Aaah"}], "entry": [{"id": "1XXX"}]}
        self.assertEqual(create_caveat_lines(mmcif), [
            "CAVEAT     1XXX    A CAVEAT",
            "CAVEAT   2 1XXX    ANOTHER",
        ])
        mock_split.assert_called_with("AAAH", 60)
    

    @patch("atomium.pdb.split_lines")
    def test_can_save_multiple_rows_to_multiple_lines(self, mock_split):
        mock_split.return_value = ["A CAVEAT", "ANOTHER"]
        mmcif = {"database_PDB_caveat": [{"text": "wow"}, {"text": "Aaah"}], "entry": [{"id": "1XXX"}]}
        self.assertEqual(create_caveat_lines(mmcif), [
            "CAVEAT     1XXX    A CAVEAT",
            "CAVEAT   2 1XXX    ANOTHER",
        ])
        mock_split.assert_called_with("WOW AAAH", 60)



class CompndLinesTests(TestCase):

    @patch("atomium.pdb.split_lines")
    def test_can_produce_minimal_compnd(self, mock_split):
        mmcif = {
            "entity": [
                {"id": "1", "type": "polymer", "pdbx_description": "?", "pdbx_ec": "?", "src_method": "nat"},
                {"id": "2", "type": "polymer", "pdbx_description": "?", "pdbx_ec": "?", "src_method": "syn"},
                {"id": "3", "type": "non-polymer"},
                {"id": "4", "type": "water"},
            ],
            "struct_asym": [
                {"id": "A", "entity_id": "1"},
                {"id": "B", "entity_id": "1"},
                {"id": "C", "entity_id": "1"},
                {"id": "D", "entity_id": "2"},
                {"id": "E", "entity_id": "2"},
                {"id": "F", "entity_id": "3"},
                {"id": "G", "entity_id": "4"},
            ],
            "atom_site": [
                {"label_asym_id": "A", "auth_asym_id": "T"},
                {"label_asym_id": "B", "auth_asym_id": "U"},
                {"label_asym_id": "C", "auth_asym_id": "V"},
                {"label_asym_id": "D", "auth_asym_id": "X"},
                {"label_asym_id": "E", "auth_asym_id": "Y"},
                {"label_asym_id": "F", "auth_asym_id": "T"},
                {"label_asym_id": "G", "auth_asym_id": "U"},
            ]
        }
        mock_split.side_effect = lambda l, n: [l]
        lines = create_compnd_lines(mmcif)
        self.assertEqual([c[0] for c in mock_split.call_args_list], [
            ("MOL_ID: 1;", 70),
            ("CHAIN: T, U, V;", 70),
            ("MOL_ID: 2;", 70),
            ("CHAIN: X, Y;", 70),
        ])
        self.assertEqual(lines, [
            "COMPND    MOL_ID: 1;",
            "COMPND   2 CHAIN: T, U, V;",
            "COMPND   3 MOL_ID: 2;",
            "COMPND   4 CHAIN: X, Y;"
        ])


    @patch("atomium.pdb.split_lines")
    def test_can_produce_full_compnd(self, mock_split):
        mmcif = {
            "entity": [
                {"id": "1", "type": "polymer", "pdbx_description": "Insulin", "pdbx_ec": "1.2.3", "src_method": "man"},
                {"id": "2", "type": "polymer", "pdbx_description": "Kinesin", "pdbx_ec": "4.5.6", "src_method": "syn"},
                {"id": "3", "type": "non-polymer"},
                {"id": "4", "type": "water"},
            ],
            "entity_name_com": [
                {"entity_id": "1", "name": "Asucrium"},
                {"entity_id": "2", "name": "Springy"},
            ],
            "struct_asym": [
                {"id": "A", "entity_id": "1"},
                {"id": "B", "entity_id": "1"},
                {"id": "C", "entity_id": "1"},
                {"id": "D", "entity_id": "2"},
                {"id": "E", "entity_id": "2"},
                {"id": "F", "entity_id": "3"},
                {"id": "G", "entity_id": "4"},
            ],
            "atom_site": [
                {"label_asym_id": "A", "auth_asym_id": "T"},
                {"label_asym_id": "B", "auth_asym_id": "U"},
                {"label_asym_id": "C", "auth_asym_id": "V"},
                {"label_asym_id": "D", "auth_asym_id": "X"},
                {"label_asym_id": "E", "auth_asym_id": "Y"},
                {"label_asym_id": "F", "auth_asym_id": "T"},
                {"label_asym_id": "G", "auth_asym_id": "U"},
            ]
        }
        mock_split.side_effect = lambda l, n: [l]
        lines = create_compnd_lines(mmcif)
        self.assertEqual([c[0] for c in mock_split.call_args_list], [
            ("MOL_ID: 1;", 70),
            ("MOLECULE: INSULIN;", 70),
            ("CHAIN: T, U, V;", 70),
            ("SYNONYM: ASUCRIUM;", 70),
            ("EC: 1.2.3;", 70),
            ("ENGINEERED: YES;", 70),
            ("MOL_ID: 2;", 70),
            ("MOLECULE: KINESIN;", 70),
            ("CHAIN: X, Y;", 70),
            ("SYNONYM: SPRINGY;", 70),
            ("EC: 4.5.6;", 70),
        ])
        self.assertEqual(lines, [
            "COMPND    MOL_ID: 1;",
            "COMPND   2 MOLECULE: INSULIN;",
            "COMPND   3 CHAIN: T, U, V;",
            "COMPND   4 SYNONYM: ASUCRIUM;",
            "COMPND   5 EC: 1.2.3;",
            "COMPND   6 ENGINEERED: YES;",
            "COMPND   7 MOL_ID: 2;",
            "COMPND   8 MOLECULE: KINESIN;",
            "COMPND   9 CHAIN: X, Y;",
            "COMPND  10 SYNONYM: SPRINGY;",
            "COMPND  11 EC: 4.5.6;",
        ])



class SourceLinesTests(TestCase):

    @patch("atomium.pdb.split_lines")
    def test_can_produce_full_source(self, mock_split):
        mmcif = {
            "entity": [
                {"id": "1", "type": "polymer", "pdbx_description": "?", "pdbx_ec": "?", "src_method": "syn"},
                {"id": "2", "type": "polymer", "pdbx_description": "?", "pdbx_ec": "?", "src_method": "man"},
                {"id": "3", "type": "polymer", "pdbx_description": "?", "pdbx_ec": "?", "src_method": "syn"},
                {"id": "4", "type": "polymer", "pdbx_description": "?", "pdbx_ec": "?", "src_method": "nat"},
                {"id": "5", "type": "non-polymer"},
                {"id": "6", "type": "water"},
            ]
        }
        mock_split.side_effect = lambda l, n: [l]
        lines = create_source_lines(mmcif)
        self.assertEqual([c[0] for c in mock_split.call_args_list], [
            ("MOL_ID: 1;", 70),
            ("SYNTHETIC: YES;", 70),
            ("MOL_ID: 3;", 70),
            ("SYNTHETIC: YES;", 70),
        ])
        self.assertEqual(lines, [
            "SOURCE    MOL_ID: 1;",
            "SOURCE   2 SYNTHETIC: YES;",
            "SOURCE   3 MOL_ID: 3;",
            "SOURCE   4 SYNTHETIC: YES;"
        ])
    

    @patch("atomium.pdb.split_lines")
    def test_can_produce_no_source(self, mock_split):
        mmcif = {
            "entity": [
                {"id": "1", "type": "polymer", "pdbx_description": "?", "pdbx_ec": "?", "src_method": "man"},
                {"id": "2", "type": "polymer", "pdbx_description": "?", "pdbx_ec": "?", "src_method": "nat"},
                {"id": "3", "type": "non-polymer"},
                {"id": "4", "type": "water"},
            ]
        }
        mock_split.side_effect = lambda l, n: [l]
        lines = create_source_lines(mmcif)
        self.assertEqual([c[0] for c in mock_split.call_args_list], [])
        self.assertEqual(lines, [])
    


class KeywdsLinesTests(TestCase):

    def test_can_handle_no_category(self):
        self.assertEqual(create_keywds_lines({}), [])
    

    def test_can_handle_no_keywds(self):
        self.assertEqual(create_keywds_lines({"struct_keywords": [{"text": "?"}]}), [])
    

    @patch("atomium.pdb.split_lines")
    def test_can_save_to_one_line(self, mock_split):
        mock_split.return_value = ["KEYWORD 1, 2", "KEYWORD 3"]
        mmcif = {"struct_keywords": [{"text": "keywords"}]}
        self.assertEqual(create_keywds_lines(mmcif), [
            "KEYWDS    KEYWORD 1, 2", "KEYWDS   1 KEYWORD 3"
        ])
        mock_split.assert_called_with("KEYWORDS", 69)
    

    @patch("atomium.pdb.split_lines")
    def test_can_save_to_multiple_lines(self, mock_split):
        mock_split.return_value = ["KEYWORD 1, 2"]
        mmcif = {"struct_keywords": [{"text": "keywords"}]}
        self.assertEqual(create_keywds_lines(mmcif), [
            "KEYWDS    KEYWORD 1, 2",
        ])
        mock_split.assert_called_with("KEYWORDS", 69)



class ExpdtaLinesTests(TestCase):

    def test_can_handle_no_table(self):
        self.assertEqual(create_expdta_lines({}), [])
    

    def test_can_handle_no_methods(self):
        mmcif = {"exptl": [{"method": "?"}]}
        self.assertEqual(create_expdta_lines(mmcif), [])
    

    @patch("atomium.pdb.split_lines")
    def test_can_handle_one_methods(self, mock_split):
        mmcif = {"exptl": [{"method": "X-ray"}, {"method": "?"}]}
        mock_split.return_value = ["X-RAY"]
        self.assertEqual(create_expdta_lines(mmcif), [
            "EXPDTA    X-RAY"
        ])
        mock_split.assert_called_with("X-RAY", 69)
    

    @patch("atomium.pdb.split_lines")
    def test_can_handle_multiple_methods(self, mock_split):
        mmcif = {"exptl": [{"method": "X-ray"}, {"method": "NMR"}]}
        mock_split.return_value = ["X-RAY; NMR", "OTHER LINE"]
        self.assertEqual(create_expdta_lines(mmcif), [
            "EXPDTA    X-RAY; NMR", "EXPDTA   1 OTHER LINE"
        ])
        mock_split.assert_called_with("X-RAY; NMR", 69)



class NummdlLinesTests(TestCase):

    def test_can_handle_no_table(self):
        self.assertEqual(create_nummdl_lines({}), ["NUMMDL    0"])
    

    def test_can_handle_one_model(self):
        mmcif = {"atom_site": [
            {"pdbx_PDB_model_num": "1"}, {"pdbx_PDB_model_num": "1"}
        ]}
        self.assertEqual(create_nummdl_lines(mmcif), [])
    

    def test_can_handle_multiple_models(self):
        mmcif = {"atom_site": [
            {"pdbx_PDB_model_num": "1"}, {"pdbx_PDB_model_num": "1"},
            {"pdbx_PDB_model_num": "2"}, {"pdbx_PDB_model_num": "2"},
            {"pdbx_PDB_model_num": "3"}, {"pdbx_PDB_model_num": "3"},
        ]}
        self.assertEqual(create_nummdl_lines(mmcif), ["NUMMDL    3"])



class MdlTypeLinesTests(TestCase):

    def test_can_handle_no_table(self):
        self.assertEqual(create_mdltyp_lines({}), [])
    

    def test_can_handle_no_values(self):
        mmcif = {"struct": [{"pdbx_model_type_details": "?"}]}
        self.assertEqual(create_mdltyp_lines(mmcif), [])
    

    @patch("atomium.pdb.split_lines")
    def test_can_parse_struct_mdltyp(self, mock_split):
        mmcif = {"struct": [{"pdbx_model_type_details": "Minimized"}]}
        mock_split.return_value = ["MINIMIZED"]
        self.assertEqual(create_mdltyp_lines(mmcif), ["MDLTYP    MINIMIZED"])
        mock_split.assert_called_with("MINIMIZED", 60)
    

    @patch("atomium.pdb.split_lines")
    def test_can_parse_coordinate_mdltyp(self, mock_split):
        mmcif = {
            "pdbx_coordinate_model": [
                {"asym_id": "A", "type": "CA atoms ONLY"},
                {"asym_id": "B", "type": "CA atoms ONLY"},
                {"asym_id": "C", "type": "P atoms ONLY"},
                {"asym_id": "D", "type": "CA atoms ONLY"},
                {"asym_id": "E", "type": "P atoms ONLY"},
                {"asym_id": "F", "type": "CA atoms ONLY"},
            ]
        }
        mock_split.return_value = ["MINIMIZED"]
        self.assertEqual(create_mdltyp_lines(mmcif), ["MDLTYP    MINIMIZED"])
        mock_split.assert_called_with("CA ATOMS ONLY, CHAIN A, B, D, F; P ATOMS ONLY, CHAIN C, E", 60)


    @patch("atomium.pdb.split_lines")
    def test_can_parse_full(self, mock_split):
        mmcif = {
            "pdbx_coordinate_model": [
                {"asym_id": "A", "type": "CA ATOMS ONLY"},
                {"asym_id": "B", "type": "CA ATOMS ONLY"},
                {"asym_id": "C", "type": "P ATOMS ONLY"},
                {"asym_id": "D", "type": "CA ATOMS ONLY"},
                {"asym_id": "E", "type": "P ATOMS ONLY"},
                {"asym_id": "F", "type": "CA ATOMS ONLY"},
            ],
            "struct": [{"pdbx_model_type_details": "Minimized"}]
        }
        mock_split.return_value = ["MINIMIZED 1", "MINIMIZED 2", "MINIMIZED 3"]
        self.assertEqual(create_mdltyp_lines(mmcif), [
            "MDLTYP    MINIMIZED 1",
            "MDLTYP   2 MINIMIZED 2",
            "MDLTYP   3 MINIMIZED 3",
        ])
        mock_split.assert_called_with("MINIMIZED; CA ATOMS ONLY, CHAIN A, B, D, F; P ATOMS ONLY, CHAIN C, E", 60)



class AuthorLinesTests(TestCase):

    def test_can_handle_no_table(self):
        self.assertEqual(create_author_lines({}), [])
    

    @patch("atomium.pdb.mmcif_names_to_pdb_names")
    @patch("atomium.pdb.split_lines")
    def test_can_handle_one_name(self, mock_split, mock_names):
        mmcif = {"audit_author": [{"name": "John"}]}
        mock_names.return_value = "X,Y,Z"
        mock_split.return_value = ["NAMES1"]
        self.assertEqual(create_author_lines(mmcif), ["AUTHOR    NAMES1"])
        mock_names.assert_called_with(["John"])
        mock_split.assert_called_with("X,Y,Z", 65)
    

    @patch("atomium.pdb.mmcif_names_to_pdb_names")
    @patch("atomium.pdb.split_lines")
    def test_can_handle_multi_name(self, mock_split, mock_names):
        mmcif = {"audit_author": [{"name": "John"}, {"name": "Flo"}]}
        mock_names.return_value = "X,Y,Z"
        mock_split.return_value = ["NAMES1", "NAMES2", "NAMES3"]
        self.assertEqual(create_author_lines(mmcif), [
            "AUTHOR    NAMES1", "AUTHOR   2 NAMES2", "AUTHOR   3 NAMES3"
        ])
        mock_names.assert_called_with(["John", "Flo"])
        mock_split.assert_called_with("X,Y,Z", 65)



class RevdatLinesTests(TestCase):

    def test_can_handle_no_table(self):
        self.assertEqual(create_revdat_lines({}), [])
    

    @patch("atomium.pdb.create_pdb_date")
    def test_can_parse_table(self, mock_date):
        mock_date.side_effect = ["01-JAN-02", "02-JUN-02", "03-FEB-03"]
        mmcif = {
            "entry": [{"id": "1XXX"}],
            "pdbx_audit_revision_history": [
                {"revision_date": "2001"}, {"revision_date": "2002"}, {"revision_date": "2003"}
            ]
        }
        self.assertEqual(create_revdat_lines(mmcif), [
            "REVDAT   3   03-FEB-03 1XXX",
            "REVDAT   2   02-JUN-02 1XXX",
            "REVDAT   1   01-JAN-02 1XXX",
        ])
    

    @patch("atomium.pdb.create_pdb_date")
    def test_can_handle_no_id(self, mock_date):
        mock_date.side_effect = ["01-JAN-02", "02-JUN-02", "03-FEB-03"]
        mmcif = {
            "pdbx_audit_revision_history": [
                {"revision_date": "2001"}, {"revision_date": "2002"}, {"revision_date": "2003"}
            ]
        }
        self.assertEqual(create_revdat_lines(mmcif), [
            "REVDAT   3   03-FEB-03     ",
            "REVDAT   2   02-JUN-02     ",
            "REVDAT   1   01-JAN-02     ",
        ])



class SprsdeLinesTests(TestCase):

    def test_can_handle_no_category(self):
        self.assertEqual(create_sprsde_lines({}), [])
    

    def test_can_handle_no_sprsde(self):
        self.assertEqual(create_sprsde_lines({
            "pdbx_database_PDB_obs_spr": [
                {"id": "OBSLTE", "pdb_id": "2XXX", "replace_pdb_id": "1XXX", "date": "2002"},
            ]
        }), [])
    

    @patch("atomium.pdb.create_pdb_date")
    def test_can_create_obslte(self, mock_date):
        mock_date.return_value = "25-JUL-06"
        self.assertEqual(create_sprsde_lines({
            "pdbx_database_PDB_obs_spr": [
                {"id": "OBSLTE", "pdb_id": "XXXX", "replace_pdb_id": "YXXX", "date": "2002"},
                {"id": "SPRSDE", "pdb_id": "2XXX", "replace_pdb_id": "1XXX", "date": "2002"},
            ]
        }), ["SPRSDE     25-JUL-06 1XXX      2XXX"])
        mock_date.assert_called_with("2002")



class JrnlLinesTests(TestCase):

    @patch("atomium.pdb.create_jrnl_auth_lines")
    @patch("atomium.pdb.create_jrnl_titl_lines")
    @patch("atomium.pdb.create_jrnl_edit_lines")
    @patch("atomium.pdb.create_jrnl_ref_lines")
    @patch("atomium.pdb.create_jrnl_publ_lines")
    @patch("atomium.pdb.create_jrnl_refn_lines")
    @patch("atomium.pdb.create_jrnl_pmid_lines")
    @patch("atomium.pdb.create_jrnl_doi_lines")
    def test_can_produce_jrnl_lines(self, *mocks):
        mmcif = {"mmcif": 1}
        for i, mock in enumerate(mocks[::-1]):
            mock.return_value = [i]
        lines = create_jrnl_lines(mmcif)
        self.assertEqual(lines, [0, 1, 2, 3, 4, 5, 6, 7])



class JrnlAuthLinesTests(TestCase):

    def test_can_handle_no_table(self):
        self.assertEqual(create_jrnl_auth_lines({}), [])
    

    @patch("atomium.pdb.mmcif_names_to_pdb_names")
    @patch("atomium.pdb.split_lines")
    def test_can_handle_one_name(self, mock_split, mock_names):
        mmcif = {"citation_author": [{"name": "John"}]}
        mock_names.return_value = "X,Y,Z"
        mock_split.return_value = ["NAMES1"]
        self.assertEqual(create_jrnl_auth_lines(mmcif), ["JRNL        AUTH   NAMES1"])
        mock_names.assert_called_with(["John"])
        mock_split.assert_called_with("X,Y,Z", 60)
    

    @patch("atomium.pdb.mmcif_names_to_pdb_names")
    @patch("atomium.pdb.split_lines")
    def test_can_handle_multi_name(self, mock_split, mock_names):
        mmcif = {"citation_author": [{"name": "John"}, {"name": "Flo"}]}
        mock_names.return_value = "X,Y,Z"
        mock_split.return_value = ["NAMES1", "NAMES2", "NAMES3"]
        self.assertEqual(create_jrnl_auth_lines(mmcif), [
            "JRNL        AUTH   NAMES1",
            "JRNL        AUTH 2 NAMES2",
            "JRNL        AUTH 3 NAMES3"
        ])
        mock_names.assert_called_with(["John", "Flo"])
        mock_split.assert_called_with("X,Y,Z", 60)



class JrnlTitlLinesTests(TestCase):

    def test_can_handle_no_category(self):
        self.assertEqual(create_jrnl_titl_lines({}), [])
    

    def test_can_handle_no_value(self):
        mmcif = {"citation": [{"title": "?"}]}
        self.assertEqual(create_jrnl_titl_lines(mmcif), [])
    

    @patch("atomium.pdb.split_lines")
    def test_can_save_single_line_title(self, mock_split):
        mmcif = {"citation": [{"title": "ttt"}]}
        mock_split.return_value = ["THE TITLE"]
        self.assertEqual(create_jrnl_titl_lines(mmcif), ["JRNL        TITL   THE TITLE"])
        mock_split.assert_called_with("TTT", 50)


    @patch("atomium.pdb.split_lines")
    def test_can_save_single_multi_line_title(self, mock_split):
        mmcif = {"citation": [{"title": "title"}]}
        mock_split.return_value = ["THE TITLE", "AND THE REST", "OF THE TITLE"]
        self.assertEqual(create_jrnl_titl_lines(mmcif), [
            "JRNL        TITL   THE TITLE",
            "JRNL        TITL 2 AND THE REST",
            "JRNL        TITL 3 OF THE TITLE",
        ])
        mock_split.assert_called_with("TITLE", 50)



class JrnlEditLinesTests(TestCase):

    def test_can_handle_no_table(self):
        self.assertEqual(create_jrnl_edit_lines({}), [])
    

    @patch("atomium.pdb.mmcif_names_to_pdb_names")
    @patch("atomium.pdb.split_lines")
    def test_can_handle_one_name(self, mock_split, mock_names):
        mmcif = {"citation_editor": [{"name": "John"}]}
        mock_names.return_value = "X,Y,Z"
        mock_split.return_value = ["NAMES1"]
        self.assertEqual(create_jrnl_edit_lines(mmcif), ["JRNL        EDIT   NAMES1"])
        mock_names.assert_called_with(["John"])
        mock_split.assert_called_with("X,Y,Z", 60)
    

    @patch("atomium.pdb.mmcif_names_to_pdb_names")
    @patch("atomium.pdb.split_lines")
    def test_can_handle_multi_name(self, mock_split, mock_names):
        mmcif = {"citation_editor": [{"name": "John"}, {"name": "Flo"}]}
        mock_names.return_value = "X,Y,Z"
        mock_split.return_value = ["NAMES1", "NAMES2", "NAMES3"]
        self.assertEqual(create_jrnl_edit_lines(mmcif), [
            "JRNL        EDIT   NAMES1",
            "JRNL        EDIT 2 NAMES2",
            "JRNL        EDIT 3 NAMES3"
        ])
        mock_names.assert_called_with(["John", "Flo"])
        mock_split.assert_called_with("X,Y,Z", 60)



class JrnlRefLinesTests(TestCase):

    def test_can_handle_no_category(self):
        self.assertEqual(create_jrnl_ref_lines({}), [])
    

    def test_can_handle_no_value(self):
        mmcif = {"citation": [
            {"journal_abbrev": "?", "journal_volume": "?", "page_first": "?", "year": "?"}
        ]}
        self.assertEqual(create_jrnl_ref_lines(mmcif), [])
    

    @patch("atomium.pdb.split_lines")
    def test_can_get_journal(self, mock_split):
        mmcif = {"citation": [
            {"journal_abbrev": "jname", "journal_volume": "?", "page_first": "?", "year": "?"}
        ]}
        mock_split.return_value = ["VOLUME NAME"]
        self.assertEqual(create_jrnl_ref_lines(mmcif), [
            "JRNL        REF    VOLUME NAME"
        ])
        mock_split.assert_called_with("JNAME", 28)
    

    @patch("atomium.pdb.split_lines")
    def test_can_get_journal_multiple_lines(self, mock_split):
        mmcif = {"citation": [
            {"journal_abbrev": "jname", "journal_volume": "?", "page_first": "?", "year": "?"}
        ]}
        mock_split.return_value = ["VOLUME NAME A", "VOLUME NAME B"]
        self.assertEqual(create_jrnl_ref_lines(mmcif), [
            "JRNL        REF    VOLUME NAME A",
            "JRNL        REF  2 VOLUME NAME B",
        ])
        mock_split.assert_called_with("JNAME", 28)
    

    @patch("atomium.pdb.split_lines")
    def test_can_get_volume(self, mock_split):
        mmcif = {"citation": [
            {"journal_abbrev": "?", "journal_volume": "277", "page_first": "?", "year": "?"}
        ]}
        mock_split.return_value = [""]
        self.assertEqual(create_jrnl_ref_lines(mmcif), [
            "JRNL        REF                                  V. 277"
        ])
        mock_split.assert_called_with("", 28)
    

    @patch("atomium.pdb.split_lines")
    def test_can_get_page(self, mock_split):
        mmcif = {"citation": [
            {"journal_abbrev": "?", "journal_volume": "?", "page_first": "100", "year": "?"}
        ]}
        mock_split.return_value = [""]
        self.assertEqual(create_jrnl_ref_lines(mmcif), [
            "JRNL        REF                                           100"
        ])
        mock_split.assert_called_with("", 28)
    

    @patch("atomium.pdb.split_lines")
    def test_can_get_year(self, mock_split):
        mmcif = {"citation": [
            {"journal_abbrev": "?", "journal_volume": "?", "page_first": "?", "year": "2002"}
        ]}
        mock_split.return_value = [""]
        self.assertEqual(create_jrnl_ref_lines(mmcif), [
            "JRNL        REF                                               2002"
        ])
        mock_split.assert_called_with("", 28)
    

    @patch("atomium.pdb.split_lines")
    def test_can_get_full_ref(self, mock_split):
        mmcif = {"citation": [
            {"journal_abbrev": "jname", "journal_volume": "277", "page_first": "28080", "year": "2002"}
        ]}
        mock_split.return_value = ["J.BIOL.CHEM.", "P.BIOL.CHEM."]
        self.assertEqual(create_jrnl_ref_lines(mmcif), [
            "JRNL        REF    J.BIOL.CHEM.                  V. 277 28080 2002",
            "JRNL        REF  2 P.BIOL.CHEM.",
        ])
        mock_split.assert_called_with("JNAME", 28)




class JrnlPublLinesTests(TestCase):

    def test_can_handle_no_category(self):
        self.assertEqual(create_jrnl_publ_lines({}), [])
    

    def test_can_handle_no_value(self):
        mmcif = {"citation": [{"book_publisher": "?"}]}
        self.assertEqual(create_jrnl_publ_lines(mmcif), [])
    

    @patch("atomium.pdb.split_lines")
    def test_can_save_single_line_publ(self, mock_split):
        mmcif = {"citation": [{"book_publisher": "ppp"}]}
        mock_split.return_value = ["THE PUBLISHER"]
        self.assertEqual(create_jrnl_publ_lines(mmcif), ["JRNL        PUBL   THE PUBLISHER"])
        mock_split.assert_called_with("PPP", 60)


    @patch("atomium.pdb.split_lines")
    def test_can_save_single_multi_line_publ(self, mock_split):
        mmcif = {"citation": [{"book_publisher": "publi"}]}
        mock_split.return_value = ["THE PUBLISHER", "AND THE REST", "OF THE PUBLISHER"]
        self.assertEqual(create_jrnl_publ_lines(mmcif), [
            "JRNL        PUBL   THE PUBLISHER",
            "JRNL        PUBL 2 AND THE REST",
            "JRNL        PUBL 3 OF THE PUBLISHER",
        ])
        mock_split.assert_called_with("PUBLI", 60)



class JrnlRefnLinesTests(TestCase):

    def test_can_handle_no_category(self):
        self.assertEqual(create_jrnl_refn_lines({}), [])
    

    def test_can_handle_no_value(self):
        mmcif = {"citation": [{"journal_id_ISSN": "?"}]}
        self.assertEqual(create_jrnl_refn_lines(mmcif), [])
    

    def test_can_produce_line(self):
        mmcif = {"citation": [{"journal_id_ISSN": "10026008"}]}
        self.assertEqual(create_jrnl_refn_lines(mmcif), [
            "JRNL        REFN                   ISSN 10026008"
        ])



class JrnlPmidLinesTests(TestCase):

    def test_can_handle_no_category(self):
        self.assertEqual(create_jrnl_pmid_lines({}), [])
    

    def test_can_handle_no_value(self):
        mmcif = {"citation": [{"pdbx_database_id_PubMed": "?"}]}
        self.assertEqual(create_jrnl_pmid_lines(mmcif), [])
    

    def test_can_produce_line(self):
        mmcif = {"citation": [{"pdbx_database_id_PubMed": "70008"}]}
        self.assertEqual(create_jrnl_pmid_lines(mmcif), [
            "JRNL        PMID   70008"
        ])



class JrnlDoiLinesTests(TestCase):

    def test_can_handle_no_category(self):
        self.assertEqual(create_jrnl_doi_lines({}), [])
    

    def test_can_handle_no_value(self):
        mmcif = {"citation": [{"pdbx_database_id_DOI": "?"}]}
        self.assertEqual(create_jrnl_doi_lines(mmcif), [])
    

    def test_can_produce_line(self):
        mmcif = {"citation": [{"pdbx_database_id_DOI": "ddd.100"}]}
        self.assertEqual(create_jrnl_doi_lines(mmcif), [
            "JRNL        DOI    ddd.100"
        ])



class RemarkLinesTests(TestCase):

    @patch("atomium.pdb.create_remark_2_lines")
    @patch("atomium.pdb.create_remark_3_lines")
    @patch("atomium.pdb.create_remark_350_lines")
    @patch("atomium.pdb.create_remark_465_lines")
    @patch("atomium.pdb.create_remark_800_lines")
    def test_can_produce_remark_lines(self, *mocks):
        mmcif = {"mmcif": 1}
        for i, mock in enumerate(mocks[::-1]):
            mock.return_value = [i]
        lines = create_remark_lines(mmcif)
        self.assertEqual(lines, [0, 1, 2, 3, 4])
        for mock in mocks:
            mock.assert_called_with(mmcif)



class Remark2LinesTests(TestCase):

    def test_can_handle_no_table(self):
        self.assertEqual(create_remark_2_lines({}), [])
    

    def test_can_handle_no_value(self):
        mmcif = {"reflns": [{"d_resolution_high": "?"}]}
        self.assertEqual(create_remark_2_lines(mmcif), [])
    

    def test_can_get_remark_2(self):
        mmcif = {"reflns": [{"d_resolution_high": "2.3"}]}
        self.assertEqual(create_remark_2_lines(mmcif), [
             "REMARK   2",  "REMARK   2 RESOLUTION.    2.30 ANGSTROMS."
        ])



class Remark3LinesTests(TestCase):

    def test_can_handle_no_table(self):
        self.assertEqual(create_remark_3_lines({}), [])
    

    @patch("atomium.pdb.create_remark_3_refinement_target_lines")
    @patch("atomium.pdb.create_remark_3_data_used_lines")
    @patch("atomium.pdb.create_remark_3_data_fit_lines")
    @patch("atomium.pdb.create_remark_3_bvalue_lines")
    @patch("atomium.pdb.create_remark_3_thermal_lines")
    @patch("atomium.pdb.create_remark_3_solvent_lines")
    def test_can_parse_remark_3(self, *mocks):
        mmcif = {"refine": 1}
        for i, mock in enumerate(mocks[::-1]):
            mock.return_value = [i]
        lines = create_remark_3_lines(mmcif)
        self.assertEqual(lines, [0, 1, 2, 3, 4, 5])
        for mock in mocks:
            mock.assert_called_with(mmcif)



class Remark3RefinementTargetLinesTests(TestCase):

    def test_can_handle_no_property(self):
        mmcif = {"refine": [{}]}
        self.assertEqual(create_remark_3_refinement_target_lines(mmcif), [
            "REMARK   3", "REMARK   3  REFINEMENT TARGET : NULL"
        ])
    

    def test_can_handle_no_value(self):
        mmcif = {"refine": [{"pdbx_stereochemistry_target_values": "?"}]}
        self.assertEqual(create_remark_3_refinement_target_lines(mmcif), [
            "REMARK   3", "REMARK   3  REFINEMENT TARGET : NULL"
        ])
    

    def test_can_get_values(self):
        mmcif = {"refine": [{"pdbx_stereochemistry_target_values": "EH"}]}
        self.assertEqual(create_remark_3_refinement_target_lines(mmcif), [
            "REMARK   3", "REMARK   3  REFINEMENT TARGET : EH"
        ])



class Remark3DataUsedLinesTests(TestCase):

    def test_can_handle_no_property(self):
        mmcif = {"refine": [{}]}
        self.assertEqual(create_remark_3_data_used_lines(mmcif), [
            "REMARK   3",
            "REMARK   3  DATA USED IN REFINEMENT.",
            "REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : NULL",
            "REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : NULL",
            "REMARK   3   DATA CUTOFF            (SIGMA(F)) : NULL",
            "REMARK   3   DATA CUTOFF HIGH         (ABS(F)) : NULL",
            "REMARK   3   DATA CUTOFF LOW          (ABS(F)) : NULL",
            "REMARK   3   COMPLETENESS (WORKING+TEST)   (%) : NULL",
            "REMARK   3   NUMBER OF REFLECTIONS             : NULL"
        ])
    

    def test_can_handle_no_value(self):
        mmcif = {"refine": [{
            "ls_d_res_high": "?", "ls_d_res_low": "?", "pdbx_data_cutoff_high_absF": "?", 
            "pdbx_data_cutoff_low_absF": "?", "ls_percent_reflns_obs": "?", "ls_number_reflns_obs": "?", 
        }]}
        self.assertEqual(create_remark_3_data_used_lines(mmcif), [
            "REMARK   3",
            "REMARK   3  DATA USED IN REFINEMENT.",
            "REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : NULL",
            "REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : NULL",
            "REMARK   3   DATA CUTOFF            (SIGMA(F)) : NULL",
            "REMARK   3   DATA CUTOFF HIGH         (ABS(F)) : NULL",
            "REMARK   3   DATA CUTOFF LOW          (ABS(F)) : NULL",
            "REMARK   3   COMPLETENESS (WORKING+TEST)   (%) : NULL",
            "REMARK   3   NUMBER OF REFLECTIONS             : NULL"
        ])
    

    def test_can_get_values(self):
        mmcif = {"refine": [{
            "ls_d_res_high": "1.2", "ls_d_res_low": "1.1", "pdbx_data_cutoff_high_absF": "22.1", 
            "pdbx_data_cutoff_low_absF": "100.3", "ls_percent_reflns_obs": "17", "ls_number_reflns_obs": "4000", 
        }]}
        self.assertEqual(create_remark_3_data_used_lines(mmcif), [
            "REMARK   3",
            "REMARK   3  DATA USED IN REFINEMENT.",
            "REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 1.2",
            "REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 1.1",
            "REMARK   3   DATA CUTOFF            (SIGMA(F)) : 100.3",
            "REMARK   3   DATA CUTOFF HIGH         (ABS(F)) : 22.1",
            "REMARK   3   DATA CUTOFF LOW          (ABS(F)) : 100.3",
            "REMARK   3   COMPLETENESS (WORKING+TEST)   (%) : 17",
            "REMARK   3   NUMBER OF REFLECTIONS             : 4000"
        ])



class Remark3DataFitLinesTests(TestCase):

    def test_can_handle_no_property(self):
        mmcif = {"refine": [{}]}
        self.assertEqual(create_remark_3_data_fit_lines(mmcif), [
            "REMARK   3",
            "REMARK   3  FIT TO DATA USED IN REFINEMENT.", 
            "REMARK   3   CROSS-VALIDATION METHOD          : NULL",
            "REMARK   3   FREE R VALUE TEST SET SELECTION  : NULL",
            "REMARK   3   R VALUE            (WORKING SET) : NULL",
            "REMARK   3   FREE R VALUE                     : NULL",
            "REMARK   3   FREE R VALUE TEST SET SIZE   (%) : NULL",
            "REMARK   3   FREE R VALUE TEST SET COUNT      : NULL",
            "REMARK   3   ESTIMATED ERROR OF FREE R VALUE  : NULL",
        ])
    

    def test_can_handle_no_values(self):
        mmcif = {"refine": [{
            "pdbx_ls_cross_valid_method": "?", "pdbx_R_Free_selection_details": "?",
            "ls_R_factor_R_work": "?", "ls_R_factor_R_free": "?", "ls_percent_reflns_R_free": "?",
            "ls_number_reflns_R_free": "?", "ls_R_factor_R_free_error": "?"
        }]}
        self.assertEqual(create_remark_3_data_fit_lines(mmcif), [
            "REMARK   3",
            "REMARK   3  FIT TO DATA USED IN REFINEMENT.", 
            "REMARK   3   CROSS-VALIDATION METHOD          : NULL",
            "REMARK   3   FREE R VALUE TEST SET SELECTION  : NULL",
            "REMARK   3   R VALUE            (WORKING SET) : NULL",
            "REMARK   3   FREE R VALUE                     : NULL",
            "REMARK   3   FREE R VALUE TEST SET SIZE   (%) : NULL",
            "REMARK   3   FREE R VALUE TEST SET COUNT      : NULL",
            "REMARK   3   ESTIMATED ERROR OF FREE R VALUE  : NULL",
        ])
    

    def test_can_get_values(self):
        mmcif = {"refine": [{
            "pdbx_ls_cross_valid_method": "THROUGH", "pdbx_R_Free_selection_details": "RAND",
            "ls_R_factor_R_work": "0.3", "ls_R_factor_R_free": "0.4", "ls_percent_reflns_R_free": "50",
            "ls_number_reflns_R_free": "4.5", "ls_R_factor_R_free_error": "3.2"
        }]}
        self.assertEqual(create_remark_3_data_fit_lines(mmcif), [
            "REMARK   3",
            "REMARK   3  FIT TO DATA USED IN REFINEMENT.", 
            "REMARK   3   CROSS-VALIDATION METHOD          : THROUGH",
            "REMARK   3   FREE R VALUE TEST SET SELECTION  : RAND",
            "REMARK   3   R VALUE            (WORKING SET) : 0.3",
            "REMARK   3   FREE R VALUE                     : 0.4",
            "REMARK   3   FREE R VALUE TEST SET SIZE   (%) : 50",
            "REMARK   3   FREE R VALUE TEST SET COUNT      : 4.5",
            "REMARK   3   ESTIMATED ERROR OF FREE R VALUE  : 3.2",
        ])



class Remark3BvalueLinesTests(TestCase):

    def test_can_handle_no_property(self):
        mmcif = {"refine": [{}]}
        self.assertEqual(create_remark_3_bvalue_lines(mmcif), [
            "REMARK   3",
            "REMARK   3  B VALUES.",
            "REMARK   3   FROM WILSON PLOT           (A**2) : NULL",
            "REMARK   3   MEAN B VALUE      (OVERALL, A**2) : NULL",
            "REMARK   3   OVERALL ANISOTROPIC B VALUE.",
            "REMARK   3    B11 (A**2) : NULL",
            "REMARK   3    B22 (A**2) : NULL",
            "REMARK   3    B33 (A**2) : NULL",
            "REMARK   3    B12 (A**2) : NULL",
            "REMARK   3    B13 (A**2) : NULL",
            "REMARK   3    B23 (A**2) : NULL"
        ])
    

    def test_can_handle_no_values(self):
        mmcif = {"reflns": [{"B_iso_Wilson_estimate": "?"}], "refine": [{
            "B_iso_mean": "?", "aniso_B[1][1]": "?", "aniso_B[2][2]": "?", "aniso_B[3][3]": "?",
            "aniso_B[1][2]": "?", "aniso_B[1][3]": "?", "aniso_B[2][3]": "?"
        }]}
        self.assertEqual(create_remark_3_bvalue_lines(mmcif), [
            "REMARK   3",
            "REMARK   3  B VALUES.",
            "REMARK   3   FROM WILSON PLOT           (A**2) : NULL",
            "REMARK   3   MEAN B VALUE      (OVERALL, A**2) : NULL",
            "REMARK   3   OVERALL ANISOTROPIC B VALUE.",
            "REMARK   3    B11 (A**2) : NULL",
            "REMARK   3    B22 (A**2) : NULL",
            "REMARK   3    B33 (A**2) : NULL",
            "REMARK   3    B12 (A**2) : NULL",
            "REMARK   3    B13 (A**2) : NULL",
            "REMARK   3    B23 (A**2) : NULL"
        ])
    

    def test_can_handle_get_values(self):
        mmcif = {"reflns": [{"B_iso_Wilson_estimate": "12.2"}], "refine": [{
            "B_iso_mean": "24.2", "aniso_B[1][1]": "1", "aniso_B[2][2]": "2", "aniso_B[3][3]": "3",
            "aniso_B[1][2]": "4", "aniso_B[1][3]": "5", "aniso_B[2][3]": "6"
        }]}
        self.assertEqual(create_remark_3_bvalue_lines(mmcif), [
            "REMARK   3",
            "REMARK   3  B VALUES.",
            "REMARK   3   FROM WILSON PLOT           (A**2) : 12.2",
            "REMARK   3   MEAN B VALUE      (OVERALL, A**2) : 24.2",
            "REMARK   3   OVERALL ANISOTROPIC B VALUE.",
            "REMARK   3    B11 (A**2) : 1",
            "REMARK   3    B22 (A**2) : 2",
            "REMARK   3    B33 (A**2) : 3",
            "REMARK   3    B12 (A**2) : 4",
            "REMARK   3    B13 (A**2) : 5",
            "REMARK   3    B23 (A**2) : 6"
        ])



class Remark3ThermalLinesTests(TestCase):

    def test_can_handle_no_property(self):
        mmcif = {"refine": [{}]}
        self.assertEqual(create_remark_3_thermal_lines(mmcif), [
            "REMARK   3",
            "REMARK   3  ISOTROPIC THERMAL MODEL : NULL",
        ])
    

    def test_can_handle_no_value(self):
        mmcif = {"refine": [{"pdbx_isotropic_thermal_model": "?"}]}
        self.assertEqual(create_remark_3_thermal_lines(mmcif), [
            "REMARK   3",
            "REMARK   3  ISOTROPIC THERMAL MODEL : NULL",
        ])
    

    def test_can_get_value(self):
        mmcif = {"refine": [{"pdbx_isotropic_thermal_model": "REST"}]}
        self.assertEqual(create_remark_3_thermal_lines(mmcif), [
            "REMARK   3",
            "REMARK   3  ISOTROPIC THERMAL MODEL : REST",
        ])



class Remark3SolventLinesTests(TestCase):

    def test_can_handle_no_property(self):
        mmcif = {"refine": [{}]}
        self.assertEqual(create_remark_3_solvent_lines(mmcif), [
            "REMARK   3",
            "REMARK   3  BULK SOLVENT MODELING.",
            "REMARK   3   METHOD USED : NULL",
            "REMARK   3   KSOL        : NULL",
            "REMARK   3   BSOL        : NULL"
        ])
    

    def test_can_handle_no_value(self):
        mmcif = {"refine": [{
            "solvent_model_details": "?",
            "solvent_model_param_ksol": "?", "solvent_model_param_bsol": "?"
        }]}
        self.assertEqual(create_remark_3_solvent_lines(mmcif), [
            "REMARK   3",
            "REMARK   3  BULK SOLVENT MODELING.",
            "REMARK   3   METHOD USED : NULL",
            "REMARK   3   KSOL        : NULL",
            "REMARK   3   BSOL        : NULL"
        ])
    

    def test_can_get_values(self):
        mmcif = {"refine": [{
            "solvent_model_details": "FLAT",
            "solvent_model_param_ksol": "0.5", "solvent_model_param_bsol": "48"
        }]}
        self.assertEqual(create_remark_3_solvent_lines(mmcif), [
            "REMARK   3",
            "REMARK   3  BULK SOLVENT MODELING.",
            "REMARK   3   METHOD USED : FLAT",
            "REMARK   3   KSOL        : 0.5",
            "REMARK   3   BSOL        : 48"
        ])



class Remark465LinesTests(TestCase):

    def test_can_handle_no_table(self):
        self.assertEqual(create_remark_465_lines({}), [])
    

    def test_can_produce_lines(self):
        mmcif = {"pdbx_unobs_or_zero_occ_residues": [
            {"auth_comp_id": "MET", "auth_asym_id": "A", "auth_seq_id": "1", "PDB_ins_code": "?"},
            {"auth_comp_id": "VAL", "auth_asym_id": "A", "auth_seq_id": "10", "PDB_ins_code": "?"},
            {"auth_comp_id": "PRO", "auth_asym_id": "A", "auth_seq_id": "10", "PDB_ins_code": "A"},
            {"auth_comp_id": "TYR", "auth_asym_id": "B", "auth_seq_id": "999", "PDB_ins_code": "?"},
            {"auth_comp_id": "CYS", "auth_asym_id": "B", "auth_seq_id": "4502", "PDB_ins_code": "?"},
        ]}
        self.assertEqual(create_remark_465_lines(mmcif), [
            "REMARK 465",
            "REMARK 465 MISSING RESIDUES",
            "REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE",
            "REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN",
            "REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)",
            "REMARK 465     MET A     1 ",
            "REMARK 465     VAL A    10 ",
            "REMARK 465     PRO A    10A",
            "REMARK 465     TYR B   999 ",
            "REMARK 465     CYS B  4502 "
        ])



class Remark800LinesTests(TestCase):

    def test_can_handle_no_table(self):
        self.assertEqual(create_remark_800_lines({}), [])
    

    def test_can_produce_lines(self):
        mmcif = {"struct_site": [
            {"id": "s1", "pdbx_evidence_code": "Software", "details": "site for res1"},
            {"id": "s2", "pdbx_evidence_code": "Thought", "details": "site for res2"},
        ]}
        self.assertEqual(create_remark_800_lines(mmcif), [
            "REMARK 800",
            "REMARK 800 SITE",
            "REMARK 800 SITE_IDENTIFIER: S1",
            "REMARK 800 EVIDENCE_CODE: SOFTWARE",
            "REMARK 800 SITE_DESCRIPTION: SITE FOR RES1",
            "REMARK 800 SITE_IDENTIFIER: S2",
            "REMARK 800 EVIDENCE_CODE: THOUGHT",
            "REMARK 800 SITE_DESCRIPTION: SITE FOR RES2"
        ])
    

    def test_can_handle_no_values(self):
        mmcif = {"struct_site": [
            {"id": "s1", "pdbx_evidence_code": "?", "details": "?"},
            {"id": "?", "pdbx_evidence_code": "Thought", "details": "site for res2"},
        ]}
        self.assertEqual(create_remark_800_lines(mmcif), [
            "REMARK 800",
            "REMARK 800 SITE",
            "REMARK 800 SITE_IDENTIFIER: S1"
        ])



class Remark350LinesTests(TestCase):

    def test_can_handle_no_category(self):
        self.assertEqual(create_remark_350_lines({}), [])


    @patch("atomium.pdb.create_single_assembly_lines")
    def test_can_parse_assemblies(self, mock_lines):
        mmcif = {"pdbx_struct_assembly": [1, 2, 3], "atom_site": [
            {"label_asym_id": "A", "auth_asym_id": "A"},
            {"label_asym_id": "A", "auth_asym_id": "A"},
            {"label_asym_id": "B", "auth_asym_id": "B"},
            {"label_asym_id": "B", "auth_asym_id": "B"},
            {"label_asym_id": "C", "auth_asym_id": "A"},
            {"label_asym_id": "C", "auth_asym_id": "A"},
            {"label_asym_id": "D", "auth_asym_id": "B"},
            {"label_asym_id": "D", "auth_asym_id": "B"},
        ]}
        mock_lines.side_effect = lambda a, m, l: [a] * 2
        self.assertEqual(create_remark_350_lines(mmcif), [
            "REMARK 350",
            "REMARK 350 COORDINATES FOR A COMPLETE MULTIMER REPRESENTING THE KNOWN",
            "REMARK 350 BIOLOGICALLY SIGNIFICANT OLIGOMERIZATION STATE OF THE",
            "REMARK 350 MOLECULE CAN BE GENERATED BY APPLYING BIOMT TRANSFORMATIONS",
            "REMARK 350 GIVEN BELOW.  BOTH NON-CRYSTALLOGRAPHIC AND",
            "REMARK 350 CRYSTALLOGRAPHIC OPERATIONS ARE GIVEN.", 1, 1, 2, 2, 3, 3
        ])
        mock_lines.assert_any_call(1, mmcif, {"A": "A", "B": "B", "C": "A", "D": "B"})
        mock_lines.assert_any_call(2, mmcif, {"A": "A", "B": "B", "C": "A", "D": "B"})
        mock_lines.assert_any_call(3, mmcif, {"A": "A", "B": "B", "C": "A", "D": "B"})



class SingleAssemblyLinesTests(TestCase):

    def test_can_get_minimal_assembly(self):
        assembly = {"id": "1"}
        mmcif = {}
        asym_lookup = {1: 2}
        self.assertEqual(create_single_assembly_lines(assembly, mmcif, asym_lookup), [
            "REMARK 350",
            "REMARK 350 BIOMOLECULE: 1"
        ])
    

    def test_can_get_assembly_with_properties(self):
        assembly = {"id": "1", "method_details": "PISA"}
        mmcif = {"pdbx_struct_assembly_prop": [
            {"biol_id": "1", "type": "MORE", "value": "2.3"},
            {"biol_id": "1", "type": "SSA (A^2)", "value": "123"},
            {"biol_id": "1", "type": "ABSA (A^2)", "value": "5674"},
            {"biol_id": "1", "type": "XYZ", "value": "999"},
            {"biol_id": "2", "type": "MORE", "value": "3.4"},
            {"biol_id": "2", "type": "SSA (A^2)", "value": "777"},
        ]}
        asym_lookup = {1: 2}
        self.assertEqual(create_single_assembly_lines(assembly, mmcif, asym_lookup), [
            "REMARK 350",
            "REMARK 350 BIOMOLECULE: 1",
            "REMARK 350 SOFTWARE USED: PISA",
            "REMARK 350 TOTAL BURIED SURFACE AREA: 5674 ANGSTROM**2",
            "REMARK 350 SURFACE AREA OF THE COMPLEX: 123 ANGSTROM**2",
            "REMARK 350 CHANGE IN SOLVENT FREE ENERGY: 2.3 KCAL/MOL"
        ])
    

    @patch("atomium.pdb.create_assembly_gen_lines")
    def test_can_produce_gen_lines(self, mock_lines):
        assembly = {"id": "1", "method_details": "?"}
        mmcif = {"pdbx_struct_assembly_gen": [
            {"assembly_id": "1", "gen": "1"},
            {"assembly_id": "1", "gen": "2"},
            {"assembly_id": "2", "gen": "3"},
            {"assembly_id": "2", "gen": "4"},
        ]}
        mock_lines.side_effect = [["LINE 1", "LINE 2"], ["LINE 3", "LINE 4"]]
        asym_lookup = {1: 2}
        self.assertEqual(create_single_assembly_lines(assembly, mmcif, asym_lookup), [
            "REMARK 350",
            "REMARK 350 BIOMOLECULE: 1",
            "LINE 1", "LINE 2", "LINE 3", "LINE 4",
        ])
    

    @patch("atomium.pdb.create_assembly_gen_lines")
    def test_can_produce_full_assembly_lines(self, mock_lines):
        assembly = {"id": "1", "method_details": "PISA"}
        mmcif = {"pdbx_struct_assembly_gen": [
            {"assembly_id": "1", "gen": "1"},
            {"assembly_id": "1", "gen": "2"},
            {"assembly_id": "2", "gen": "3"},
            {"assembly_id": "2", "gen": "4"},
        ], "pdbx_struct_assembly_prop": [
            {"biol_id": "1", "type": "MORE", "value": "2.3"},
            {"biol_id": "1", "type": "SSA (A^2)", "value": "123"},
            {"biol_id": "1", "type": "ABSA (A^2)", "value": "5674"},
            {"biol_id": "1", "type": "XYZ", "value": "999"},
            {"biol_id": "2", "type": "MORE", "value": "3.4"},
            {"biol_id": "2", "type": "SSA (A^2)", "value": "777"},
        ]}
        mock_lines.side_effect = [["LINE 1", "LINE 2"], ["LINE 3", "LINE 4"]]
        asym_lookup = {1: 2}
        self.assertEqual(create_single_assembly_lines(assembly, mmcif, asym_lookup), [
            "REMARK 350",
            "REMARK 350 BIOMOLECULE: 1",
            "REMARK 350 SOFTWARE USED: PISA",
            "REMARK 350 TOTAL BURIED SURFACE AREA: 5674 ANGSTROM**2",
            "REMARK 350 SURFACE AREA OF THE COMPLEX: 123 ANGSTROM**2",
            "REMARK 350 CHANGE IN SOLVENT FREE ENERGY: 2.3 KCAL/MOL",
            "LINE 1", "LINE 2", "LINE 3", "LINE 4",
        ])



class AssemblyGenLinesTests(TestCase):

    @patch("atomium.pdb.split_lines")
    @patch("atomium.pdb.get_operations")
    def test_can_get_single_chain_line(self, mock_ops, mock_lines):
        gen = {"asym_id_list": "A,B,C,D,E", "oper_expression": "EXP"}
        mmcif = {1: 2}
        asym_lookup = {"A": "A", "B": "B", "C": "C", "D": "A", "E": "B"}
        mock_lines.return_value = ["A, B, C"]
        mock_ops.return_value = [
            [[1, 0, 0, 1], [0, -1, 0, 2], [0, 0, -1, 1], [1, 0, 0, 1]],
            [[0.2234, -1009.2, 5, -2], [4.2334, -8.993, 1.999887766554433, 4.04], [0, 0, -1, 1], [1, 0, 0, 1]],
        ]
        self.assertEqual(create_assembly_gen_lines(gen, mmcif, asym_lookup), [
            "REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, B, C,",
            "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        1.000000",
            "REMARK 350   BIOMT2   1  0.000000 -1.000000  0.000000        2.000000",
            "REMARK 350   BIOMT3   1  0.000000  0.000000 -1.000000        1.000000",
            "REMARK 350   BIOMT1   2  0.223400 -1009.200  5.000000       -2.000000",
            "REMARK 350   BIOMT2   2  4.233400 -8.993000  1.999887        4.040000",
            "REMARK 350   BIOMT3   2  0.000000  0.000000 -1.000000        1.000000"
        ])
        mock_lines.assert_called_with("A, B, C", 27)
        mock_ops.assert_called_with({1: 2}, gen)
    

    @patch("atomium.pdb.split_lines")
    @patch("atomium.pdb.get_operations")
    def test_can_get_single_chain_line(self, mock_ops, mock_lines):
        gen = {"asym_id_list": "A,B,C,D,E", "oper_expression": "EXP"}
        mmcif = {1: 2}
        asym_lookup = {"A": "A", "B": "B", "C": "C", "D": "A", "E": "B"}
        mock_lines.return_value = ["A, B, C", "D, E"]
        mock_ops.return_value = [
            [[1, 0, 0, 1], [0, -1, 0, 2], [0, 0, -1, 1], [1, 0, 0, 1]],
            [[0.2234, -1009.2, 5, -2], [4.2334, -8.993, 1.999887766554433, 4.04], [0, 0, -1, 1], [1, 0, 0, 1]],
        ]
        self.assertEqual(create_assembly_gen_lines(gen, mmcif, asym_lookup), [
            "REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, B, C,",
            "REMARK 350                    AND CHAINS: D, E,",
            "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        1.000000",
            "REMARK 350   BIOMT2   1  0.000000 -1.000000  0.000000        2.000000",
            "REMARK 350   BIOMT3   1  0.000000  0.000000 -1.000000        1.000000",
            "REMARK 350   BIOMT1   2  0.223400 -1009.200  5.000000       -2.000000",
            "REMARK 350   BIOMT2   2  4.233400 -8.993000  1.999887        4.040000",
            "REMARK 350   BIOMT3   2  0.000000  0.000000 -1.000000        1.000000"
        ])
        mock_lines.assert_called_with("A, B, C", 27)
        mock_ops.assert_called_with({1: 2}, gen)



class DbrefLinesTests(TestCase):

    def test_can_handle_no_table(self):
        self.assertEqual(create_dbref_lines({}), [])
    

    def test_can_save_short_sbref(self):
        mmcif = {"struct_ref": [
            {"id": "1", "db_name": "UNP", "db_code": "PYRF_METTH", "pdbx_db_accession": "O26232"},
            {"id": "2", "db_name": "REF", "db_code": "PYRF_SEC", "pdbx_db_accession": "P26232"},
        ], "struct_ref_seq": [{
            "ref_id": "1", "pdbx_strand_id": "A", "pdbx_auth_seq_align_beg": "10",
            "pdbx_seq_align_beg_ins_code": "F", "pdbx_seq_align_end_ins_code": "?",
            "pdbx_auth_seq_align_end": "20", "db_align_beg": "11", "db_align_end": "21",
            "pdbx_db_align_beg_ins_code": "?", "pdbx_db_align_end_ins_code": "?"
        }, {
            "ref_id": "1", "pdbx_strand_id": "B", "pdbx_auth_seq_align_beg": "100",
            "pdbx_seq_align_beg_ins_code": "?", "pdbx_seq_align_end_ins_code": "G",
            "pdbx_auth_seq_align_end": "200", "db_align_beg": "101", "db_align_end": "201",
            "pdbx_db_align_beg_ins_code": "?", "pdbx_db_align_end_ins_code": "?"
        }, {
            "ref_id": "2", "pdbx_strand_id": "C", "pdbx_auth_seq_align_beg": "1000",
            "pdbx_seq_align_beg_ins_code": "?", "pdbx_seq_align_end_ins_code": "?",
            "pdbx_auth_seq_align_end": "2000", "db_align_beg": "1001", "db_align_end": "2001",
            "pdbx_db_align_beg_ins_code": "P", "pdbx_db_align_end_ins_code": "?"
        }, {
            "ref_id": "2", "pdbx_strand_id": "D", "pdbx_auth_seq_align_beg": "6000",
            "pdbx_seq_align_beg_ins_code": "?", "pdbx_seq_align_end_ins_code": "?",
            "pdbx_auth_seq_align_end": "7000", "db_align_beg": "6001", "db_align_end": "7001",
            "pdbx_db_align_beg_ins_code": "?", "pdbx_db_align_end_ins_code": "Z"
        }]}
        self.assertEqual(create_dbref_lines(mmcif), [
            "DBREF       A   10F   20  UNP    O26232   PYRF_METTH      11     21 ",
            "DBREF       B  100   200G UNP    O26232   PYRF_METTH     101    201 ",
            "DBREF       C 1000  2000  REF    P26232   PYRF_SEC      1001P  2001 ",
            "DBREF       D 6000  7000  REF    P26232   PYRF_SEC      6001   7001Z"
        ])
    

    @patch("atomium.pdb.create_dbrefn_lines")
    def test_can_save_long_dbref(self, mock_create):
        mock_create.side_effect = [["L1", "L2"], ["L3", "L4"]]
        mmcif = {"struct_ref": [
            {"id": "1", "db_name": "UNP", "db_code": "A0A1B8Y853_XENTR", "pdbx_db_accession": "A0A1B8Y853"},
            {"id": "2", "db_name": "REF", "db_code": "PYRF_SEC", "pdbx_db_accession": "P26232"},
        ], "struct_ref_seq": [{
            "ref_id": "1", "pdbx_strand_id": "A", "pdbx_auth_seq_align_beg": "10",
            "pdbx_seq_align_beg_ins_code": "Y", "pdbx_seq_align_end_ins_code": "L",
            "pdbx_auth_seq_align_end": "20", "db_align_beg": "11", "db_align_end": "21",
            "pdbx_db_align_beg_ins_code": "O", "pdbx_db_align_end_ins_code": "X"
        }, {
            "ref_id": "1", "pdbx_strand_id": "B", "pdbx_auth_seq_align_beg": "100",
            "pdbx_seq_align_beg_ins_code": "?", "pdbx_seq_align_end_ins_code": "?",
            "pdbx_auth_seq_align_end": "200", "db_align_beg": "101", "db_align_end": "201",
            "pdbx_db_align_beg_ins_code": "?", "pdbx_db_align_end_ins_code": "?"
        }, {
            "ref_id": "2", "pdbx_strand_id": "C", "pdbx_auth_seq_align_beg": "1000",
            "pdbx_seq_align_beg_ins_code": "?", "pdbx_seq_align_end_ins_code": "?",
            "pdbx_auth_seq_align_end": "2000", "db_align_beg": "1001", "db_align_end": "2001",
            "pdbx_db_align_beg_ins_code": "P", "pdbx_db_align_end_ins_code": "?"
        }, {
            "ref_id": "2", "pdbx_strand_id": "D", "pdbx_auth_seq_align_beg": "6000",
            "pdbx_seq_align_beg_ins_code": "?", "pdbx_seq_align_end_ins_code": "?",
            "pdbx_auth_seq_align_end": "7000", "db_align_beg": "6001", "db_align_end": "7001",
            "pdbx_db_align_beg_ins_code": "?", "pdbx_db_align_end_ins_code": "Z"
        }]}
        self.assertEqual(create_dbref_lines(mmcif), [
            "L1", "L2", "L3", "L4",
            "DBREF       C 1000  2000  REF    P26232   PYRF_SEC      1001P  2001 ",
            "DBREF       D 6000  7000  REF    P26232   PYRF_SEC      6001   7001Z"
        ])
        mock_create.assert_any_call("", mmcif["struct_ref"][0], mmcif["struct_ref_seq"][0])
        mock_create.assert_any_call("", mmcif["struct_ref"][0], mmcif["struct_ref_seq"][1])
    
    

class DbrefnLinesTests(TestCase):

    def test_full_values(self):
        seq = {
            "ref_id": "1", "pdbx_strand_id": "A", "pdbx_auth_seq_align_beg": "10",
            "pdbx_seq_align_beg_ins_code": "Y", "pdbx_seq_align_end_ins_code": "L",
            "pdbx_auth_seq_align_end": "20", "db_align_beg": "11", "db_align_end": "21",
            "pdbx_db_align_beg_ins_code": "O", "pdbx_db_align_end_ins_code": "X"
        }
        ref = {"id": "1", "db_name": "UNP", "db_code": "A0A1B8Y853_XENTR", "pdbx_db_accession": "A0A1B8Y853"}
        self.assertEqual(create_dbrefn_lines("1xxx", ref, seq), [
            "DBREF1 1xxx A   10Y   20L UNP                  A0A1B8Y853_XENTR",
            "DBREF2 1xxx A     A0A1B8Y853                         11O         21X",
        ])
    

    def test_blank_values(self):
        seq = {
            "ref_id": "1", "pdbx_strand_id": "B", "pdbx_auth_seq_align_beg": "100",
            "pdbx_seq_align_beg_ins_code": "?", "pdbx_seq_align_end_ins_code": "?",
            "pdbx_auth_seq_align_end": "200", "db_align_beg": "101", "db_align_end": "201",
            "pdbx_db_align_beg_ins_code": "?", "pdbx_db_align_end_ins_code": "?"
        }
        ref = {"id": "1", "db_name": "UNP", "db_code": "A0A1B8Y853_XENTR", "pdbx_db_accession": "A0A1B8Y853"}
        self.assertEqual(create_dbrefn_lines("1xxx", ref, seq), [
            "DBREF1 1xxx B  100   200  UNP                  A0A1B8Y853_XENTR",
            "DBREF2 1xxx B     A0A1B8Y853                        101         201 ",
        ])



class SeqadvLinesTests(TestCase):

    def test_can_handle_no_table(self):
        self.assertEqual(create_seqadv_lines({}), [])
    

    def test_can_create_lines(self):
        mmcif = {
            "struct_ref_seq_dif": [{
                "pdbx_pdb_id_code": "1LOL", "mon_id": "LEU",
                "pdbx_pdb_strand_id": "A", "pdbx_pdb_ins_code": "N",
                "pdbx_seq_db_name": "UNP", "pdbx_seq_db_accession_code": "O26232",
                "db_mon_id": "MET", "pdbx_seq_db_seq_num": "1",
                "details": "SEE REMARK 999", "pdbx_auth_seq_num": "4",
            }, {
                "pdbx_pdb_id_code": "1LOL", "mon_id": "TYR",
                "pdbx_pdb_strand_id": "B", "pdbx_pdb_ins_code": "?",
                "pdbx_seq_db_name": "REF", "pdbx_seq_db_accession_code": "R12",
                "db_mon_id": "PRO", "pdbx_seq_db_seq_num": "?",
                "details": "?", "pdbx_auth_seq_num": "7",
            }]
        }
        self.assertEqual(create_seqadv_lines(mmcif), [
            "SEQADV 1LOL LEU A    4N UNP  O26232    MET     1 SEE REMARK 999",
            "SEQADV 1LOL TYR B    7  REF  R12       PRO",
        ])



class SeqresLinesTests(TestCase):

    def test_can_handle_no_table(self):
        self.assertEqual(create_seqres_lines({}), [])
    

    def test_can_create_lines(self):
        mmcif = {
            "entity_poly": [
                {"entity_id": "1", "pdbx_strand_id": "A,B"},
                {"entity_id": "2", "pdbx_strand_id": "C"},
            ],
            "entity_poly_seq": [
                {"entity_id": "1", "mon_id": "HIS"},
                {"entity_id": "1", "mon_id": "PRO"},
                {"entity_id": "1", "mon_id": "MET"},
                {"entity_id": "1", "mon_id": "TYR"},
                {"entity_id": "1", "mon_id": "HIS"},
                {"entity_id": "1", "mon_id": "MET"},
                {"entity_id": "1", "mon_id": "HIS"},
                {"entity_id": "1", "mon_id": "TRP"},
                {"entity_id": "1", "mon_id": "VAL"},
                {"entity_id": "1", "mon_id": "CYS"},
                {"entity_id": "1", "mon_id": "THR"},
                {"entity_id": "1", "mon_id": "VAL"},
                {"entity_id": "1", "mon_id": "PRO"},
                {"entity_id": "1", "mon_id": "PRO"},
                {"entity_id": "1", "mon_id": "ASP"},
                {"entity_id": "1", "mon_id": "TYR"},
                {"entity_id": "1", "mon_id": "ILE"},
                {"entity_id": "1", "mon_id": "TYR"},
                {"entity_id": "1", "mon_id": "HIS"},
                {"entity_id": "1", "mon_id": "MET"},
                {"entity_id": "1", "mon_id": "HIS"},
                {"entity_id": "1", "mon_id": "TRP"},
                {"entity_id": "1", "mon_id": "VAL"},
                {"entity_id": "1", "mon_id": "HIS"},
                {"entity_id": "1", "mon_id": "MET"},
                {"entity_id": "1", "mon_id": "HIS"},
                {"entity_id": "1", "mon_id": "TRP"},
                {"entity_id": "1", "mon_id": "VAL"},
                {"entity_id": "1", "mon_id": "TYR"},
                {"entity_id": "1", "mon_id": "HIS"},
                {"entity_id": "1", "mon_id": "MET"},
                {"entity_id": "2", "mon_id": "HIS"},
                {"entity_id": "2", "mon_id": "TRP"},
                {"entity_id": "2", "mon_id": "VAL"},
                {"entity_id": "2", "mon_id": "MET"},
                {"entity_id": "2", "mon_id": "HIS"},
                {"entity_id": "2", "mon_id": "HIS"},
                {"entity_id": "2", "mon_id": "VAL"},
                {"entity_id": "2", "mon_id": "TRP"},
            ]
        }
        self.assertEqual(create_seqres_lines(mmcif), [
            "SEQRES   1 A   31  HIS PRO MET TYR HIS MET HIS TRP VAL CYS THR VAL PRO",
            "SEQRES   2 A   31  PRO ASP TYR ILE TYR HIS MET HIS TRP VAL HIS MET HIS",
            "SEQRES   3 A   31  TRP VAL TYR HIS MET",
            "SEQRES   1 B   31  HIS PRO MET TYR HIS MET HIS TRP VAL CYS THR VAL PRO",
            "SEQRES   2 B   31  PRO ASP TYR ILE TYR HIS MET HIS TRP VAL HIS MET HIS",
            "SEQRES   3 B   31  TRP VAL TYR HIS MET",
            "SEQRES   1 C    8  HIS TRP VAL MET HIS HIS VAL TRP",
        ])



class ModresLinesTests(TestCase):

    def test_can_handle_no_table(self):
        self.assertEqual(create_modres_lines({}), [])
    

    def test_can_create_lines(self):
        mmcif = {
            "entry": [{"id": "1XXX"}],
            "pdbx_struct_mod_residue": [{
                "auth_asym_id": "T", "auth_comp_id": "4SU", "auth_seq_id": "8",
                "PDB_ins_code": "N", "parent_comp_id": "U", "details": "4-THIOURIDINE-5'-MONOPHOSPHATE",
            }, {
                "auth_asym_id": "L", "auth_comp_id": "5SU", "auth_seq_id": "12",
                "PDB_ins_code": "?", "parent_comp_id": "P", "details": "?",
            }]
        }
        self.assertEqual(create_modres_lines(mmcif), [
            "MODRES 1XXX 4SU T    8N   U  4-THIOURIDINE-5'-MONOPHOSPHATE",
            "MODRES 1XXX 5SU L   12    P"
        ])
    

    def test_can_create_lines_without_entry(self):
        mmcif = {
            "pdbx_struct_mod_residue": [{
                "auth_asym_id": "T", "auth_comp_id": "4SU", "auth_seq_id": "8",
                "PDB_ins_code": "N", "parent_comp_id": "U", "details": "4-THIOURIDINE-5'-MONOPHOSPHATE",
            }, {
                "auth_asym_id": "L", "auth_comp_id": "5SU", "auth_seq_id": "12",
                "PDB_ins_code": "?", "parent_comp_id": "P", "details": "?",
            }]
        }
        self.assertEqual(create_modres_lines(mmcif), [
            "MODRES      4SU T    8N   U  4-THIOURIDINE-5'-MONOPHOSPHATE",
            "MODRES      5SU L   12    P"
        ])



class HetLinesTests(TestCase):

    @patch("atomium.pdb.get_sig_counts")
    def test_can_handle_no_hets(self, mock_counts):
        mock_counts.return_value = {}
        self.assertEqual(create_het_lines({}), [])
        mock_counts.assert_called_with({})
    

    @patch("atomium.pdb.get_sig_counts")
    def test_can_handle_hets(self, mock_counts):
        mock_counts.return_value = {("A", "10", "", "XYZ"): 2, ("B", "12", "P", "BU2"): 3}
        self.assertEqual(create_het_lines("MMCIF"), [
            "HET    BU2  B  12P      3",
            "HET    XYZ  A  10       2",
        ])
        mock_counts.assert_called_with("MMCIF")



class HetnamLinesTests(TestCase):

    def test_can_handle_no_hetnams(self):
        self.assertEqual(create_hetnam_lines({}), [])
    

    @patch("atomium.pdb.split_lines")
    def test_can_get_get_hetnams(self, mock_split):
        mock_split.side_effect = [["The first name"], ["The second name", "over multiple", "lines"]]
        mmcif = {
            "entity": [
                {"id": "1", "type": "polymer"},
                {"id": "2", "type": "polymer"},
                {"id": "3", "type": "non-polymer"},
                {"id": "4", "type": "non-polymer"},
                {"id": "5", "type": "water"},
            ],
            "pdbx_entity_nonpoly": [
                {"entity_id": "3", "comp_id": "XYZ"},
                {"entity_id": "4", "comp_id": "BU2"},
                {"entity_id": "5", "comp_id": "HOH"},
            ],
            "chem_comp": [
                {"mon_nstd_flag": "y", "name": "modium", "id": "MOD"},
                {"mon_nstd_flag": "n", "name": "?", "id": "MYS"},
                {"mon_nstd_flag": ".", "name": "water", "id": "HOH"},
                {"mon_nstd_flag": "n", "name": "xanthosine", "id": "XYZ"},
                {"mon_nstd_flag": ".", "name": "butanol", "id": "BU2"},
                {"mon_nstd_flag": "n", "name": "mysterio", "id": "MY2"},
            ]
        }
        self.assertEqual(create_hetnam_lines(mmcif), [
            "HETNAM     XYZ The first name",
            "HETNAM     BU2 The second name",
            "HETNAM   2 BU2  over multiple",
            "HETNAM   3 BU2  lines",
        ])
        mock_split.assert_any_call("xanthosine", 55)
        mock_split.assert_any_call("butanol", 55)



class HetsynLinesTests(TestCase):

    def test_can_handle_no_hetsyns(self):
        self.assertEqual(create_hetsyn_lines({}), [])
    

    @patch("atomium.pdb.split_lines")
    def test_can_get_get_hetsyns(self, mock_split):
        mock_split.side_effect = [["The first synonym"], ["The second synonym", "over multiple", "lines"]]
        mmcif = {
            "entity": [
                {"id": "1", "type": "polymer"},
                {"id": "2", "type": "polymer"},
                {"id": "3", "type": "non-polymer"},
                {"id": "4", "type": "non-polymer"},
                {"id": "5", "type": "water"},
            ],
            "pdbx_entity_nonpoly": [
                {"entity_id": "3", "comp_id": "XYZ"},
                {"entity_id": "4", "comp_id": "BU2"},
                {"entity_id": "5", "comp_id": "HOH"},
            ],
            "chem_comp": [
                {"mon_nstd_flag": "y", "pdbx_synonyms": "modium", "id": "MOD"},
                {"mon_nstd_flag": "n", "pdbx_synonyms": "?", "id": "MYS"},
                {"mon_nstd_flag": ".", "pdbx_synonyms": "water", "id": "HOH"},
                {"mon_nstd_flag": "n", "pdbx_synonyms": "xanthosine, xanth", "id": "XYZ"},
                {"mon_nstd_flag": ".", "pdbx_synonyms": "butanol", "id": "BU2"},
                {"mon_nstd_flag": "n", "pdbx_synonyms": "mysterio", "id": "MY2"},
            ]
        }
        self.assertEqual(create_hetsyn_lines(mmcif), [
            "HETSYN     XYZ The first synonym",
            "HETSYN     BU2 The second synonym",
            "HETSYN   2 BU2  over multiple",
            "HETSYN   3 BU2  lines",
        ])
        mock_split.assert_any_call("xanthosine; xanth", 55)
        mock_split.assert_any_call("butanol", 55)



class FormulLinesTests(TestCase):

    @patch("atomium.pdb.get_sig_counts")
    @patch("atomium.pdb.split_lines")
    def test_can_make_formul_lines(self, mock_split, mock_sigs):
        mock_sigs.return_value = {
            ("A", "10", "", "BU2"): 3,
            ("A", "10", "", "HOH"): 5,
            ("A", "10", "", "MYS"): 2,
            ("A", "10", "", "XYZ"): 12,
        }
        mock_split.side_effect = [["formul1"], ["formul2"], ["formul3", "which needs", "3 lines"]]
        mmcif = {
            "entity": [
                {"id": "1", "type": "polymer"},
                {"id": "2", "type": "polymer"},
                {"id": "3", "type": "non-polymer"},
                {"id": "4", "type": "non-polymer"},
                {"id": "5", "type": "non-polymer"},
                {"id": "6", "type": "water"},
            ],
            "struct_asym": [
                {"entity_id": "1"},
                {"entity_id": "1"},
                {"entity_id": "2"},
                {"entity_id": "3"},
                {"entity_id": "4"},
                {"entity_id": "3"},
                {"entity_id": "4"},
                {"entity_id": "5"},
                {"entity_id": "6"},
                {"entity_id": "6"},
            ],
            "pdbx_entity_nonpoly": [
                {"entity_id": "3", "comp_id": "XYZ"},
                {"entity_id": "4", "comp_id": "BU2"},
                {"entity_id": "5", "comp_id": "MYS"},
                {"entity_id": "6", "comp_id": "HOH"},
            ],
            "chem_comp": [
                {"mon_nstd_flag": ".", "id": "MOD", "formula": "C2"},
                {"mon_nstd_flag": "n", "id": "MYS", "formula": "?"},
                {"mon_nstd_flag": ".", "id": "HOH", "formula": "C4"},
                {"mon_nstd_flag": "n", "id": "XYZ", "formula": "C5"},
                {"mon_nstd_flag": ".", "id": "BU2", "formula": "C6"},
                {"mon_nstd_flag": "n", "id": "MY2", "formula": "C7"},
            ]
        }
        self.assertEqual(create_formul_lines(mmcif), [
            "FORMUL   4  BU2    formul1",
            "FORMUL   6  HOH   *formul2",
            "FORMUL   3  XYZ    formul3",
            "FORMUL   3  XYZ  2 which needs",
            "FORMUL   3  XYZ  3 3 lines",
        ])
        mock_sigs.assert_called_with(mmcif, include_water=True, representative=True)
        mock_split.assert_any_call("5(C4)", 50)
        mock_split.assert_any_call("12(C5)", 50)
        mock_split.assert_any_call("3(C6)", 50)



class SigCountTests(TestCase):

    def setUp(self):
        self.mmcif = {
            "entity": [
                {"id": "1", "type": "polymer"},
                {"id": "2", "type": "non-polymer"},
                {"id": "3", "type": "non-polymer"},
                {"id": "4", "type": "water"},
            ],
            "chem_comp": [
                {"id": "HIS", "mon_nstd_flag": "y"},
                {"id": "XYZ", "mon_nstd_flag": "."},
                {"id": "BU2", "mon_nstd_flag": "."},
                {"id": "HOH", "mon_nstd_flag": "."},
            ],
            "atom_site": [
                {"auth_asym_id": "A", "auth_seq_id": "1", "pdbx_PDB_ins_code": "?", "auth_comp_id": "HIS", "label_entity_id": "1"},
                {"auth_asym_id": "A", "auth_seq_id": "1", "pdbx_PDB_ins_code": "?", "auth_comp_id": "HIS", "label_entity_id": "1"},
                {"auth_asym_id": "A", "auth_seq_id": "10", "pdbx_PDB_ins_code": "?", "auth_comp_id": "XYZ", "label_entity_id": "2"},
                {"auth_asym_id": "A", "auth_seq_id": "10", "pdbx_PDB_ins_code": "?", "auth_comp_id": "XYZ", "label_entity_id": "2"},
                {"auth_asym_id": "B", "auth_seq_id": "20", "pdbx_PDB_ins_code": "?", "auth_comp_id": "XYZ", "label_entity_id": "2"},
                {"auth_asym_id": "B", "auth_seq_id": "20", "pdbx_PDB_ins_code": "?", "auth_comp_id": "XYZ", "label_entity_id": "2"},
                {"auth_asym_id": "A", "auth_seq_id": "11", "pdbx_PDB_ins_code": "?", "auth_comp_id": "BU2", "label_entity_id": "3"},
                {"auth_asym_id": "A", "auth_seq_id": "11", "pdbx_PDB_ins_code": "?", "auth_comp_id": "BU2", "label_entity_id": "3"},
                {"auth_asym_id": "B", "auth_seq_id": "21", "pdbx_PDB_ins_code": "?", "auth_comp_id": "BU2", "label_entity_id": "3"},
                {"auth_asym_id": "B", "auth_seq_id": "21", "pdbx_PDB_ins_code": "?", "auth_comp_id": "BU2", "label_entity_id": "3"},
                {"auth_asym_id": "B", "auth_seq_id": "21", "pdbx_PDB_ins_code": "?", "auth_comp_id": "BU2", "label_entity_id": "3"},
                {"auth_asym_id": "B", "auth_seq_id": "21", "pdbx_PDB_ins_code": "A", "auth_comp_id": "BU2", "label_entity_id": "3"},
                {"auth_asym_id": "B", "auth_seq_id": "21", "pdbx_PDB_ins_code": "A", "auth_comp_id": "BU2", "label_entity_id": "3"},
                {"auth_asym_id": "A", "auth_seq_id": "100", "pdbx_PDB_ins_code": "?", "auth_comp_id": "HOH", "label_entity_id": "4"},
                {"auth_asym_id": "A", "auth_seq_id": "100", "pdbx_PDB_ins_code": "?", "auth_comp_id": "HOH", "label_entity_id": "4"},
                {"auth_asym_id": "A", "auth_seq_id": "101", "pdbx_PDB_ins_code": "?", "auth_comp_id": "HOH", "label_entity_id": "4"},
                {"auth_asym_id": "A", "auth_seq_id": "101", "pdbx_PDB_ins_code": "?", "auth_comp_id": "HOH", "label_entity_id": "4"},
                {"auth_asym_id": "B", "auth_seq_id": "200", "pdbx_PDB_ins_code": "?", "auth_comp_id": "HOH", "label_entity_id": "4"},
                {"auth_asym_id": "B", "auth_seq_id": "200", "pdbx_PDB_ins_code": "?", "auth_comp_id": "HOH", "label_entity_id": "4"},
                {"auth_asym_id": "B", "auth_seq_id": "201", "pdbx_PDB_ins_code": "?", "auth_comp_id": "HOH", "label_entity_id": "4"},
                {"auth_asym_id": "B", "auth_seq_id": "201", "pdbx_PDB_ins_code": "?", "auth_comp_id": "HOH", "label_entity_id": "4"},
            ]
        }

    def test_can_get_individual_molecules_without_water(self):
        self.assertEqual(get_sig_counts(self.mmcif), {
            ("A", "10", "", "XYZ"): 2,
            ("B", "20", "", "XYZ"): 2,
            ("A", "11", "", "BU2"): 2,
            ("B", "21", "", "BU2"): 3,
            ("B", "21", "A", "BU2"): 2,
        })
    

    def test_can_get_representative_molecules_with_water(self):
        self.assertEqual(get_sig_counts(self.mmcif, include_water=True, representative=True), {
            ("A", "10", "", "XYZ"): 1,
            ("B", "20", "", "XYZ"): 1,
            ("A", "11", "", "BU2"): 1,
            ("B", "21", "", "BU2"): 1,
            ("B", "21", "A", "BU2"): 1,
            ("A", "100", "", "HOH"): 1,
            ("A", "101", "", "HOH"): 1,
            ("B", "200", "", "HOH"): 1,
            ("B", "201", "", "HOH"): 1,
        })



class HelixLinesTests(TestCase):

    def test_can_handle_no_helix(self):
        self.assertEqual(create_helix_lines({}), [])
    

    def test_can_produce_helix_lines(self):
        mmcif = {
            "struct_conf": [{
                "conf_type_id": "HELX_P", "id": "HELX_P1", "pdbx_PDB_helix_id": "1",
                "beg_label_comp_id": "ILE", "beg_label_asym_id": "A",
                "beg_label_seq_id": "2", "pdbx_beg_PDB_ins_code": "?",
                "end_label_comp_id": "CYS", "end_label_asym_id": "A",
                "end_label_seq_id": "6", "pdbx_end_PDB_ins_code": "?",
                "beg_auth_comp_id": "ILE", "beg_auth_asym_id": "A",
                "beg_auth_seq_id": "2", "end_auth_comp_id": "CYS",
                "end_auth_asym_id": "A", "end_auth_seq_id": "6",
                "pdbx_PDB_helix_class": "1", "details": "?",
                "pdbx_PDB_helix_length": "5",
            }, {
                "conf_type_id": "HELX_P", "id": "HELX_P12", "pdbx_PDB_helix_id": "12",
                "beg_label_comp_id": "VAL", "beg_label_asym_id": "H",
                "beg_label_seq_id": "2", "pdbx_beg_PDB_ins_code": "P",
                "end_label_comp_id": "ARG", "end_label_asym_id": "H",
                "end_label_seq_id": "22", "pdbx_end_PDB_ins_code": "L",
                "beg_auth_comp_id": "VAL", "beg_auth_asym_id": "H",
                "beg_auth_seq_id": "2", "end_auth_comp_id": "ARG",
                "end_auth_asym_id": "H", "end_auth_seq_id": "22",
                "pdbx_PDB_helix_class": "1", "details": "?",
                "pdbx_PDB_helix_length": "21",
            }]
        }
        self.assertEqual(create_helix_lines(mmcif), [
            "HELIX    1   1 ILE A    2  CYS A    6  1                                   5",
            "HELIX   12  12 VAL H    2P ARG H   22L 1                                  21",
        ])



class SheetLineTests(TestCase):

    def test_can_handle_no_range(self):
        self.assertEqual(create_sheet_lines({}), [])
    

    def test_can_create_sheet_lines(self):
        mmcif = {
            "struct_sheet_order": [
                {"sheet_id": "A", "range_id_1": "1", "range_id_2": "2", "sense": "parallel"},
                {"sheet_id": "B", "range_id_1": "1", "range_id_2": "2", "sense": "anti-parallel"},
            ],
            "pdbx_struct_sheet_hbond": [{
                "sheet_id": "A", "range_id_1": "1", "range_id_2": "2",
                "range_1_PDB_ins_code": "?", "range_1_auth_atom_id": "O",
                "range_1_auth_comp_id": "PHE", "range_1_auth_asym_id": "B",
                "range_1_auth_seq_id": "24", "range_2_PDB_ins_code": "?",
                "range_2_auth_atom_id": "N", "range_2_auth_comp_id": "TYR",
                "range_2_auth_asym_id": "D", "range_2_auth_seq_id": "26",
            }, {
                "sheet_id": "B", "range_id_1": "1", "range_id_2": "2",
                "range_1_PDB_ins_code": "Z", "range_1_auth_atom_id": "O",
                "range_1_auth_comp_id": "PHE", "range_1_auth_asym_id": "F",
                "range_1_auth_seq_id": "24", "range_2_PDB_ins_code": "Q",
                "range_2_auth_atom_id": "N", "range_2_auth_comp_id": "TYR",
                "range_2_auth_asym_id": "H", "range_2_auth_seq_id": "26",
            }],
            "struct_sheet_range": [{
                "sheet_id": "A", "id": "1", "beg_auth_comp_id": "PHE", "beg_auth_asym_id": "B",
                "beg_auth_seq_id": "24", "pdbx_beg_PDB_ins_code": "?", "end_auth_comp_id": "TYR",
                "end_auth_asym_id": "B", "end_auth_seq_id": "26", "pdbx_end_PDB_ins_code": "?"
            }, {
                "sheet_id": "A", "id": "2", "beg_auth_comp_id": "PHE", "beg_auth_asym_id": "D",
                "beg_auth_seq_id": "24", "pdbx_beg_PDB_ins_code": "?", "end_auth_comp_id": "TYR",
                "end_auth_asym_id": "D", "end_auth_seq_id": "26", "pdbx_end_PDB_ins_code": "?"
            }, {
                "sheet_id": "B", "id": "1", "beg_auth_comp_id": "PHE", "beg_auth_asym_id": "F",
                "beg_auth_seq_id": "24", "pdbx_beg_PDB_ins_code": "?", "end_auth_comp_id": "TYR",
                "end_auth_asym_id": "F", "end_auth_seq_id": "26", "pdbx_end_PDB_ins_code": "?"
            }, {
                "sheet_id": "B", "id": "2", "beg_auth_comp_id": "PHE", "beg_auth_asym_id": "H",
                "beg_auth_seq_id": "24", "pdbx_beg_PDB_ins_code": "L", "end_auth_comp_id": "TYR",
                "end_auth_asym_id": "H", "end_auth_seq_id": "26", "pdbx_end_PDB_ins_code": "N"
            }]
        }
        self.assertEqual(create_sheet_lines(mmcif), [
            "SHEET    1   A 2 PHE B  24  TYR B  26  0",
            "SHEET    2   A 2 PHE D  24  TYR D  26  1  N  TYR D  26   O  PHE B  24",
            "SHEET    1   B 2 PHE F  24  TYR F  26  0",
            "SHEET    2   B 2 PHE H  24L TYR H  26N-1  N  TYR H  26Q  O  PHE F  24Z"
        ])



class SsbondLinesTests(TestCase):

    def test_can_handle_no_ssbond(self):
        self.assertEqual(create_ssbond_lines({}), [])
    

    def test_can_save_ssbond_lines(self):
        mmcif = {
            "struct_conn": [{
                "id": "disulf1", "conn_type_id": "disulf",
                "pdbx_ptnr1_PDB_ins_code": "?", "ptnr1_symmetry": "1555",
                "pdbx_ptnr2_PDB_ins_code": "?", "ptnr1_auth_asym_id": "A",
                "ptnr1_auth_comp_id": "CYS", "ptnr1_auth_seq_id": "6",
                "ptnr2_auth_asym_id": "A", "ptnr2_auth_comp_id": "CYS",
                "ptnr2_auth_seq_id": "11", "ptnr2_symmetry": "1555",
                "pdbx_dist_value": "2.02",
            }, {
                "id": "disulf1", "conn_type_id": "disulf", 
                "pdbx_ptnr1_PDB_ins_code": "P", "pdbx_ptnr1_standard_comp_id": "?",
                "ptnr1_symmetry": "1554", "pdbx_ptnr2_PDB_ins_code": "Z",
                "ptnr1_auth_asym_id": "B", "ptnr1_auth_comp_id": "MET",
                "ptnr1_auth_seq_id": "99", "ptnr2_auth_asym_id": "C",
                "ptnr2_auth_comp_id": "THR", "ptnr2_auth_seq_id": "101",
                "ptnr2_symmetry": "1556", "pdbx_dist_value": "2.55",
            }, {
                "id": "covale29", "conn_type_id": "covale",
                "pdbx_leaving_atom_flag": "?", "pdbx_PDB_id": "?",
                "ptnr1_label_asym_id": "V", "ptnr1_label_comp_id": "ZN",
                "ptnr1_label_seq_id": ".", "ptnr1_label_atom_id": "ZN",
            }]
        }
        self.assertEqual(create_ssbond_lines(mmcif), [
            "SSBOND   1 CYS A    6    CYS A   11                          1555   1555  2.02",
            "SSBOND   2 MET B   99P   THR C  101Z                         1554   1556  2.55"
        ])



class LinkLinesTests(TestCase):

    def test_can_handle_no_links(self):
        self.assertEqual(create_link_lines({}), [])
    

    def test_can_produce_link_lines(self):
        mmcif = {
            "struct_conn": [{
                "id": "covale1", "conn_type_id": "covale", "ptnr1_label_atom_id": "C",
                "pdbx_ptnr1_label_alt_id": "?", "pdbx_ptnr1_PDB_ins_code": "?",
                "ptnr1_symmetry": "1555", "ptnr2_label_atom_id": "N",
                "pdbx_ptnr2_label_alt_id": "?", "pdbx_ptnr2_PDB_ins_code": "?",
                "ptnr1_auth_asym_id": "A", "ptnr1_auth_comp_id": "ACE",
                "ptnr1_auth_seq_id": "0", "ptnr2_auth_asym_id": "A",
                "ptnr2_auth_comp_id": "HIS", "ptnr2_auth_seq_id": "1",
                "ptnr2_symmetry": "1555", "pdbx_dist_value": "1.38",
            }, {
                "id": "covale7", "conn_type_id": "covale", "ptnr1_label_atom_id": "CE",
                "pdbx_ptnr1_label_alt_id": "?", "pdbx_ptnr1_PDB_ins_code": "E",
                "ptnr1_symmetry": "1554", "ptnr2_label_atom_id": "CE",
                "pdbx_ptnr2_label_alt_id": "?", "pdbx_ptnr2_PDB_ins_code": "K",
                "ptnr1_auth_asym_id": "B", "ptnr1_auth_comp_id": "MK8",
                "ptnr1_auth_seq_id": "2", "ptnr2_auth_asym_id": "C",
                "ptnr2_auth_comp_id": "MK8", "ptnr2_auth_seq_id": "6",
                "ptnr2_symmetry": "1556", "pdbx_dist_value": "1.39",
            }, {
                "id": "covale29", "conn_type_id": "disulf",
                "pdbx_leaving_atom_flag": "?", "pdbx_PDB_id": "?",
                "ptnr1_label_asym_id": "V", "ptnr1_label_comp_id": "ZN",
                "ptnr1_label_seq_id": ".", "ptnr1_label_atom_id": "ZN",
            }]
        }
        self.assertEqual(create_link_lines(mmcif), [
            "LINK           C ACE A   0                   N HIS A   1     1555   1555  1.38",
            "LINK          CE MK8 B   2E                 CE MK8 C   6K    1554   1556  1.39"
        ])



class CispepLinesTests(TestCase):

    def test_can_handle_no_cispep(self):
        self.assertEqual(create_cispep_lines({}), [])
    

    def test_can_save_cispep_lines(self):
        mmcif = {
            "struct_mon_prot_cis": [{
                "pdbx_id": "1",  "pdbx_PDB_ins_code": "?", "auth_comp_id": "ASN",
                "auth_seq_id": "77", "auth_asym_id": "A", "pdbx_PDB_ins_code_2": "?",
                "pdbx_auth_comp_id_2": "PRO", "pdbx_auth_seq_id_2": "78",
                "pdbx_auth_asym_id_2": "A", "pdbx_PDB_model_num": "1",
                "pdbx_omega_angle": "1.58",
            }, {
                "pdbx_id": "2", "pdbx_PDB_ins_code": "N", "auth_comp_id": "ASN",
                "auth_seq_id": "88", "auth_asym_id": "C", "pdbx_PDB_ins_code_2": "K",
                "pdbx_auth_comp_id_2": "PRO", "pdbx_auth_seq_id_2": "79",
                "pdbx_auth_asym_id_2": "C", "pdbx_PDB_model_num": "1",
                "pdbx_omega_angle": "-1.32",
            }]
        }
        self.assertEqual(create_cispep_lines(mmcif), [
            "CISPEP   1 ASN A   77    PRO A   78          0         1.58",
            "CISPEP   2 ASN C   88N   PRO C   79K         0        -1.32"
        ])



class SiteLinesTests(TestCase):

    def test_can_handle_no_site(self):
        self.assertEqual(create_site_lines({}), [])
    

    def test_can_produce_site_lines(self):
        mmcif = {
            "struct_site_gen": [{
                "site_id": "AC1", "pdbx_num_res": "6", "pdbx_auth_ins_code": "?",
                "auth_comp_id": "HOH", "auth_asym_id": "A", "auth_seq_id": "1577"
            }, {
                "site_id": "AC1", "pdbx_num_res": "6", "pdbx_auth_ins_code": "P",
                "auth_comp_id": "PRO", "auth_asym_id": "A", "auth_seq_id": "1578"
            }, {
                "site_id": "AC1", "pdbx_num_res": "6", "pdbx_auth_ins_code": "?",
                "auth_comp_id": "VAL", "auth_asym_id": "A", "auth_seq_id": "1579"
            }, {
               "site_id": "AC1", "pdbx_num_res": "6", "pdbx_auth_ins_code": "?",
                "auth_comp_id": "CYS", "auth_asym_id": "A", "auth_seq_id": "1560"
            }, {
                "site_id": "AC1", "pdbx_num_res": "6", "pdbx_auth_ins_code": "?",
                "auth_comp_id": "GLU", "auth_asym_id": "B", "auth_seq_id": "1561"
            }, {
                "site_id": "AC1", "pdbx_num_res": "6", "pdbx_auth_ins_code": "?",
                "auth_comp_id": "ASP", "auth_asym_id": "B", "auth_seq_id": "8"
            }, {
                "site_id": "AC2", "pdbx_num_res": "2", "pdbx_auth_ins_code": "?",
                "auth_comp_id": "HOH", "auth_asym_id": "A", "auth_seq_id": "22"
            }, {
                "site_id": "AC2", "pdbx_num_res": "2", "pdbx_auth_ins_code": "?",
                "auth_comp_id": "GLN", "auth_asym_id": "C", "auth_seq_id": "33"
            }]
        }
        self.assertEqual(create_site_lines(mmcif), [
            "SITE     1 AC1  6 HOH A1577  PRO A1578P VAL A1579  CYS A1560",
            "SITE     2 AC1  6 GLU B1561  ASP B   8",
            "SITE     1 AC2  2 HOH A  22  GLN C  33"
        ])



class Cryst1LinesTests(TestCase):

    def test_can_handle_no_tables(self):
        self.assertEqual(create_cryst1_line({}), [])
    

    def test_can_handle_missing_values(self):
        mmcif = {
            "cell": [{
                "entry_id": "?", "length_a": "?", "length_b": "?",
                "length_c": "?", "angle_alpha": "?", "angle_beta": "?",
                "angle_gamma": "?", "Z_pdb": "?", "pdbx_unique_axis": "?",
            }],
            "symmetry": [{"space_group_name_H-M": "?"}]
        }
        self.assertEqual(create_cryst1_line(mmcif), [])

    
    def test_can_handle_some_missing_values(self):
        mmcif = {
            "cell": [{
                "entry_id": "?", "length_a": "57.50", "length_b": "?",
                "length_c": "?", "angle_alpha": "?", "angle_beta": "?",
                "angle_gamma": "?", "Z_pdb": "?", "pdbx_unique_axis": "?",
            }],
            "symmetry": [{"space_group_name_H-M": "?"}]
        }
        self.assertEqual(create_cryst1_line(mmcif), ["CRYST1    57.50"])
    

    def test_can_produce_cryst1_lines(self):
        mmcif = {
            "cell": [{
                "entry_id": "1LOL", "length_a": "57.570", "length_b": "55.482",
                "length_c": "66.129", "angle_alpha": "90.00", "angle_beta": "94.28",
                "angle_gamma": "45.00", "Z_pdb": "4"
            }],
            "symmetry": [{"space_group_name_H-M": "P 1 21 1"}]
        }
        self.assertEqual(create_cryst1_line(mmcif), [
            "CRYST1   57.570   55.482   66.129  90.00  94.28  45.00 P 1 21 1      4"
        ])



class OrigxLinesTests(TestCase):

    def test_can_handle_no_table(self):
        self.assertEqual(create_origxn_lines({}), [])
    

    def test_can_create_origx_lines(self):
        mmcif = {
            "database_PDB_matrix": [{
                "entry_id": "1LOL", "origx[1][1]": "1.000000",
                "origx[1][2]": "0.000001", "origx[1][3]": "0.000002",
                "origx[2][1]": "0.000003", "origx[2][2]": "9.000000",
                "origx[2][3]": "0.000005", "origx[3][1]": "0.000004",
                "origx[3][2]": "0.000007", "origx[3][3]": "10.000000",
                "origx_vector[1]": "0.40000", "origx_vector[2]": "0.05000",
                "origx_vector[3]": "0.00060",
            }]
        }
        self.assertEqual(create_origxn_lines(mmcif), [
            "ORIGX1      1.000000  0.000001  0.000002        0.40000",
            "ORIGX2      0.000003  9.000000  0.000005        0.05000",
            "ORIGX3      0.000004  0.000007 10.000000        0.00060"
        ])



class ScalenLinesTests(TestCase):

    def test_can_handle_no_table(self):
        self.assertEqual(create_scalen_lines({}), [])
    

    def test_can_create_scalen_lines(self):
        mmcif = {
            "atom_sites": [{
                "entry_id": "1LOL", "fract_transf_matrix[1][1]": "0.017370",
                "fract_transf_matrix[1][2]": "0.000000",
                "fract_transf_matrix[1][3]": "0.001301",
                "fract_transf_matrix[2][1]": "0.000000",
                "fract_transf_matrix[2][2]": "0.018024",
                "fract_transf_matrix[2][3]": "0.000000",
                "fract_transf_matrix[3][1]": "0.040003",
                "fract_transf_matrix[3][2]": "0.000090",
                "fract_transf_matrix[3][3]": "0.015164",
                "fract_transf_vector[1]": "0.10000",
                "fract_transf_vector[2]": "0.20000",
                "fract_transf_vector[3]": "0.30000",
            }]
        }
        self.assertEqual(create_scalen_lines(mmcif), [
            "SCALE1      0.017370  0.000000  0.001301        0.10000",
            "SCALE2      0.000000  0.018024  0.000000        0.20000",
            "SCALE3      0.040003  0.000090  0.015164        0.30000"
        ])



class MtrixnLinesTests(TestCase):

    def test_can_handle_no_table(self):
        self.assertEqual(create_mtrixn_lines({}), [])
    

    def test_can_create_mtrixn_lines(self):
        mmcif = {
            "struct_ncs_oper": [{
                "id": "1", "code": "given", "details": "?",
                "matrix[1][1]": "-0.745000", "matrix[1][2]": "0.008000", "matrix[1][3]": "0.667000",
                "matrix[2][1]": "0.020000", "matrix[2][2]": "-0.999000", "matrix[2][3]": "0.035000",
                "matrix[3][1]": "0.667000", "matrix[3][2]": "0.040000", "matrix[3][3]": "0.744000",
                "vector[1]": "14.56700", "vector[2]": "27.10700", "vector[3]": "-6.04300"
            }, {
                "id": "2", "code": "?", "details": "?",
                "matrix[1][1]": "-0.746000", "matrix[1][2]": "-0.029000", "matrix[1][3]": "0.665000",
                "matrix[2][1]": "0.045000", "matrix[2][2]": "-0.999000", "matrix[2][3]": "0.007000",
                "matrix[3][1]": "0.664000", "matrix[3][2]": "0.035000", "matrix[3][3]": "0.747000",
                "vector[1]": "14.68300", "vector[2]": "27.24500", "vector[3]": "-5.93200"
            }]
        }
        self.assertEqual(create_mtrixn_lines(mmcif), [
            "MTRIX1   1 -0.745000  0.008000  0.667000       14.56700    1",
            "MTRIX2   1  0.020000 -0.999000  0.035000       27.10700    1",
            "MTRIX3   1  0.667000  0.040000  0.744000       -6.04300    1",
            "MTRIX1   2 -0.746000 -0.029000  0.665000       14.68300     ",
            "MTRIX2   2  0.045000 -0.999000  0.007000       27.24500     ",
            "MTRIX3   2  0.664000  0.035000  0.747000       -5.93200     "
        ])



class AtomLinesTests(TestCase):

    def test_can_handle_no_table(self):
        self.assertEqual(create_atom_lines({}), [])
    

    @patch("atomium.pdb.create_atom_line")
    @patch("atomium.pdb.create_aniso_line")
    def test_can_get_lines(self, mock_aniso, mock_atom):
        mock_atom.side_effect = lambda a, i: "ATOM " + str(i) + " " + a["label_asym_id"] * 15
        mock_aniso.return_value = None
        mmcif = {
            "entity": [
                {"id": "1", "type": "polymer"},
                {"id": "2", "type": "polymer"},
                {"id": "3", "type": "non-polymer"},
                {"id": "4", "type": "water"},
            ],
            "atom_site": [
                {"id": "1", "pdbx_PDB_model_num": "1", "label_asym_id": "A", "label_entity_id": "1"},
                {"id": "2", "pdbx_PDB_model_num": "1", "label_asym_id": "A", "label_entity_id": "1"},
                {"id": "3", "pdbx_PDB_model_num": "1", "label_asym_id": "B", "label_entity_id": "1"},
                {"id": "4", "pdbx_PDB_model_num": "1", "label_asym_id": "B", "label_entity_id": "1"},
                {"id": "5", "pdbx_PDB_model_num": "1", "label_asym_id": "C", "label_entity_id": "2"},
                {"id": "6", "pdbx_PDB_model_num": "1", "label_asym_id": "C", "label_entity_id": "2"},
                {"id": "7", "pdbx_PDB_model_num": "1", "label_asym_id": "D", "label_entity_id": "3"},
                {"id": "8", "pdbx_PDB_model_num": "1", "label_asym_id": "D", "label_entity_id": "3"},
                {"id": "9", "pdbx_PDB_model_num": "1", "label_asym_id": "E", "label_entity_id": "4"},
                {"id": "10", "pdbx_PDB_model_num": "1", "label_asym_id": "F", "label_entity_id": "4"}
            ]
        }
        self.assertEqual(create_atom_lines(mmcif), [
            "ATOM 1 AAAAAAAAAAAAAAA", 
            "ATOM 2 AAAAAAAAAAAAAAA", 
            "TER       3      AAAAA",
            "ATOM 4 BBBBBBBBBBBBBBB", 
            "ATOM 5 BBBBBBBBBBBBBBB", 
            "TER       6      BBBBB",
            "ATOM 7 CCCCCCCCCCCCCCC", 
            "ATOM 8 CCCCCCCCCCCCCCC", 
            "TER       9      CCCCC",
            "ATOM 10 DDDDDDDDDDDDDDD",
            "ATOM 11 DDDDDDDDDDDDDDD",
            "ATOM 12 EEEEEEEEEEEEEEE",
            "ATOM 13 FFFFFFFFFFFFFFF",
        ])
        self.assertEqual([c[0] for c in mock_atom.call_args_list], [
            ({"id": "1", "pdbx_PDB_model_num": "1", "label_asym_id": "A", "label_entity_id": "1"}, 1),
            ({"id": "2", "pdbx_PDB_model_num": "1", "label_asym_id": "A", "label_entity_id": "1"}, 2),
            ({"id": "3", "pdbx_PDB_model_num": "1", "label_asym_id": "B", "label_entity_id": "1"}, 4),
            ({"id": "4", "pdbx_PDB_model_num": "1", "label_asym_id": "B", "label_entity_id": "1"}, 5),
            ({"id": "5", "pdbx_PDB_model_num": "1", "label_asym_id": "C", "label_entity_id": "2"}, 7),
            ({"id": "6", "pdbx_PDB_model_num": "1", "label_asym_id": "C", "label_entity_id": "2"}, 8),
            ({"id": "7", "pdbx_PDB_model_num": "1", "label_asym_id": "D", "label_entity_id": "3"}, 10),
            ({"id": "8", "pdbx_PDB_model_num": "1", "label_asym_id": "D", "label_entity_id": "3"}, 11),
            ({"id": "9", "pdbx_PDB_model_num": "1", "label_asym_id": "E", "label_entity_id": "4"}, 12),
            ({"id": "10", "pdbx_PDB_model_num": "1", "label_asym_id": "F", "label_entity_id": "4"}, 13),
        ])
        self.assertEqual([c[0] for c in mock_aniso.call_args_list], [
            ({"id": "1", "pdbx_PDB_model_num": "1", "label_asym_id": "A", "label_entity_id": "1"}, None, 1),
            ({"id": "2", "pdbx_PDB_model_num": "1", "label_asym_id": "A", "label_entity_id": "1"}, None, 2),
            ({"id": "3", "pdbx_PDB_model_num": "1", "label_asym_id": "B", "label_entity_id": "1"}, None, 4),
            ({"id": "4", "pdbx_PDB_model_num": "1", "label_asym_id": "B", "label_entity_id": "1"}, None, 5),
            ({"id": "5", "pdbx_PDB_model_num": "1", "label_asym_id": "C", "label_entity_id": "2"}, None, 7),
            ({"id": "6", "pdbx_PDB_model_num": "1", "label_asym_id": "C", "label_entity_id": "2"}, None, 8),
            ({"id": "7", "pdbx_PDB_model_num": "1", "label_asym_id": "D", "label_entity_id": "3"}, None, 10),
            ({"id": "8", "pdbx_PDB_model_num": "1", "label_asym_id": "D", "label_entity_id": "3"}, None, 11),
            ({"id": "9", "pdbx_PDB_model_num": "1", "label_asym_id": "E", "label_entity_id": "4"}, None, 12),
            ({"id": "10", "pdbx_PDB_model_num": "1", "label_asym_id": "F", "label_entity_id": "4"}, None, 13),
        ])
    

    @patch("atomium.pdb.create_atom_line")
    @patch("atomium.pdb.create_aniso_line")
    def test_can_get_lines_with_aniso(self, mock_aniso, mock_atom):
        mock_atom.side_effect = lambda a, i: "ATOM " + str(i) + " " + a["label_asym_id"] * 15
        mock_aniso.side_effect = lambda a, _, i: "ANISO " + str(i) + " " + a["label_asym_id"] * 15
        mmcif = {
            "entity": [
                {"id": "1", "type": "polymer"},
                {"id": "2", "type": "polymer"},
                {"id": "3", "type": "non-polymer"},
                {"id": "4", "type": "water"},
            ],
            "atom_site": [
                {"id": "1", "pdbx_PDB_model_num": "1", "label_asym_id": "A", "label_entity_id": "1"},
                {"id": "2", "pdbx_PDB_model_num": "1", "label_asym_id": "A", "label_entity_id": "1"},
                {"id": "3", "pdbx_PDB_model_num": "1", "label_asym_id": "B", "label_entity_id": "1"},
                {"id": "4", "pdbx_PDB_model_num": "1", "label_asym_id": "B", "label_entity_id": "1"},
                {"id": "5", "pdbx_PDB_model_num": "1", "label_asym_id": "C", "label_entity_id": "2"},
                {"id": "6", "pdbx_PDB_model_num": "1", "label_asym_id": "C", "label_entity_id": "2"},
                {"id": "7", "pdbx_PDB_model_num": "1", "label_asym_id": "D", "label_entity_id": "3"},
                {"id": "8", "pdbx_PDB_model_num": "1", "label_asym_id": "D", "label_entity_id": "3"},
                {"id": "9", "pdbx_PDB_model_num": "1", "label_asym_id": "E", "label_entity_id": "4"},
                {"id": "10", "pdbx_PDB_model_num": "1", "label_asym_id": "F", "label_entity_id": "4"}
            ],
            "atom_site_anisotrop": [
                {"id": "1", "aniso": "10"},
                {"id": "2", "aniso": "20"},
                {"id": "8", "aniso": "80"},
            ]
        }
        self.assertEqual(create_atom_lines(mmcif), [
            "ATOM 1 AAAAAAAAAAAAAAA", "ANISO 1 AAAAAAAAAAAAAAA",
            "ATOM 2 AAAAAAAAAAAAAAA", "ANISO 2 AAAAAAAAAAAAAAA",
            "TER       3      AAAAAA",
            "ATOM 4 BBBBBBBBBBBBBBB", "ANISO 4 BBBBBBBBBBBBBBB",
            "ATOM 5 BBBBBBBBBBBBBBB", "ANISO 5 BBBBBBBBBBBBBBB",
            "TER       6      BBBBBB",
            "ATOM 7 CCCCCCCCCCCCCCC", "ANISO 7 CCCCCCCCCCCCCCC",
            "ATOM 8 CCCCCCCCCCCCCCC", "ANISO 8 CCCCCCCCCCCCCCC",
            "TER       9      CCCCCC",
            "ATOM 10 DDDDDDDDDDDDDDD", "ANISO 10 DDDDDDDDDDDDDDD",
            "ATOM 11 DDDDDDDDDDDDDDD", "ANISO 11 DDDDDDDDDDDDDDD",
            "ATOM 12 EEEEEEEEEEEEEEE", "ANISO 12 EEEEEEEEEEEEEEE",
            "ATOM 13 FFFFFFFFFFFFFFF", "ANISO 13 FFFFFFFFFFFFFFF",
        ])
        self.assertEqual([c[0] for c in mock_atom.call_args_list], [
            ({"id": "1", "pdbx_PDB_model_num": "1", "label_asym_id": "A", "label_entity_id": "1"}, 1),
            ({"id": "2", "pdbx_PDB_model_num": "1", "label_asym_id": "A", "label_entity_id": "1"}, 2),
            ({"id": "3", "pdbx_PDB_model_num": "1", "label_asym_id": "B", "label_entity_id": "1"}, 4),
            ({"id": "4", "pdbx_PDB_model_num": "1", "label_asym_id": "B", "label_entity_id": "1"}, 5),
            ({"id": "5", "pdbx_PDB_model_num": "1", "label_asym_id": "C", "label_entity_id": "2"}, 7),
            ({"id": "6", "pdbx_PDB_model_num": "1", "label_asym_id": "C", "label_entity_id": "2"}, 8),
            ({"id": "7", "pdbx_PDB_model_num": "1", "label_asym_id": "D", "label_entity_id": "3"}, 10),
            ({"id": "8", "pdbx_PDB_model_num": "1", "label_asym_id": "D", "label_entity_id": "3"}, 11),
            ({"id": "9", "pdbx_PDB_model_num": "1", "label_asym_id": "E", "label_entity_id": "4"}, 12),
            ({"id": "10", "pdbx_PDB_model_num": "1", "label_asym_id": "F", "label_entity_id": "4"}, 13),
        ])
        self.assertEqual([c[0] for c in mock_aniso.call_args_list], [
            ({"id": "1", "pdbx_PDB_model_num": "1", "label_asym_id": "A", "label_entity_id": "1"}, {"id": "1", "aniso": "10"}, 1),
            ({"id": "2", "pdbx_PDB_model_num": "1", "label_asym_id": "A", "label_entity_id": "1"}, {"id": "2", "aniso": "20"}, 2),
            ({"id": "3", "pdbx_PDB_model_num": "1", "label_asym_id": "B", "label_entity_id": "1"}, None, 4),
            ({"id": "4", "pdbx_PDB_model_num": "1", "label_asym_id": "B", "label_entity_id": "1"}, None, 5),
            ({"id": "5", "pdbx_PDB_model_num": "1", "label_asym_id": "C", "label_entity_id": "2"}, None, 7),
            ({"id": "6", "pdbx_PDB_model_num": "1", "label_asym_id": "C", "label_entity_id": "2"}, None, 8),
            ({"id": "7", "pdbx_PDB_model_num": "1", "label_asym_id": "D", "label_entity_id": "3"}, None, 10),
            ({"id": "8", "pdbx_PDB_model_num": "1", "label_asym_id": "D", "label_entity_id": "3"}, {"id": "8", "aniso": "80"}, 11),
            ({"id": "9", "pdbx_PDB_model_num": "1", "label_asym_id": "E", "label_entity_id": "4"}, None, 12),
            ({"id": "10", "pdbx_PDB_model_num": "1", "label_asym_id": "F", "label_entity_id": "4"}, None, 13),
        ])
    

    @patch("atomium.pdb.create_atom_line")
    @patch("atomium.pdb.create_aniso_line")
    def test_can_get_lines_with_models(self, mock_aniso, mock_atom):
        mock_atom.side_effect = lambda a, i: "ATOM " + str(i) + " " + a["label_asym_id"] * 15
        mock_aniso.return_value = None
        mmcif = {
            "entity": [
                {"id": "1", "type": "polymer"},
                {"id": "2", "type": "polymer"},
                {"id": "3", "type": "non-polymer"},
                {"id": "4", "type": "water"},
            ],
            "atom_site": [
                {"id": "1", "pdbx_PDB_model_num": "1", "label_asym_id": "A", "label_entity_id": "1"},
                {"id": "2", "pdbx_PDB_model_num": "1", "label_asym_id": "A", "label_entity_id": "1"},
                {"id": "3", "pdbx_PDB_model_num": "1", "label_asym_id": "B", "label_entity_id": "1"},
                {"id": "4", "pdbx_PDB_model_num": "1", "label_asym_id": "B", "label_entity_id": "1"},
                {"id": "7", "pdbx_PDB_model_num": "1", "label_asym_id": "D", "label_entity_id": "3"},
                {"id": "8", "pdbx_PDB_model_num": "1", "label_asym_id": "D", "label_entity_id": "3"},
                {"id": "1", "pdbx_PDB_model_num": "2", "label_asym_id": "A", "label_entity_id": "1"},
                {"id": "2", "pdbx_PDB_model_num": "2", "label_asym_id": "A", "label_entity_id": "1"},
                {"id": "3", "pdbx_PDB_model_num": "2", "label_asym_id": "B", "label_entity_id": "1"},
                {"id": "4", "pdbx_PDB_model_num": "2", "label_asym_id": "B", "label_entity_id": "1"},
                {"id": "7", "pdbx_PDB_model_num": "2", "label_asym_id": "D", "label_entity_id": "3"},
                {"id": "8", "pdbx_PDB_model_num": "2", "label_asym_id": "D", "label_entity_id": "3"},
                {"id": "1", "pdbx_PDB_model_num": "3", "label_asym_id": "A", "label_entity_id": "1"},
                {"id": "2", "pdbx_PDB_model_num": "3", "label_asym_id": "A", "label_entity_id": "1"},
                {"id": "3", "pdbx_PDB_model_num": "3", "label_asym_id": "B", "label_entity_id": "1"},
                {"id": "4", "pdbx_PDB_model_num": "3", "label_asym_id": "B", "label_entity_id": "1"},
                {"id": "7", "pdbx_PDB_model_num": "3", "label_asym_id": "D", "label_entity_id": "3"},
                {"id": "8", "pdbx_PDB_model_num": "3", "label_asym_id": "D", "label_entity_id": "3"},
            ]
        }
        self.assertEqual(create_atom_lines(mmcif), [
            "MODEL        1",
            "ATOM 1 AAAAAAAAAAAAAAA", 
            "ATOM 2 AAAAAAAAAAAAAAA", 
            "TER       3      AAAAA",
            "ATOM 4 BBBBBBBBBBBBBBB", 
            "ATOM 5 BBBBBBBBBBBBBBB", 
            "TER       6      BBBBB",
            "ATOM 7 DDDDDDDDDDDDDDD", 
            "ATOM 8 DDDDDDDDDDDDDDD", 
            "ENDMDL",
            "MODEL        2",
            "ATOM 1 AAAAAAAAAAAAAAA", 
            "ATOM 2 AAAAAAAAAAAAAAA", 
            "TER       3      AAAAA",
            "ATOM 4 BBBBBBBBBBBBBBB", 
            "ATOM 5 BBBBBBBBBBBBBBB", 
            "TER       6      BBBBB",
            "ATOM 7 DDDDDDDDDDDDDDD", 
            "ATOM 8 DDDDDDDDDDDDDDD", 
            "ENDMDL",
            "MODEL        3",
            "ATOM 1 AAAAAAAAAAAAAAA", 
            "ATOM 2 AAAAAAAAAAAAAAA", 
            "TER       3      AAAAA",
            "ATOM 4 BBBBBBBBBBBBBBB", 
            "ATOM 5 BBBBBBBBBBBBBBB", 
            "TER       6      BBBBB",
            "ATOM 7 DDDDDDDDDDDDDDD", 
            "ATOM 8 DDDDDDDDDDDDDDD", 
            "ENDMDL",
        ])
        self.assertEqual([c[0] for c in mock_atom.call_args_list], [
            ({"id": "1", "pdbx_PDB_model_num": "1", "label_asym_id": "A", "label_entity_id": "1"}, 1),
            ({"id": "2", "pdbx_PDB_model_num": "1", "label_asym_id": "A", "label_entity_id": "1"}, 2),
            ({"id": "3", "pdbx_PDB_model_num": "1", "label_asym_id": "B", "label_entity_id": "1"}, 4),
            ({"id": "4", "pdbx_PDB_model_num": "1", "label_asym_id": "B", "label_entity_id": "1"}, 5),
            ({"id": "7", "pdbx_PDB_model_num": "1", "label_asym_id": "D", "label_entity_id": "3"}, 7),
            ({"id": "8", "pdbx_PDB_model_num": "1", "label_asym_id": "D", "label_entity_id": "3"}, 8),
            ({"id": "1", "pdbx_PDB_model_num": "2", "label_asym_id": "A", "label_entity_id": "1"}, 1),
            ({"id": "2", "pdbx_PDB_model_num": "2", "label_asym_id": "A", "label_entity_id": "1"}, 2),
            ({"id": "3", "pdbx_PDB_model_num": "2", "label_asym_id": "B", "label_entity_id": "1"}, 4),
            ({"id": "4", "pdbx_PDB_model_num": "2", "label_asym_id": "B", "label_entity_id": "1"}, 5),
            ({"id": "7", "pdbx_PDB_model_num": "2", "label_asym_id": "D", "label_entity_id": "3"}, 7),
            ({"id": "8", "pdbx_PDB_model_num": "2", "label_asym_id": "D", "label_entity_id": "3"}, 8),
            ({"id": "1", "pdbx_PDB_model_num": "3", "label_asym_id": "A", "label_entity_id": "1"}, 1),
            ({"id": "2", "pdbx_PDB_model_num": "3", "label_asym_id": "A", "label_entity_id": "1"}, 2),
            ({"id": "3", "pdbx_PDB_model_num": "3", "label_asym_id": "B", "label_entity_id": "1"}, 4),
            ({"id": "4", "pdbx_PDB_model_num": "3", "label_asym_id": "B", "label_entity_id": "1"}, 5),
            ({"id": "7", "pdbx_PDB_model_num": "3", "label_asym_id": "D", "label_entity_id": "3"}, 7),
            ({"id": "8", "pdbx_PDB_model_num": "3", "label_asym_id": "D", "label_entity_id": "3"}, 8),
        ])
        self.assertEqual([c[0] for c in mock_aniso.call_args_list], [
            ({"id": "1", "pdbx_PDB_model_num": "1", "label_asym_id": "A", "label_entity_id": "1"}, None, 1),
            ({"id": "2", "pdbx_PDB_model_num": "1", "label_asym_id": "A", "label_entity_id": "1"}, None, 2),
            ({"id": "3", "pdbx_PDB_model_num": "1", "label_asym_id": "B", "label_entity_id": "1"}, None, 4),
            ({"id": "4", "pdbx_PDB_model_num": "1", "label_asym_id": "B", "label_entity_id": "1"}, None, 5),
            ({"id": "7", "pdbx_PDB_model_num": "1", "label_asym_id": "D", "label_entity_id": "3"}, None, 7),
            ({"id": "8", "pdbx_PDB_model_num": "1", "label_asym_id": "D", "label_entity_id": "3"}, None, 8),
            ({"id": "1", "pdbx_PDB_model_num": "2", "label_asym_id": "A", "label_entity_id": "1"}, None, 1),
            ({"id": "2", "pdbx_PDB_model_num": "2", "label_asym_id": "A", "label_entity_id": "1"}, None, 2),
            ({"id": "3", "pdbx_PDB_model_num": "2", "label_asym_id": "B", "label_entity_id": "1"}, None, 4),
            ({"id": "4", "pdbx_PDB_model_num": "2", "label_asym_id": "B", "label_entity_id": "1"}, None, 5),
            ({"id": "7", "pdbx_PDB_model_num": "2", "label_asym_id": "D", "label_entity_id": "3"}, None, 7),
            ({"id": "8", "pdbx_PDB_model_num": "2", "label_asym_id": "D", "label_entity_id": "3"}, None, 8),
            ({"id": "1", "pdbx_PDB_model_num": "3", "label_asym_id": "A", "label_entity_id": "1"}, None, 1),
            ({"id": "2", "pdbx_PDB_model_num": "3", "label_asym_id": "A", "label_entity_id": "1"}, None, 2),
            ({"id": "3", "pdbx_PDB_model_num": "3", "label_asym_id": "B", "label_entity_id": "1"}, None, 4),
            ({"id": "4", "pdbx_PDB_model_num": "3", "label_asym_id": "B", "label_entity_id": "1"}, None, 5),
            ({"id": "7", "pdbx_PDB_model_num": "3", "label_asym_id": "D", "label_entity_id": "3"}, None, 7),
            ({"id": "8", "pdbx_PDB_model_num": "3", "label_asym_id": "D", "label_entity_id": "3"}, None, 8),
        ])



class AtomLineCreation(TestCase):

    def test_atom_full_values(self):
        self.assertEqual(create_atom_line({
            "group_PDB": "ATOM", "id": "46479", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": "Q", "label_comp_id": "ALA",
            "label_asym_id": "N", "label_entity_id": "14", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "D", "Cartn_x": "-148.648", "Cartn_y": "29.717",
            "Cartn_z": "-30.517", "occupancy": "1.00", "B_iso_or_equiv": "21.86",
            "pdbx_formal_charge": "1-", "auth_seq_id": "1", "auth_comp_id": "ALA",
            "auth_asym_id": "N", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        }, 20), "ATOM     20  N   ALA N   1D   -148.648  29.717 -30.517  1.00 21.86           N-1")
    

    def test_hetatm_empty_values(self):
        self.assertEqual(create_atom_line({
            "group_PDB": "HETATM", "id": "51464", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "JA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-124.828", "Cartn_y": "42.224",
            "Cartn_z": "13.312", "occupancy": "1.00", "B_iso_or_equiv": "128.95",
            "pdbx_formal_charge": "?", "auth_seq_id": "1548", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        }, 20), "HETATM   20  MG  MG  A1548    -124.828  42.224  13.312  1.00128.95          MG  ")



class AnisouLineCreation(TestCase):

    def test_can_handle_no_ainso(self):
        self.assertEqual(create_aniso_line({
            "group_PDB": "ATOM", "id": "46479", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": "Q", "label_comp_id": "ALA",
            "label_asym_id": "N", "label_entity_id": "14", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "D", "Cartn_x": "-148.648", "Cartn_y": "29.717",
            "Cartn_z": "-30.517", "occupancy": "1.00", "B_iso_or_equiv": "21.86",
            "pdbx_formal_charge": "1-", "auth_seq_id": "1", "auth_comp_id": "ALA",
            "auth_asym_id": "N", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        }, None, 20), None)


    def test_atom_full_values(self):
        self.assertEqual(create_aniso_line({
            "group_PDB": "ATOM", "id": "46479", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": "Q", "label_comp_id": "ALA",
            "label_asym_id": "N", "label_entity_id": "14", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "D", "Cartn_x": "-148.648", "Cartn_y": "29.717",
            "Cartn_z": "-30.517", "occupancy": "1.00", "B_iso_or_equiv": "21.86",
            "pdbx_formal_charge": "1-", "auth_seq_id": "1", "auth_comp_id": "ALA",
            "auth_asym_id": "N", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        }, {
            "id": "1319", "type_symbol": "C", "pdbx_label_atom_id": "C4",
            "pdbx_label_alt_id": ".", "pdbx_label_comp_id": "DG",
            "pdbx_label_asym_id": "C", "pdbx_label_seq_id": "16",
            "pdbx_PDB_ins_code": "?", "U[1][1]": "0.4728", "U[2][2]": "0.6687",
            "U[3][3]": "0.5867", "U[1][2]": "0.2574", "U[1][3]": "-0.0037",
            "U[2][3]": "-0.1117", "pdbx_auth_seq_id": "15", "pdbx_auth_comp_id": "DG",
            "pdbx_auth_asym_id": "B", "pdbx_auth_atom_id ": "C4"
        }, 20), "ANISOU   20  N   ALA N   1D    4728   6687   5867   2574    -37  -1117       N-1")



class PdbDateTests(TestCase):

    def test_can_handle_empty_date(self):
        self.assertEqual(create_pdb_date(""), None)
        self.assertEqual(create_pdb_date("   "), None)
    

    def test_can_get_date(self):
        self.assertEqual(create_pdb_date("2015-01-01"), "01-JAN-15")
        self.assertEqual(create_pdb_date("2001-11-28"), "28-NOV-01")



class SplitLinesTests(TestCase):

    def test_can_split_lines_under_length(self):
        self.assertEqual(split_lines("xxxxxxxx", 9), ["xxxxxxxx"])
    

    def test_can_split_lines_over_length_with_spaces(self):
        self.assertEqual(split_lines("abc def ghi jkl", 8), ["abc def", "ghi jkl"])
    

    def test_can_handle_no_spaces_within_limit(self):
        self.assertEqual(split_lines("abcdefghi jkl", 8), ["abcdefgh", "i jkl"])
    

    def test_can_split_lines_over_length_without_spaces(self):
        self.assertEqual(split_lines("abc,def,ghi,jkl", 8), ["abc,def,", "ghi,jkl"])



class MmcifNamesToPdbNamesTests(TestCase):

    def test_can_convert_names(self):
        self.assertEqual(mmcif_names_to_pdb_names(["Wu, N.", "Pai, E.F."]), "N.WU,E.F.PAI")