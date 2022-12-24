from unittest import TestCase
from unittest.mock import patch
from atomium.pdb import *

class PdbStringToDictTests(TestCase):

    @patch("atomium.pdb.parse_metadata")
    @patch("atomium.pdb.get_entities")
    @patch("atomium.pdb.build_structure_categories")
    def test_can_get_dict(self, mock_build, mock_get, mock_parse):
        mock_parse.return_value = {"mmcif": 1}
        mock_get.return_value = ["polymers", "nonpolymers"]
        d = pdb_string_to_mmcif_dict("filestring\n\n")
        self.assertEqual(d, {"mmcif": 1})
        mock_parse.assert_called_with("filestring\n\n")
        mock_get.assert_called_with("filestring\n\n")
        mock_build.assert_called_with("filestring\n\n", "polymers", "nonpolymers", {"mmcif": 1})



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
    @patch("atomium.pdb.parse_scalen")
    @patch("atomium.pdb.parse_mtrixn")
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
        for mock in mocks[1:7]:
            mock.assert_called_with("polymers", {"mmcif": 1})
        for mock in mocks[7:9]:
            mock.assert_called_with("nonpolymers", {"mmcif": 1})
        mocks[9].assert_called_with("filestring", {"mmcif": 1})
        mocks[10].assert_called_with("filestring", "polymers", "nonpolymers", {"mmcif": 1})
        mocks[11].assert_called_with("filestring", {"mmcif": 1})
        mocks[12].assert_called_with({"mmcif": 1})
        mocks[13].assert_called_with({"mmcif": 1})
    


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
        self.assertEqual(mmcif, {
            "entry": [{"id": "1XXX"}], "exptl": [{"entry_id": "1XXX", "method": "?"}]
        })
    

    def test_can_parse_expdta(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        parse_expdta("EXPDTA    NEUTRON DIFFRACTION; X-RAY DIFFRACTION  ", mmcif)
        self.assertEqual(mmcif, {
            "entry": [{"id": "1XXX"}],
            "exptl": [{"entry_id": "1XXX", "method": "NEUTRON DIFFRACTION; X-RAY DIFFRACTION"}]
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


    @patch("atomium.pdb.process_names")
    def test_can_parse_author(self, mock_names):
        mock_names.return_value = ["Berry, M.B"]
        mmcif = {}
        parse_author("AUTHOR    M.B.BERRY  ", mmcif)
        self.assertEqual(mmcif, {
            "audit_author": [{"name": "Berry, M.B", "pdbx_ordinal": "1"}]
        })
        mock_names.assert_called_with(["M.B.BERRY"])
    

    @patch("atomium.pdb.process_names")
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
        self.patch1 = patch("atomium.pdb.process_names")
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



class JournalidsParsingTests(TestCase):

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
    @patch("atomium.pdb.parse_remark_200")
    @patch("atomium.pdb.parse_remark_465")
    @patch("atomium.pdb.parse_remark_470")
    @patch("atomium.pdb.parse_remark_480")
    @patch("atomium.pdb.parse_remark_800")
    def test_can_parse_remark_records(self, *mocks):
        parse_remarks("filestring", {"mmicf": 1})
        for mock in mocks:
            mock.assert_called_with("filestring", {"mmicf": 1})



class Remark2ParsingTests(TestCase):

    def test_can_handle_no_remark_2(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        parse_remark_2("", mmcif)
        self.assertEqual(mmcif, {"entry": [{"id": "1XXX"}], "reflns": [{
            "entry_id": "1XXX", "B_iso_Wilson_estimate": "?", "data_reduction_details": "?",
            "data_reduction_method": "?", "d_resolution_high": "?", "d_resolution_low": "?",
            "details": "?", "limit_h_max": "?", "limit_h_min": "?", "limit_k_max": "?", "limit_k_min": "?",
            "limit_l_max": "?", "limit_l_min": "?", "number_all": "?", "number_obs": "?",
            "observed_criterion": "?", "observed_criterion_F_max": "?",
            "observed_criterion_F_min": "?", "observed_criterion_I_max": "?",
            "observed_criterion_I_min": "?", "observed_criterion_sigma_F": "?",
            "observed_criterion_sigma_I": "?", "percent_possible_obs": "?", "R_free_details": "?",
            "Rmerge_F_all": "?", "Rmerge_F_obs": "?", "Friedel_coverage": "?", "number_gt": "?",
            "threshold_expression": "?", "pdbx_redundancy": "?", "pdbx_Rmerge_I_obs": "?",
            "pdbx_Rmerge_I_all": "?", "pdbx_Rsym_value": "?", "pdbx_netI_over_av_sigmaI": "?",
            "pdbx_netI_over_sigmaI": "?", "pdbx_res_netI_over_av_sigmaI_2": "?",
            "pdbx_res_netI_over_sigmaI_2": "?", "pdbx_chi_squared": "?",
            "pdbx_scaling_rejects": "?", "pdbx_d_res_high_opt": "?", "pdbx_d_res_low_opt": "?",
            "pdbx_d_res_opt_method": "?", "phase_calculation_details": "?", "pdbx_Rrim_I_all": "?",
            "pdbx_Rpim_I_all": "?", "pdbx_d_opt": "?", "pdbx_number_measured_all": "?",
            "pdbx_diffrn_id": "?", "pdbx_ordinal": "?", "pdbx_CC_half": "?", "pdbx_R_split": "?"
        }]})
    

    def test_can_get_resolution(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        parse_remark_2("REMARK   2   \nREMARK   2 RESOLUTION.    1.74 ANGSTROMS.", mmcif)
        self.assertEqual(mmcif, {"entry": [{"id": "1XXX"}], "reflns": [{
            "entry_id": "1XXX", "B_iso_Wilson_estimate": "?", "data_reduction_details": "?",
            "data_reduction_method": "?", "d_resolution_high": "1.74", "d_resolution_low": "?",
            "details": "?", "limit_h_max": "?", "limit_h_min": "?", "limit_k_max": "?", "limit_k_min": "?",
            "limit_l_max": "?", "limit_l_min": "?", "number_all": "?", "number_obs": "?",
            "observed_criterion": "?", "observed_criterion_F_max": "?",
            "observed_criterion_F_min": "?", "observed_criterion_I_max": "?",
            "observed_criterion_I_min": "?", "observed_criterion_sigma_F": "?",
            "observed_criterion_sigma_I": "?", "percent_possible_obs": "?", "R_free_details": "?",
            "Rmerge_F_all": "?", "Rmerge_F_obs": "?", "Friedel_coverage": "?", "number_gt": "?",
            "threshold_expression": "?", "pdbx_redundancy": "?", "pdbx_Rmerge_I_obs": "?",
            "pdbx_Rmerge_I_all": "?", "pdbx_Rsym_value": "?", "pdbx_netI_over_av_sigmaI": "?",
            "pdbx_netI_over_sigmaI": "?", "pdbx_res_netI_over_av_sigmaI_2": "?",
            "pdbx_res_netI_over_sigmaI_2": "?", "pdbx_chi_squared": "?",
            "pdbx_scaling_rejects": "?", "pdbx_d_res_high_opt": "?", "pdbx_d_res_low_opt": "?",
            "pdbx_d_res_opt_method": "?", "phase_calculation_details": "?", "pdbx_Rrim_I_all": "?",
            "pdbx_Rpim_I_all": "?", "pdbx_d_opt": "?", "pdbx_number_measured_all": "?",
            "pdbx_diffrn_id": "?", "pdbx_ordinal": "?", "pdbx_CC_half": "?", "pdbx_R_split": "?"
        }]})
    

    def test_can_get_na_resolution(self):
        mmcif = {"entry": [{"id": "1XXX"}]}
        parse_remark_2("REMARK   2   \nREMARK   2 RESOLUTION.   NOT APPLICABLE. ", mmcif)
        self.assertEqual(mmcif, {"entry": [{"id": "1XXX"}], "reflns": [{
            "entry_id": "1XXX", "B_iso_Wilson_estimate": "?", "data_reduction_details": "?",
            "data_reduction_method": "?", "d_resolution_high": "?", "d_resolution_low": "?",
            "details": "?", "limit_h_max": "?", "limit_h_min": "?", "limit_k_max": "?", "limit_k_min": "?",
            "limit_l_max": "?", "limit_l_min": "?", "number_all": "?", "number_obs": "?",
            "observed_criterion": "?", "observed_criterion_F_max": "?",
            "observed_criterion_F_min": "?", "observed_criterion_I_max": "?",
            "observed_criterion_I_min": "?", "observed_criterion_sigma_F": "?",
            "observed_criterion_sigma_I": "?", "percent_possible_obs": "?", "R_free_details": "?",
            "Rmerge_F_all": "?", "Rmerge_F_obs": "?", "Friedel_coverage": "?", "number_gt": "?",
            "threshold_expression": "?", "pdbx_redundancy": "?", "pdbx_Rmerge_I_obs": "?",
            "pdbx_Rmerge_I_all": "?", "pdbx_Rsym_value": "?", "pdbx_netI_over_av_sigmaI": "?",
            "pdbx_netI_over_sigmaI": "?", "pdbx_res_netI_over_av_sigmaI_2": "?",
            "pdbx_res_netI_over_sigmaI_2": "?", "pdbx_chi_squared": "?",
            "pdbx_scaling_rejects": "?", "pdbx_d_res_high_opt": "?", "pdbx_d_res_low_opt": "?",
            "pdbx_d_res_opt_method": "?", "phase_calculation_details": "?", "pdbx_Rrim_I_all": "?",
            "pdbx_Rpim_I_all": "?", "pdbx_d_opt": "?", "pdbx_number_measured_all": "?",
            "pdbx_diffrn_id": "?", "pdbx_ordinal": "?", "pdbx_CC_half": "?", "pdbx_R_split": "?"
        }]})



class Remark3ParsingTests(TestCase):

    def test_can_handle_no_remark_3(self):
        mmcif = {"reflns": [{}]}
        parse_remark_3("", mmcif)
        self.assertEqual(mmcif, {"reflns": [{}]})