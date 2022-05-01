from unittest import TestCase
import atomium

class DictParsingTests(TestCase):

    def test_1lol_pdb(self):
        d = atomium.open("tests/integration/files/1lol.pdb", dictionary=True)
        
        self.assertEqual(d["entry"], [{"id": "1LOL"}])
        self.assertEqual(d["struct_keywords"], [
            {"entry_id": "1LOL", "pdbx_keywords": "LYASE", "text": "TIM BARREL, LYASE"}
        ])
        self.assertEqual(d["pdbx_database_status"], [{
            "status_code": "REL", "entry_id": "1LOL",
            "recvd_initial_deposition_date": "2002-05-06"
        }])

        self.assertEqual(d["struct"], [{
            "entry_id": "1LOL", "pdbx_CASP_flag": "?", "pdbx_descriptor": "?",
            "pdbx_model_details": "?", "pdbx_model_type_details": "?",
            "title": "CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE COMPLEX WITH XMP"
        }])

        self.assertEqual(d["exptl"], [{"entry_id": "1LOL", "method": "X-RAY DIFFRACTION"}])

        self.assertEqual(d["audit_author"], [
            {"name": "Wu, N", "pdbx_ordinal": "1"},
            {"name": "Pai, E.F", "pdbx_ordinal": "2"},
        ])

        self.assertEqual(d["pdbx_audit_revision_history"], [{
            "ordinal": "1", "data_content_type": "Structure model",
            "major_revision": "1", "minor_revision": "0", "revision_date": "2002-08-07"
        }, {
            "ordinal": "2", "data_content_type": "Structure model", 
            "major_revision": "2","minor_revision": "0", "revision_date": "2002-08-14"
        }, {
            "ordinal": "3", "data_content_type": "Structure model",
            "major_revision": "3", "minor_revision": "0", "revision_date": "2003-04-01"
        }, {
            "ordinal": "4", "data_content_type": "Structure model",
            "major_revision": "4", "minor_revision": "0", "revision_date": "2009-02-24"
        }])

        self.assertEqual(d["citation"], [{
            "id": "primary",
            "title": "CRYSTAL STRUCTURES OF INHIBITOR COMPLEXES REVEAL AN ALTERNATE BINDING MODE IN OROTIDINE-5'-MONOPHOSPHATE DECARBOXYLASE.",
            "journal_abbrev": "J.Biol.Chem.", "journal_volume": "277",
            "page_first": "28080", "page_last": "?", "year": "2002", "journal_id_ASTM": "?",
            "country": "?", "journal_id_ISSN": "0021-9258", "journal_id_CSD": "?",
            "book_publisher": "?", "pdbx_database_id_PubMed": "12011084",
            "pdbx_database_id_DOI": "10.1074/jbc.m202362200"
        }])

        self.assertEqual(d["citation_author"], [
            {"citation_id": "primary", "name": "Wu, N", "pdbx_ordinal": "1"},
            {"citation_id": "primary", "name": "Pai, E.F", "pdbx_ordinal": "2"}
        ])
        self.assertEqual(d["reflns"], [{
            "B_iso_Wilson_estimate": "?", "entry_id": "1LOL", "data_reduction_details": "?",
            "data_reduction_method": "?", "d_resolution_high": "1.90", "d_resolution_low": "27.07",
            "details": "?", "limit_h_max": "?", "limit_h_min": "?", "limit_k_max": "?",
            "limit_k_min": "?", "limit_l_max": "?", "limit_l_min": "?", "number_all": "?",
            "number_obs": "?", "observed_criterion": "?", "observed_criterion_F_max": "?",
            "observed_criterion_F_min": "?", "observed_criterion_I_max": "?",
            "observed_criterion_I_min": "?", "observed_criterion_sigma_F": "?",
            "observed_criterion_sigma_I": "?", "percent_possible_obs": "?", "R_free_details": "?",
            "Rmerge_F_all": "?", "Rmerge_F_obs": "?", "Friedel_coverage": "?", "number_gt": "?",
            "threshold_expression": "?", "pdbx_redundancy": "?", "pdbx_Rmerge_I_obs": "?",
            "pdbx_Rmerge_I_all": "?", "pdbx_Rsym_value": "?", "pdbx_netI_over_av_sigmaI": "?",
            "pdbx_netI_over_sigmaI": "?", "pdbx_res_netI_over_av_sigmaI_2": "?",
            "pdbx_res_netI_over_sigmaI_2": "?", "pdbx_chi_squared": "?", "pdbx_scaling_rejects": "?",
            "pdbx_d_res_high_opt": "?", "pdbx_d_res_low_opt": "?", "pdbx_d_res_opt_method": "?",
            "phase_calculation_details": "?", "pdbx_Rrim_I_all": "?", "pdbx_Rpim_I_all": "?",
            "pdbx_d_opt": "?", "pdbx_number_measured_all": "?", "pdbx_diffrn_id": "?",
            "pdbx_ordinal": "?", "pdbx_CC_half": "?", "pdbx_R_split": "?",
            "ls_R_factor_R_work": "0.193", "ls_R_factor_R_free": "0.229",
            "ls_percent_reflns_R_free": "4.900", "ls_number_reflns_R_free": "1583"
        }])

        self.assertEqual(d["diffrn"], [
            {"id": "1", "crystal_id": "1", "ambient_temp": "100", "ambient_temp_details": "?"}
        ])
        self.assertEqual(d["diffrn_detector"], [{
            "details": "?", "detector": "CCD", "diffrn_id": "1",
            "pdbx_collection_date": "2001-09-12", "type": "ADSC QUANTUM 4"
        }])
    

    def test_1CK8_pdb(self):
        # Tests OBSLTE, CAVEAT, AUTHOR
        d = atomium.open("tests/integration/files/1ck8.pdb", dictionary=True)
        self.assertEqual(d["pdbx_database_PDB_obs_spr"], [{
            "id": "OBSLTE", "date": "2006-07-25", "pdb_id": "1T0K",
            "replace_pdb_id": "1CK8", "details": "?"
        }])
        self.assertEqual(d["database_PDB_caveat"], [{
            "text": "THERE ARE CHIRALITY ERRORS IN C-ALPHA CENTERS"
        }])
        self.assertEqual(d["audit_author"], [
            {"name": "Mao, H", "pdbx_ordinal": "1"},
            {"name": "White, S.A", "pdbx_ordinal": "2"},
            {"name": "Willamson, J.R", "pdbx_ordinal": "3"},
        ])
    

    def test_4HHB_pdb(self):
        # Tests SPRSDE, CAVEAT, REVDAT
        d = atomium.open("tests/integration/files/4hhb.pdb", dictionary=True)
        self.assertEqual(d["pdbx_database_PDB_obs_spr"], [{
            "id": "SPRSDE", "date": "1984-07-17", "pdb_id": "4HHB",
            "replace_pdb_id": "1HHB", "details": "?"
        }])
        self.assertEqual(d["database_PDB_caveat"], [{
            "text": "THR A 137 HAS WRONG CHIRALITY AT ATOM CB "
            "THR B 12 HAS WRONG CHIRALITY AT ATOM CB "
            "THR B 50 HAS WRONG CHIRALITY AT ATOM CB "
            "ASN C 78 HAS WRONG CHIRALITY AT ATOM CA "
            "THR C 118 HAS WRONG CHIRALITY AT ATOM CB "
            "HIS D 2 HAS WRONG CHIRALITY AT ATOM CA "
            "SER D 72 HAS WRONG CHIRALITY AT ATOM CA "
            "ASP D 73 HAS WRONG CHIRALITY AT ATOM CA "
            "LEU D 78 HAS WRONG CHIRALITY AT ATOM CA "
            "LYS D 144 HAS WRONG CHIRALITY AT ATOM CA"
        }])
        self.assertEqual(d["pdbx_audit_revision_history"], [{
           "ordinal": "1", "data_content_type": "Structure model",
           "major_revision": "1", "minor_revision": "0", "revision_date": "1984-07-17"
        }, {
            "ordinal": "2", "data_content_type": "Structure model",
            "major_revision": "2", "minor_revision": "0", "revision_date": "1989-10-15"
        }, {
            "ordinal": "3", "data_content_type": "Structure model",
            "major_revision": "3", "minor_revision": "0", "revision_date": "2003-04-01"
        }, {
            "ordinal": "4", "data_content_type": "Structure model",
            "major_revision": "4", "minor_revision": "0", "revision_date": "2009-02-24"
        }, {
            "ordinal": "5", "data_content_type": "Structure model",
            "major_revision": "5", "minor_revision": "0", "revision_date": "2011-07-13"
        }, {
            "ordinal": "6", "data_content_type": "Structure model",
            "major_revision": "6", "minor_revision": "0", "revision_date": "2020-06-17"
        }, {
            "ordinal": "7", "data_content_type": "Structure model",
            "major_revision": "7", "minor_revision": "0", "revision_date": "2021-03-31"
        }])
    

    def test_1GDJ_pdb(self):
        # Tests SPRSDE
        d = atomium.open("tests/integration/files/1gdj.pdb", dictionary=True)
        self.assertEqual(d["pdbx_database_PDB_obs_spr"], [{
            "id": "SPRSDE", "date": "1995-02-27", "pdb_id": "1GDJ",
            "replace_pdb_id": "1LH4 2LH4", "details": "?"
        }])
    

    def test_4CWU_pdb(self):
        # Tests OBSLTE, SPRSDE
        d = atomium.open("tests/integration/files/4cwu.pdb", dictionary=True)
        self.assertEqual(d["pdbx_database_PDB_obs_spr"], [{
            "id": "SPRSDE", "date": "2014-08-06", "pdb_id": "4CWU",
            "replace_pdb_id": "1VSZ", "details": "?"
        }, {
            "id": "OBSLTE", "date": "2018-05-02", "pdb_id": "6CGV",
            "replace_pdb_id": "4CWU", "details": "?"
        }])
    

    def test_1VUA_pdb(self):
        # Tests SPLIT
        d = atomium.open("tests/integration/files/1vua.pdb", dictionary=True)
        self.assertEqual(d["pdbx_database_related"], [{
            "db_name": "PDB", "db_id": code, "content_type": "split",
            "details": f"Split {n}"
        } for n, code in enumerate([
            "1VU4", "1VU5", "1VU6", "1VU7", "1VU8", "1VU9", "1VUA", "1VUC",
            "1VUD", "1VUE", "1VUF", "1VUG", "1VUH", "1VUI", "1VUJ", "1VUK",
            "1VUL", "1VUM", "1VUN", "1VUO", "1VUP", "1VUQ", "1VUR", "1VUS", "1VUT"
        ], start=1)])
    

    def test_2CSE_pdb(self):
        # Tests MDLTYP
        d = atomium.open("tests/integration/files/2cse.pdb", dictionary=True)
        self.assertEqual(d["struct"], [{
            "entry_id": "2CSE", "pdbx_CASP_flag": "?", "pdbx_descriptor": "?",
            "pdbx_model_details": "?",
            "pdbx_model_type_details": "?",
            "title": "FEATURES OF REOVIRUS OUTER-CAPSID PROTEIN MU1 REVEALED BY "
            "ELECTRON AND IMAGE RECONSTRUCTION OF THE VIRION AT 7.0-A RESOLUTION"
        }])
        self.assertEqual(d["pdbx_coordinate_model"], [{
            "asym_id": char, "type": "CA ATOMS ONLY"
        } for char in "ABCSDEFMNOGHIPQRJKLTUVWXYZ1"
        ])
    

    def test_1GIY_pdb(self):
        # Tests MDLTYP
        d = atomium.open("tests/integration/files/1giy.pdb", dictionary=True)
        self.assertEqual(d["struct"][0]["pdbx_model_type_details"], "?")
        self.assertEqual(d["pdbx_coordinate_model"], [{
            "asym_id": char, "type": "CA ATOMS ONLY"
        } for char in "CDEFGHIJKLMNOPQRSTUVWX"] + [{
            "asym_id": char, "type": "P ATOMS ONLY"
        } for char in "AB"])
    

    def test_2JRQ_pdb(self):
        # Tests MDLTYP
        d = atomium.open("tests/integration/files/2jrq.pdb", dictionary=True)
        self.assertEqual(d["struct"], [{
            "entry_id": "2JRQ", "pdbx_CASP_flag": "?", "pdbx_descriptor": "?",
            "pdbx_model_details": "?",
            "pdbx_model_type_details": "MINIMIZED AVERAGE",
            "title": "NMR SOLUTION STRUCTURE OF THE ANTICODON OF E. COLI TRNA-VAL3 WITH 1 MODIFICATION (CMO5U34)"
        }])
        self.assertNotIn("pdbx_coordinate_model", d)
    

    def test_1COJ_pdb(self):
        # Tests JRNL
        d = atomium.open("tests/integration/files/1coj.pdb", dictionary=True)
        self.assertEqual(d["citation"], [{
            "id": "primary",
            "title": "THE CRYSTAL STRUCTURE OF AN FE-SUPEROXIDE DISMUTASE FROM THE HYPERTHERMOPHILE AQUIFEX PYROPHILUS AT 1.9 A RESOLUTION: STRUCTURAL BASIS FOR THERMOSTABILITY.",
            "journal_abbrev": "J.Mol.Biol.", "journal_volume": "270",
            "page_first": "259", "page_last": "?", "year": "1997",
            "journal_id_ASTM": "?", "country": "?", "journal_id_ISSN": "0022-2836",
            "journal_id_CSD": "?", "book_publisher": "U.K. : ACADEMIC PRESS",
            "pdbx_database_id_PubMed": "9236127",
            "pdbx_database_id_DOI": "10.1006/jmbi.1997.1105"
        }])
        self.assertEqual(d["citation_author"], [
            {"citation_id": "primary", "name": "Lim, J.H", "pdbx_ordinal": "1"},
            {"citation_id": "primary", "name": "Yu, Y.G", "pdbx_ordinal": "2"},
            {"citation_id": "primary", "name": "Han, Y.S", "pdbx_ordinal": "3"},
            {"citation_id": "primary", "name": "Cho, S", "pdbx_ordinal": "4"},
            {"citation_id": "primary", "name": "Ahn, B.Y", "pdbx_ordinal": "5"},
            {"citation_id": "primary", "name": "Kim, S.H", "pdbx_ordinal": "6"},
            {"citation_id": "primary", "name": "Cho, Y", "pdbx_ordinal": "7"}
        ])
        self.assertEqual(d["citation_editor"], [
            {"citation_id": "primary", "name": "Peter Wright", "pdbx_ordinal": "1"},
        ])