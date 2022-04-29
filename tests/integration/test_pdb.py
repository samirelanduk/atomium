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