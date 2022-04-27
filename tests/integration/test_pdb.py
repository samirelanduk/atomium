from datetime import date
from unittest import TestCase
import atomium

class DictParsingTests(TestCase):

    def test_1lol_pdb(self):
        d = atomium.open("tests/integration/files/1lol.pdb", dictionary=True)
        
        self.assertEqual(d["entry"], [{"id": "1LOL"}])
        self.assertEqual(d["struct_keywords"], [
            {"entry_id": "1LOL", "pdbx_keywords": "LYASE", "text": "?"}
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