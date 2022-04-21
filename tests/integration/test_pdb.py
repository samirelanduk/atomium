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