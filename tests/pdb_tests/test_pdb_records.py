from unittest import TestCase
from atomium.pdb.pdbfile import PdbRecord

class PdbRecordCreationTests(TestCase):

    def test_can_create_pdb_record(self):
        record = PdbRecord("RECORD XXX YYY ZZZ 01")
        self.assertEqual(record._text, "RECORD XXX YYY ZZZ 01")
