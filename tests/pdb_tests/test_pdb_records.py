from unittest import TestCase
from atomium.pdb.pdbfile import PdbRecord

class PdbRecordCreationTests(TestCase):

    def test_can_create_pdb_record(self):
        record = PdbRecord("RECORD XXX YYY ZZZ 01")
        self.assertEqual(record._text, "RECORD XXX YYY ZZZ 01")


    def test_pdb_record_needs_string(self):
        with self.assertRaises(TypeError):
            PdbRecord(100)


    def test_pdb_record_must_be_less_than_80_chars(self):
        PdbRecord("." * 80)
        with self.assertRaises(ValueError):
            PdbRecord("." * 81)



class PdbRecordTextTests(TestCase):

    def test_can_get_pdb_record_text(self):
        record = PdbRecord("RECORD XXX YYY ZZZ 01")
        self.assertIs(record._text, record.text())


    def test_can_update_record_text(self):
        record = PdbRecord("RECORD XXX YYY ZZZ 01")
        record.text("RECORD AAA BBB CCC 02")
        self.assertEqual(record._text, "RECORD AAA BBB CCC 02")


    def test_record_must_be_str(self):
        record = PdbRecord("RECORD XXX YYY ZZZ 01")
        with self.assertRaises(TypeError):
            record.text(100)


    def test_record_must_be_less_than_80_chars(self):
        record = PdbRecord("RECORD XXX YYY ZZZ 01")
        with self.assertRaises(ValueError):
            record.text("." * 81)
        record.text("." * 80)



class PdbRecordNameTests(TestCase):

    def test_can_get_record_name_from_text(self):
        record = PdbRecord("RECORD XXX YYY ZZZ 01")
        self.assertEqual(record.name(), "RECORD")


    def test_record_name_will_truncate(self):
        record = PdbRecord("REC    XXX YYY ZZZ 01")
        self.assertEqual(record.name(), "REC")
