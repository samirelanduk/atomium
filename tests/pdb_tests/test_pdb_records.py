from unittest import TestCase
from molecupy.pdb.pdbfile import PdbRecord

class RecordCreationTests(TestCase):

    def test_can_create_record(self):
        record = PdbRecord("TEST   123  123.8    HYT")
        self.assertEqual(record._name, "TEST")
        self.assertTrue(record._text.startswith("TEST   123  123.8    HYT"))
        self.assertTrue(record._content.startswith(" 123  123.8    HYT"))
        self.assertEqual(len(record._text), 80)
        self.assertEqual(len(record._content), 74)


    def test_text_must_be_string(self):
        with self.assertRaises(TypeError):
            PdbRecord(100)


    def test_cannot_provide_empty_string(self):
        with self.assertRaises(ValueError):
            PdbRecord("")


    def test_repr(self):
        record = PdbRecord("TEST   123  123.8    HYT")
        self.assertEqual(str(record), "<PdbRecord (TEST)>")



class RecordPropertyTests(TestCase):

    def test_record_properties(self):
        record = PdbRecord("TEST   123  123.8    HYT")
        self.assertEqual(record._name, record.name())
        self.assertEqual(record._text, record.text())
        self.assertTrue(record._content, record.content())



class RecordAccessTests(TestCase):

    def setUp(self):
        self.record = PdbRecord("TEST   123  123.8    HYT")


    def test_can_get_individual_characters(self):
        self.assertEqual(self.record[0], "T")
        self.assertEqual(self.record[21], "H")


    def test_can_get_strings_from_record(self):
        self.assertEqual(self.record[1:4], "EST")
        self.assertEqual(self.record[21:24], "HYT")


    def test_record_indexes_will_strip_strings(self):
        self.assertEqual(self.record[0:7], "TEST")
        self.assertEqual(self.record[19:24], "HYT")
        self.assertEqual(self.record[19:34], "HYT")


    def test_records_can_covert_to_int(self):
        self.assertEqual(self.record[5:11], 123)


    def test_records_can_covert_to_float(self):
        self.assertEqual(self.record[10:19], 123.8)


    def test_can_force_record_to_return_string(self):
        self.assertEqual(self.record.get_as_string(5, 11), "123")
        self.assertEqual(self.record.get_as_string(10, 19), "123.8")


    def test_empty_sections_return_none(self):
        self.assertIs(self.record[17:21], None)
        self.assertIs(self.record.get_as_string(17, 21), None)
