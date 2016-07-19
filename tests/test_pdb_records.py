from unittest import TestCase
from molecupy.pdbfile import PdbRecord

class RecordCreationTests(TestCase):

    def test_can_create_record(self):
        record = PdbRecord("TEST   123  123.8    HYT", 23)
        self.assertEqual(record._number, 23)
        self.assertEqual(record._name, "TEST")
        self.assertTrue(record._text.startswith("TEST   123  123.8    HYT"))
        self.assertTrue(record._content.startswith(" 123  123.8    HYT"))
        self.assertEqual(len(record._text), 80)
        self.assertEqual(len(record._content), 74)


    def test_number_must_be_int(self):
        with self.assertRaises(TypeError):
            PdbRecord("TEST   123  123.8    HYT", "23")
        with self.assertRaises(TypeError):
            PdbRecord("TEST   123  123.8    HYT", 23.5)


    def test_text_must_be_string(self):
        with self.assertRaises(TypeError):
            PdbRecord(100, 23)


    def test_cannot_provide_empty_string(self):
        with self.assertRaises(ValueError):
            PdbRecord("", 23)


    def test_repr(self):
        record = PdbRecord("TEST   123  123.8    HYT", 23)
        self.assertEqual(str(record), "<PdbRecord 23 (TEST)>")
