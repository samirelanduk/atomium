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
