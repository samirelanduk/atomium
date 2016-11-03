from unittest import TestCase
import unittest.mock
from molecupy.pdb.pdbfile import PdbRecord, PdbFile

class RecordCreationTests(TestCase):

    def test_can_create_record(self):
        record = PdbRecord("TEST   123  123.8    HYT")
        self.assertEqual(record._name, "TEST")
        self.assertTrue(record._text.startswith("TEST   123  123.8    HYT"))
        self.assertTrue(record._content.startswith(" 123  123.8    HYT"))
        self.assertEqual(len(record._text), 80)
        self.assertEqual(len(record._content), 74)
        self.assertEqual(record._pdb_file, None)


    def test_text_must_be_string(self):
        with self.assertRaises(TypeError):
            PdbRecord(100)


    def test_cannot_provide_empty_string(self):
        with self.assertRaises(ValueError):
            PdbRecord("")


    def test_repr(self):
        record = PdbRecord("TEST   123  123.8    HYT")
        self.assertEqual(str(record), "<PdbRecord (TEST)>")


    def test_can_create_record_with_pdb_file(self):
        pdb_file = unittest.mock.Mock(PdbFile)
        record = PdbRecord("TEST   123  123.8    HYT", pdb_file)
        self.assertEqual(record._pdb_file, pdb_file)


    def test_pdb_file_must_be_pdbfile(self):
        with self.assertRaises(TypeError):
            record = PdbRecord("TEST   123  123.8    HYT", "pdb_file")



class RecordPropertyTests(TestCase):

    def setUp(self):
        self.record = PdbRecord("TEST   123  123.8    HYT")


    def test_record_properties(self):
        self.assertEqual(self.record._name, self.record.name())
        self.assertEqual(self.record._text, self.record.text())
        self.assertEqual(self.record._content, self.record.content())
        self.assertEqual(self.record._pdb_file, self.record.pdb_file())


    def test_can_modify_record_name(self):
        self.record.name("NOVEL")
        self.assertEqual(self.record.name(), "NOVEL")
        self.assertEqual(
         self.record.content(),
         " 123  123.8    HYT" + (" " * 56)
        )
        self.assertEqual(
         self.record.text(),
         "NOVEL  123  123.8    HYT" + (" " * 56)
        )


    def test_name_must_be_str(self):
        with self.assertRaises(TypeError):
            self.record.name(123456)


    def test_record_name_cannot_be_made_longer_than_6_chars(self):
        with self.assertRaises(ValueError):
            self.record.name("1234567")


    def test_record_number_when_associated_file_is_none(self):
        self.assertEqual(self.record.number(), None)


    def test_record_number_when_file_is_associated(self):
        pdb_file = unittest.mock.Mock(PdbFile)
        record2 = PdbRecord("TEST   129  123.8    HYT", pdb_file)
        record3 = PdbRecord("TEST   133  123.8    HYT", pdb_file)
        record4 = PdbRecord("TEST   523  123.8    HYT", pdb_file)
        pdb_file.records.return_value = [record2, record3, record4]
        self.assertEqual(record2.number(), 1)
        self.assertEqual(record3.number(), 2)
        self.assertEqual(record4.number(), 3)


    def test_can_add_pdbfile_association(self):
        self.assertEqual(self.record.pdb_file(), None)
        pdb_file = unittest.mock.Mock(PdbFile)
        self.record.pdb_file(pdb_file)
        self.assertEqual(self.record.pdb_file(), pdb_file)


    def test_can_only_add_pdbfile_as_pdb_file(self):
        with self.assertRaises(TypeError):
            self.record.pdb_file("pdb_file")




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
