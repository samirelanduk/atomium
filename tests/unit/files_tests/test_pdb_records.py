from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.files.pdbfile import PdbRecord

class PdbRecordCreationTests(TestCase):

    def test_can_create_pdb_record(self):
        record = PdbRecord("RECORD XXX YYY ZZZ 01")
        self.assertEqual(record._text, "RECORD XXX YYY ZZZ 01" + " " * 59)
        self.assertEqual(record._number, None)


    def test_pdb_record_needs_string(self):
        with self.assertRaises(TypeError):
            PdbRecord(100)


    def test_pdb_record_must_be_less_than_80_chars(self):
        PdbRecord("." * 80)
        with self.assertRaises(ValueError):
            PdbRecord("." * 81)



class PdbRecordReprTests(TestCase):

    @patch("atomium.files.pdbfile.PdbRecord.name")
    def test_pdb_record_repr(self, mock_name):
        mock_name.return_value = "MOCKNM"
        self.assertEqual(
         str(PdbRecord("RECORD XXX YYY ZZZ 01")),
         "<MOCKNM record>"
        )



class PdbRecordEqualityTests(TestCase):

    def test_equal_records(self):
        record1 = PdbRecord("RECORD XXX YYY ZZZ 01")
        record2 = PdbRecord("RECORD XXX YYY ZZZ 01")
        self.assertEqual(record1, record2)


    def test_unequal_records(self):
        record1 = PdbRecord("RECORD AAA YYY ZZZ 01")
        record2 = PdbRecord("RECORD XXX YYY ZZZ 01")
        self.assertNotEqual(record1, record2)


    def test_records_not_equal_other_objects(self):
        record1 = PdbRecord("RECORD AAA YYY ZZZ 01")
        self.assertNotEqual(record1, {})



class PdbRecordLengthTests(TestCase):

    def test_pdb_record_length(self):
        record = PdbRecord("RECORD XXX YYY ZZZ 01")
        self.assertEqual(len(record), len(record._text.strip()))



class PdbRecordContainerTests(TestCase):

    def test_pdb_record_container(self):
        record = PdbRecord("RECORD XXX YYY ZZZ 01")
        self.assertIn("XXX", record)
        self.assertNotIn("AAA", record)



class PdbRecordIterableTests(TestCase):

    def test_can_iterate_through_pdb_record(self):
        record = PdbRecord("RECORD XXX YYY ZZZ 01")
        chars = []
        for char in record:
            chars.append(char)
        self.assertEqual(record._text, "".join(chars))



class PdbRecordIndexingTests(TestCase):

    def test_can_get_character_at_index_of_record(self):
        record = PdbRecord("RECORDXXXYYYZZZ")
        for index in range(len(record._text)):
            if record._text[index].strip():
                self.assertEqual(record[index], record._text[index])


    def test_record_negative_index(self):
        record = PdbRecord("RECORD XXX YYY ZZZ 01")
        self.assertEqual(record[-1], None)
        self.assertEqual(record[-10:], None)
        record = PdbRecord("." * 80)
        self.assertEqual(record[-1], ".")
        self.assertEqual(record[-10:], "..........")


    def test_can_get_slice_of_record(self):
        record = PdbRecord("RECORD XXX YYY ZZZ 01")
        self.assertEqual(record[2:8], "CORD X")


    def test_record_index_will_strip(self):
        record = PdbRecord("RECORD XXX YYY ZZZ 01")
        self.assertEqual(record[6:11], "XXX")


    def test_record_index_returns_None_if_empty(self):
        record = PdbRecord("RECORD XXX  YYY ZZZ 01")
        self.assertEqual(record[6], None)
        self.assertEqual(record[10:12], None)


    def test_record_index_can_return_int(self):
        record = PdbRecord("RECORD XXX  YYY ZZZ 01")
        self.assertEqual(record[21], 1)
        self.assertEqual(record[20:22], 1)
        self.assertIsInstance(record[21], int)


    def test_record_index_can_return_number_when_number_is_zero(self):
        record = PdbRecord("RECORD XXX  YYY ZZZ 0.0")
        self.assertEqual(record[20], 0)
        self.assertIsInstance(record[20], int)
        self.assertEqual(record[20:23], 0.0)
        self.assertIsInstance(record[20:23], float)


    def test_record_index_can_return_float(self):
        record = PdbRecord("RECORD 1.23 YYY ZZZ 01")
        self.assertEqual(record[6:12], 1.23)



class PdbRecordForceStringReturnTests(TestCase):

    def test_can_force_record_to_return_string_index(self):
        record = PdbRecord("RECORD 1.23 YYY ZZZ 01")
        self.assertEqual(record.get_as_string(21), "1")


    def test_can_force_record_to_return_string_slice(self):
        record = PdbRecord("RECORD 1.23 YYY ZZZ 01")
        self.assertEqual(record.get_as_string(19, 22), " 01")



class PdbRecordNumberTests(TestCase):

    def test_record_number(self):
        record = PdbRecord("RECORD 1.23 YYY ZZZ 01")
        record._number = 10105
        self.assertEqual(record.number(), 10105)



class PdbRecordTextTests(TestCase):

    def test_can_get_pdb_record_text(self):
        record = PdbRecord("RECORD XXX YYY ZZZ 01")
        self.assertIs(record._text, record.text())


    def test_can_update_record_text(self):
        record = PdbRecord("RECORD XXX YYY ZZZ 01")
        record.text("RECORD AAA BBB CCC 02")
        self.assertEqual(record._text, "RECORD AAA BBB CCC 02" + " " * 59)


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


    def test_can_update_record_name(self):
        record = PdbRecord("RECORD XXX YYY ZZZ 01")
        record.name("DROCER")
        self.assertEqual(record._text, "DROCER XXX YYY ZZZ 01" + " " * 59)


    def test_record_name_update_will_pad(self):
        record = PdbRecord("RECORD XXX YYY ZZZ 01")
        record.name("DRO")
        self.assertEqual(record._text, "DRO    XXX YYY ZZZ 01" + " " * 59)


    def test_record_name_must_be_str(self):
        record = PdbRecord("RECORD XXX YYY ZZZ 01")
        with self.assertRaises(TypeError) as e:
            record.name(100)


    def test_record_name_cannot_be_more_than_6_chars(self):
        record = PdbRecord("RECORD XXX YYY ZZZ 01")
        with self.assertRaises(ValueError):
            record.name("DROCER.")
