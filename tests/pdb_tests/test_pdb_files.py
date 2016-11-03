import random
from unittest import TestCase
from molecupy.pdb.pdbfile import PdbFile, PdbRecord

class PdbFileTest(TestCase):

    def setUp(self):
        self.file_string = "\n".join([
         "HEADER    LYASE                                   06-MAY-02   1LOL",
         "TITLE     CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE",
         "TITLE    2 COMPLEX WITH XMP",
         "",
         "",
         "COMPND    MOL_ID: 1;",
         "COMPND   2 MOLECULE: OROTIDINE 5'-MONOPHOSPHATE DECARBOXYLASE;",
         "COMPND   3 CHAIN: A, B;",
         "COMPND   4 SYNONYM: OMP DECARBOXYLASE, OMPDCASE, OMPDECASE;",
         "COMPND   5 EC: 4.1.1.23;",
         "COMPND   6 ENGINEERED: YES"
        ])



class PdbFileCreationTests(PdbFileTest):

    def test_can_create_file_from_string(self):
        pdb_file = PdbFile(self.file_string)
        self.assertEqual(pdb_file._file_string, self.file_string)
        self.assertEqual(len(pdb_file._records), 9)
        self.assertEqual(pdb_file._records[0].name(), "HEADER")
        self.assertEqual(pdb_file._records[-1].name(), "COMPND")
        self.assertTrue(pdb_file._records[0].content().startswith("    LYASE"))
        self.assertTrue(pdb_file._records[-1].content().startswith("   6 ENGINEERED"))


    def test_repr(self):
        pdb_file = PdbFile(self.file_string)
        self.assertEqual(str(pdb_file), "<PdbFile (9 Records)>")


    def test_special_characters_are_stripped_out(self):
        new_file = []
        for char in self.file_string:
            new_file.append(char)
            weird_char = chr(random.choice(list(range(31)) + list(range(127, 255))))
            new_file.append(weird_char if weird_char != "\n" else "")
        new_file = "".join(new_file)
        pdb_file = PdbFile(new_file)
        self.assertEqual(pdb_file._file_string, self.file_string)


    def test_can_create_file_without_string(self):
        pdb_file = PdbFile()
        self.assertEqual(len(pdb_file.records()), 0)



class PdbFilePropertiesTests(PdbFileTest):

    def test_pdb_file_properties(self):
        pdb_file = PdbFile(self.file_string)
        self.assertEqual(pdb_file.file_string(), pdb_file._file_string)
        self.assertEqual(pdb_file.records(), pdb_file._records)


    def test_can_add_records(self):
        pdb_file = PdbFile()
        self.assertEqual(len(pdb_file.records()), 0)
        record1 = PdbRecord("TEST   123  123.8    HYT", 1)
        pdb_file.add_record(record1)
        self.assertEqual(len(pdb_file.records()), 1)
        self.assertEqual(pdb_file.records()[-1].name(), "TEST")
        record2 = PdbRecord("TEST2  123  123.8    HYT", 2)
        pdb_file.add_record(record2)
        self.assertEqual(len(pdb_file.records()), 2)
        self.assertEqual(pdb_file.records()[-1].name(), "TEST2")


    def test_can_only_add_records(self):
        pdb_file = PdbFile()
        with self.assertRaises(TypeError):
            pdb_file.add_record("TEST2  123  123.8    HYT")



class PdbFileRecordTests(PdbFileTest):

    def test_can_get_record_by_name(self):
        pdb_file = PdbFile(self.file_string)
        self.assertEqual(
         pdb_file.get_record_by_name("COMPND"),
         pdb_file.records()[3]
        )
        self.assertEqual(
         pdb_file.get_record_by_name("XXX"),
         None
        )


    def test_can_only_get_record_by_string(self):
        pdb_file = PdbFile(self.file_string)
        with self.assertRaises(TypeError):
            pdb_file.get_record_by_name(100)


    def test_can_get_multiple_records_by_name(self):
        pdb_file = PdbFile(self.file_string)
        self.assertEqual(
         pdb_file.get_records_by_name("TITLE"),
         pdb_file.records()[1:3]
        )
        self.assertEqual(
         pdb_file.get_records_by_name("XXX"),
         []
        )


    def test_can_only_get_records_by_string(self):
        pdb_file = PdbFile(self.file_string)
        with self.assertRaises(TypeError):
            pdb_file.get_records_by_name(100)



class PdbFileToStringTests(PdbFileTest):

    def test_can_turn_pdb_file_back_to_string(self):
        pdb_file = PdbFile(self.file_string)
        file_lines = self.file_string.split("\n")
        file_lines = [line + (" " * (80-len(line))) for line in file_lines if line]
        target_string = "\n".join(file_lines)
        self.assertEqual(pdb_file.convert_to_string(), target_string)
