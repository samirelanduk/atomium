import random
from unittest import TestCase
from molecupy.pdbfile import PdbFile, PdbRecord

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


    def check_pdb_file(self, pdb_file):
        self.assertIsInstance(pdb_file, PdbFile)
        self.assertIsInstance(pdb_file.records, list)
        self.assertIsInstance(pdb_file.file_contents, str)
        self.assertRegex(
         str(pdb_file),
         r"<PdbFile \((\d+) records\)>"
        )



class FileCreationTests(PdbFileTest):

    def test_can_create_pdb_file(self):
        pdb_file = PdbFile(self.file_string)
        self.check_pdb_file(pdb_file)
        self.assertEqual(pdb_file.file_contents, self.file_string)


    def test_special_characters_are_stripped_out(self):
        new_file = []
        for char in self.file_string:
            new_file.append(char)
            weird_char = chr(random.choice(list(range(31)) + list(range(127, 255))))
            new_file.append(weird_char if weird_char != "\n" else "")
        new_file = "".join(new_file)
        pdb_file = PdbFile(new_file)
        self.assertEqual(pdb_file.file_contents, self.file_string)



class FileRecordTests(PdbFileTest):

    def test_files_have_records(self):
        pdb_file = PdbFile(self.file_string)
        self.assertNotEqual(pdb_file.records, [])
        for record in pdb_file.records:
            self.assertIsInstance(record, PdbRecord)


    def test_empty_lines_are_ignored(self):
        pdb_file = PdbFile(self.file_string)
        self.assertEqual(len(pdb_file.records), 9)


    def test_can_get_record_by_name(self):
        pdb_file = PdbFile(self.file_string)
        self.assertEqual(
         pdb_file.get_record_by_name("COMPND"),
         pdb_file.records[3]
        )
        self.assertEqual(
         pdb_file.get_record_by_name("XXX"),
         None
        )


    def test_can_get_multiple_records_by_name(self):
        pdb_file = PdbFile(self.file_string)
        self.assertEqual(
         pdb_file.get_records_by_name("TITLE"),
         pdb_file.records[1:3]
        )
        self.assertEqual(
         pdb_file.get_records_by_name("XXX"),
         []
        )
