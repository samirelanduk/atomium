from unittest import TestCase
from molecupy.pdbfile import PdbFile

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

    def test_can_create_file(self):
        pdb_file = PdbFile(self.file_string)
        self.assertEqual(pdb_file._file_string, self.file_string)
        self.assertEqual(len(pdb_file._records), 9)
        self.assertEqual(pdb_file._records[0].name(), "HEADER")
        self.assertEqual(pdb_file._records[-1].name(), "COMPND")
        self.assertTrue(pdb_file._records[0].content().startswith("    LYASE"))
        self.assertTrue(pdb_file._records[-1].content().startswith("   6 ENGINEERED"))
