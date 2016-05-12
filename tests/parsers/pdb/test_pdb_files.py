from unittest import TestCase
from molecupy.parsers.pdb.pdb_file import PdbFile

class PdbFileTest(TestCase):

    def setUp(self):
        self.file_string = \
"""HEADER    LYASE                                   06-MAY-02   1LOL
TITLE     CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE
TITLE    2 COMPLEX WITH XMP


COMPND    MOL_ID: 1;
COMPND   2 MOLECULE: OROTIDINE 5'-MONOPHOSPHATE DECARBOXYLASE;
COMPND   3 CHAIN: A, B;
COMPND   4 SYNONYM: OMP DECARBOXYLASE, OMPDCASE, OMPDECASE;
COMPND   5 EC: 4.1.1.23;
COMPND   6 ENGINEERED: YES"""


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
