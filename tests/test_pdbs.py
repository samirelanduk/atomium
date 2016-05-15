import datetime
from unittest import TestCase
from molecupy.pdbfile import PdbFile
from molecupy.pdbdatafile import PdbDataFile
from molecupy.pdb import Pdb

class PdbTest(TestCase):

    def setUp(self):
        self.empty = Pdb(PdbDataFile(PdbFile("")))



class PdbPropertiesTests(PdbTest):

    def test_has_pdb_file(self):
        self.assertIsInstance(self.empty.data_file, PdbDataFile)


    def test_repr(self):
        self.assertRegex(
         str(self.empty),
         r"<Pdb \(([^\s]{4})\)>"
        )


    def test_other_properties(self):
        self.assertEqual(
         self.empty.classification,
         self.empty.data_file.classification
        )
        self.assertEqual(
         self.empty.deposition_date,
         self.empty.data_file.deposition_date
        )
        self.assertEqual(
         self.empty.pdb_code,
         self.empty.data_file.pdb_code
        )
        self.assertEqual(
         self.empty.is_obsolete,
         self.empty.data_file.is_obsolete
        )
        self.assertEqual(
         self.empty.obsolete_date,
         self.empty.data_file.obsolete_date
        )
        self.assertEqual(
         self.empty.replacement_code,
         self.empty.data_file.replacement_code
        )
        self.assertEqual(
         self.empty.title,
         self.empty.data_file.title
        )
        self.assertEqual(
         self.empty.split_codes,
         self.empty.data_file.split_codes
        )
        self.assertEqual(
         self.empty.caveat,
         self.empty.data_file.caveat
        )
        self.assertEqual(
         self.empty.keywords,
         self.empty.data_file.keywords
        )
        self.assertEqual(
         self.empty.experimental_techniques,
         self.empty.data_file.experimental_techniques
        )
        self.assertEqual(
         self.empty.model_count,
         self.empty.data_file.model_count
        )
        self.assertEqual(
         self.empty.model_annotations,
         self.empty.data_file.model_annotations
        )
        self.assertEqual(
         self.empty.revisions,
         self.empty.data_file.revisions
        )
        self.assertEqual(
         self.empty.supercedes,
         self.empty.data_file.supercedes
        )
        self.assertEqual(
         self.empty.supercede_date,
         self.empty.data_file.supercede_date
        )
        self.assertEqual(
         self.empty.journal,
         self.empty.data_file.journal
        )
