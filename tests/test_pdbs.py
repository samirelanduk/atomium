import datetime
from unittest import TestCase
from molecupy.pdbfile import PdbFile
from molecupy.pdbdatafile import PdbDataFile
from molecupy.pdb import Pdb
from molecupy.structures import PdbModel
from molecupy.exceptions import *
from molecupy import pdb_from_string, get_pdb_remotely, get_pdb_from_file

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



class ModelTests(PdbTest):

    def setUp(self):
        PdbTest.setUp(self)
        self.multi_model = Pdb(PdbDataFile(PdbFile(
         "MODEL        1\n"
         "ATOM    107  N   GLY A  13      12.681  37.302 -25.211 1.000 15.56           N\n"
         "ENDMDL\n"
         "MODEL        2\n"
         "ATOM    107  N   GLY A  13      12.681  37.302 -25.211 1.000 15.56           N\n"
         "ENDMDL"
        )))
        self.single_model = Pdb(PdbDataFile(PdbFile(
         "ATOM    107  N   GLY A  13      12.681  37.302 -25.211 1.000 15.56           N"
        )))


    def test_single_model_case(self):
        self.assertEqual(len(self.single_model.models), 1)
        self.assertEqual(len(self.empty.models), 1)
        self.assertIsInstance(self.single_model.models[0], PdbModel)
        self.assertIsInstance(self.single_model.models[0], PdbModel)


    def test_multi_models(self):
        self.assertEqual(len(self.multi_model.models), 2)
        self.assertIsInstance(self.multi_model.models[0], PdbModel)
        self.assertIsInstance(self.multi_model.models[1], PdbModel)
        self.assertIsNot(self.multi_model.models[0], self.multi_model.models[1])


    def test_model(self):
        self.assertIs(self.empty.model, self.empty.models[0])
        self.assertIs(self.single_model.model, self.single_model.models[0])
        self.assertIs(self.multi_model.model, self.multi_model.models[0])



class PdbFromStringTests(TestCase):

    def test_can_get_pdb_from_string(self):
        pdb = pdb_from_string("TITLE     CRYSTAL")
        self.assertIsInstance(pdb, Pdb)
        self.assertEqual(pdb.title, "CRYSTAL")



class PdbFromFileTests(TestCase):

    def test_can_get_pdb_from_file(self):
        pdb = get_pdb_from_file("tests/1AOI.pdb")
        self.assertIsInstance(pdb, Pdb)
        self.assertEqual(pdb.classification, "DNA BINDING PROTEIN/DNA")



class PdbFromRemotelyTests(TestCase):

    def test_can_get_pdb_remotely(self):
        pdb = get_pdb_remotely("1NVQ")
        self.assertIsInstance(pdb, Pdb)
        self.assertEqual(pdb.pdb_code, "1NVQ")


    def test_invalid_code_raises_error(self):
        with self.assertRaises(InvalidPdbCodeError):
            pdb = get_pdb_remotely("XXXX")
