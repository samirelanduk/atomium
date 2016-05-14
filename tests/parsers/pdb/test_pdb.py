import datetime
from unittest import TestCase
from molecupy.parsers.pdb.pdb_file import PdbFile
from molecupy.parsers.pdb.pdb_data_file import PdbDataFile
from molecupy.parsers.pdb.pdb import Pdb
from molecupy.macromolecules import MacroModel

class PdbTest(TestCase):

    def setUp(self):
        with open("tests/parsers/pdb/test.pdb") as f:
           pdb_file = PdbFile(f.read())
           data_file = PdbDataFile(pdb_file)
           self.pdb = Pdb(data_file)
           self.empty_pdb = Pdb(PdbDataFile(PdbFile("")))


    def check_valid_pdb(self, pdb):
        self.assertIsInstance(pdb, Pdb)
        self.assertIsInstance(pdb.data_file, PdbDataFile)
        self.assertRegex(
         str(pdb),
         r"<Pdb \(([^\s]{4})\)>"
        )


    def check_property(self, property_, type_):
        if property_ is not None:
            self.assertIsInstance(property_, type_)



class PdbCreationTests(PdbTest):

    def test_can_make_pdb(self):
        self.check_valid_pdb(self.pdb)
        self.check_valid_pdb(self.empty_pdb)



class PdbPropertiesTest(PdbTest):

    def test_pdb_has_properties(self):
        self.check_property(self.pdb.classification, str)
        self.check_property(self.pdb.deposition_date, datetime.date)
        self.check_property(self.pdb.pdb_code, str)
        self.check_property(self.pdb.is_obsolete, bool)
        self.check_property(self.pdb.obsolete_date, datetime.date)
        self.check_property(self.pdb.replacement_code, str)
        self.check_property(self.pdb.title, str)
        self.check_property(self.pdb.split_codes, list)
        self.check_property(self.pdb.caveat, str)
        self.check_property(self.pdb.keywords, list)
        self.check_property(self.pdb.experimental_techniques, list)
        self.check_property(self.pdb.model_num, int)
        self.check_property(self.pdb.model_annotations, list)
        self.check_property(self.pdb.revisions, list)
        self.check_property(self.pdb.supercedes, list)
        self.check_property(self.pdb.supercede_date, datetime.date)
        self.check_property(self.pdb.journal, dict)



class PdbModelTests(PdbTest):

    def setUp(self):
        data = """
MODEL        1
ATOM      1  N   PRO A   1       4.524   9.887  -0.667  1.00  0.00           N
ENDMDL
MODEL        2
ATOM      1  N   PRO A   1       8.041  10.551  -0.115  1.00  0.00           N
ENDMDL
MODEL        2
ATOM      1  N   PRO A   1       8.041  10.551  -0.115  1.00  0.00           N
ENDMDL"""
        self.multi_model = Pdb(PdbDataFile(PdbFile(data)))
        data = """
        ATOM      1  N   PRO A   1       4.524   9.887  -0.667  1.00  0.00           N"""
        self.single_model = Pdb(PdbDataFile(PdbFile(data)))


    def test_pdbs_have_models(self):
        self.assertIsInstance(self.multi_model.models, list)
        self.assertIsInstance(self.single_model.models, list)
        self.assertEqual(len(self.multi_model.models), 3)
        self.assertEqual(len(self.single_model.models), 1)
        for model in self.multi_model.models:
            self.assertIsInstance(model, MacroModel)
        for model in self.single_model.models:
            self.assertIsInstance(model, MacroModel)


    def test_pdbs_have_model(self):
        for pdb in (self.multi_model, self.single_model):
            self.assertIs(pdb.model, pdb.models[0])
