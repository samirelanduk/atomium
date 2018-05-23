from datetime import datetime
from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.files.pdb import Pdb
from atomium.structures.models import Model

class PdbCreationTests(TestCase):

    def test_can_create_pdb(self):
        pdb = Pdb()
        self.assertEqual(pdb._models, [])
        self.assertEqual(pdb._code, None)
        self.assertEqual(pdb._deposition_date, None)
        self.assertEqual(pdb._title, None)
        self.assertEqual(pdb._resolution, None)
        self.assertEqual(pdb._organism, None)
        self.assertEqual(pdb._expression_system, None)
        self.assertEqual(pdb._technique, None)
        self.assertEqual(pdb._classification, None)
        self.assertEqual(pdb._rfactor, None)
        self.assertEqual(pdb._keywords, [])



class PdbReprTests(TestCase):

    def test_pdb_repr_no_models(self):
        pdb = Pdb()
        self.assertEqual(str(pdb), "<Pdb (0 models)>")


    def test_pdb_repr_one_model(self):
        pdb = Pdb()
        pdb._models = ["1"]
        self.assertEqual(str(pdb), "<Pdb (1 model)>")


    def test_pdb_repr_multiple_models(self):
        pdb = Pdb()
        pdb._models = ["1", "2", "3"]
        self.assertEqual(str(pdb), "<Pdb (3 models)>")


    def test_pdb_repr_with_code(self):
        pdb = Pdb()
        pdb._code = "1XXX"
        pdb._models = ["1", "2", "3"]
        self.assertEqual(str(pdb), "<Pdb 1XXX (3 models)>")



class PdbModelsTests(TestCase):

    def test_can_get_pdb_models(self):
        pdb = Pdb()
        pdb._models = ["1", "2", "3"]
        self.assertEqual(pdb.models, ("1", "2", "3"))



class PdbModelTests(TestCase):

    def test_model_gets_first_model(self):
        pdb = Pdb()
        pdb._models = ["1", "2", "3"]
        self.assertEqual(pdb.model, "1")


    def test_can_get_no_model(self):
        pdb = Pdb()
        self.assertIsNone(pdb.model)



class PdbCodeTests(TestCase):

    def test_can_get_pdb_code(self):
        pdb = Pdb()
        pdb._code = "1xxx"
        self.assertIs(pdb._code, pdb.code)


    def test_can_update_code(self):
        pdb = Pdb()
        pdb._code = "1xxx"
        pdb.code = "2yyy"
        self.assertEqual(pdb._code, "2yyy")


    def test_code_must_be_str(self):
        pdb = Pdb()
        with self.assertRaises(TypeError):
            pdb.code = 100


    def test_code_must_be_valid(self):
        pdb = Pdb()
        with self.assertRaises(ValueError):
            pdb.code = "1xxxx"
        with self.assertRaises(ValueError):
            pdb.code = "1xx"



class PdbDateTests(TestCase):

    def test_can_get_pdb_date(self):
        pdb = Pdb()
        pdb._deposition_date = "date"
        self.assertIs(pdb._deposition_date, pdb.deposition_date)


    def test_can_update_date(self):
        pdb = Pdb()
        pdb._deposition_date = "date"
        pdb.deposition_date = datetime(2017, 9, 21).date()
        self.assertEqual(pdb._deposition_date, datetime(2017, 9, 21).date())


    def test_date_must_be_date(self):
        pdb = Pdb()
        with self.assertRaises(TypeError):
            pdb.deposition_date = "date"



class PdbTitleTests(TestCase):

    def test_can_get_pdb_title(self):
        pdb = Pdb()
        pdb._title = "TTT"
        self.assertIs(pdb._title, pdb.title)


    def test_can_update_title(self):
        pdb = Pdb()
        pdb._title = "TTT"
        pdb.title = "TTTTTTT"
        self.assertEqual(pdb._title, "TTTTTTT")


    def test_title_must_be_str(self):
        pdb = Pdb()
        with self.assertRaises(TypeError):
            pdb.title = 100



class PdbOrganismTests(TestCase):

    def test_can_get_pdb_organism(self):
        pdb = Pdb()
        pdb._organism = "GGG SSS"
        self.assertIs(pdb._organism, pdb.organism)


    def test_can_update_organism(self):
        pdb = Pdb()
        pdb._organism = "GGG SSS"
        pdb.organism = "GG SSSSSS"
        self.assertEqual(pdb._organism, "GG SSSSSS")


    def test_organism_must_be_str(self):
        pdb = Pdb()
        with self.assertRaises(TypeError):
            pdb.organism = 100



class PdbExpressionSystemTests(TestCase):

    def test_can_get_pdb_expression_system(self):
        pdb = Pdb()
        pdb._expression_system= "GGG SSS"
        self.assertIs(pdb._expression_system, pdb.expression_system)


    def test_can_update_expression_system(self):
        pdb = Pdb()
        pdb._expression_system = "GGG SSS"
        pdb.expression_system = "GG SSSSSS"
        self.assertEqual(pdb._expression_system, "GG SSSSSS")


    def test_expression_system_must_be_str(self):
        pdb = Pdb()
        with self.assertRaises(TypeError):
            pdb.expression_system = 100



class PdbTechniqueTests(TestCase):

    def test_can_get_pdb_technique(self):
        pdb = Pdb()
        pdb._technique = "GGG SSS"
        self.assertIs(pdb._technique, pdb.technique)


    def test_can_update_technique(self):
        pdb = Pdb()
        pdb._technique = "GGG SSS"
        pdb.technique = "GG SSSSSS"
        self.assertEqual(pdb._technique, "GG SSSSSS")


    def test_technique_must_be_str(self):
        pdb = Pdb()
        with self.assertRaises(TypeError):
            pdb.technique = 100



class PdbClassificationTests(TestCase):

    def test_can_get_pdb_classification(self):
        pdb = Pdb()
        pdb._classification = "GGG SSS"
        self.assertIs(pdb._classification, pdb.classification)


    def test_can_update_classification(self):
        pdb = Pdb()
        pdb._classification = "GGG SSS"
        pdb.classification = "GG SSSSSS"
        self.assertEqual(pdb._classification, "GG SSSSSS")


    def test_classification_must_be_str(self):
        pdb = Pdb()
        with self.assertRaises(TypeError):
            pdb.classification = 100



class PdbResolutionTests(TestCase):

    def test_can_get_pdb_resolution(self):
        pdb = Pdb()
        pdb._resolution = 1.2
        self.assertIs(pdb._resolution, pdb.resolution)


    def test_can_update_resolution(self):
        pdb = Pdb()
        pdb._resolution = 1.2
        pdb.resolution = 1.5
        self.assertEqual(pdb._resolution, 1.5)


    def test_resolution_must_be_number(self):
        pdb = Pdb()
        with self.assertRaises(TypeError):
            pdb.resolution = "100"
        pdb.resolution = 3



class PdbRfactorTests(TestCase):

    def test_can_get_pdb_rfactor(self):
        pdb = Pdb()
        pdb._rfactor = 1.2
        self.assertIs(pdb._rfactor, pdb.rfactor)


    def test_can_update_rfactor(self):
        pdb = Pdb()
        pdb._rfactor = 1.2
        pdb.rfactor = 1.5
        self.assertEqual(pdb._rfactor, 1.5)


    def test_rfactor_must_be_number(self):
        pdb = Pdb()
        with self.assertRaises(TypeError):
            pdb.rfactor = "100"
        pdb.rfactor = 3



class PdbTechniqueTests(TestCase):

    def test_can_get_pdb_keywords(self):
        pdb = Pdb()
        pdb._keywords = ["a", "b"]
        self.assertIs(pdb._keywords, pdb.keywords)



class PdbToStringTests(TestCase):

    @patch("atomium.files.pdb2pdbdict.pdb_to_pdb_dict")
    @patch("atomium.files.pdbdict2pdbstring.pdb_dict_to_pdb_string")
    def test_can_get_string_from_pdb(self, mock_string, mock_dict):
        pdb = Pdb()
        pdb_dict = Mock()
        mock_string.return_value = "filecontents"
        mock_dict.return_value = pdb_dict
        s = pdb.to_file_string()
        mock_dict.assert_called_with(pdb)
        mock_string.assert_called_with(pdb_dict)
        self.assertEqual(s, "filecontents")



class PdbToFileTests(TestCase):

    @patch("atomium.files.utilities.string_to_file")
    @patch("atomium.files.pdb.Pdb.to_file_string")
    def test_can_save_pdb_to_file(self, mock_string, mock_save):
        pdb = Pdb()
        mock_string.return_value = "filestring"
        pdb.save("test.pdb")
        mock_save.assert_called_with("filestring", "test.pdb")
