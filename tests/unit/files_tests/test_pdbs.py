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
        self.assertEqual(pdb.models(), ("1", "2", "3"))



class PdbModelTests(TestCase):

    def test_model_gets_first_model(self):
        pdb = Pdb()
        pdb._models = ["1", "2", "3"]
        self.assertEqual(pdb.model(), "1")


    def test_can_get_no_model(self):
        pdb = Pdb()
        self.assertIsNone(pdb.model())



class PdbCodeTests(TestCase):

    def test_can_get_pdb_code(self):
        pdb = Pdb()
        pdb._code = "1xxx"
        self.assertIs(pdb._code, pdb.code())


    def test_can_update_code(self):
        pdb = Pdb()
        pdb._code = "1xxx"
        pdb.code("2yyy")
        self.assertEqual(pdb._code, "2yyy")


    def test_code_must_be_str(self):
        pdb = Pdb()
        with self.assertRaises(TypeError):
            pdb.code(100)


    def test_code_must_be_valid(self):
        pdb = Pdb()
        with self.assertRaises(ValueError):
            pdb.code("1xxxx")
        with self.assertRaises(ValueError):
            pdb.code("1xx")



class PdbDateTests(TestCase):

    def test_can_get_pdb_date(self):
        pdb = Pdb()
        pdb._deposition_date = "date"
        self.assertIs(pdb._deposition_date, pdb.deposition_date())


    def test_can_update_date(self):
        pdb = Pdb()
        pdb._deposition_date = "date"
        pdb.deposition_date(datetime(2017, 9, 21).date())
        self.assertEqual(pdb._deposition_date, datetime(2017, 9, 21).date())


    def test_date_must_be_date(self):
        pdb = Pdb()
        with self.assertRaises(TypeError):
            pdb.deposition_date("date")



class PdbTitleTests(TestCase):

    def test_can_get_pdb_title(self):
        pdb = Pdb()
        pdb._title = "TTT"
        self.assertIs(pdb._title, pdb.title())


    def test_can_update_title(self):
        pdb = Pdb()
        pdb._title = "TTT"
        pdb.title("TTTTTTT")
        self.assertEqual(pdb._title, "TTTTTTT")


    def test_title_must_be_str(self):
        pdb = Pdb()
        with self.assertRaises(TypeError):
            pdb.title(100)



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
