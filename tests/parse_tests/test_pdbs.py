from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.parse.pdb import Pdb
from atomium.structures.models import Model

class PdbCreationTests(TestCase):

    def test_can_create_pdb(self):
        pdb = Pdb()
        self.assertIsNone(pdb._model)



class PdbReprTests(TestCase):

    def test_pdb_repr(self):
        pdb = Pdb()
        self.assertEqual(str(pdb), "<Pdb>")



class PdbModelTests(TestCase):

    def test_model_property(self):
        pdb = Pdb()
        pdb._model = "snarglefargle"
        self.assertIs(pdb._model, pdb.model())


    def test_can_change_model(self):
        model = Mock(Model)
        pdb = Pdb()
        pdb._model = "snarglefargle"
        pdb.model(model)
        self.assertIs(pdb._model, model)


    def test_pdb_model_must_be_model(self):
        pdb = Pdb()
        pdb._model = "snarglefargle"
        with self.assertRaises(TypeError):
            pdb.model(100)
