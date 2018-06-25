from datetime import datetime
from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock
from atomium.files.pdb import Pdb
from atomium.models.molecules import Model

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
        self.assertEqual(pdb._rfree, None)
        self.assertEqual(pdb._rcount, None)
        self.assertEqual(pdb._keywords, [])
        self.assertEqual(pdb._biomolecules, [])



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



class PdbFreeRfactorTests(TestCase):

    def test_can_get_pdb_rfree(self):
        pdb = Pdb()
        pdb._rfree = 1.2
        self.assertIs(pdb._rfree, pdb.rfree)


    def test_can_update_rfactor(self):
        pdb = Pdb()
        pdb._rfree = 1.2
        pdb.rfree = 1.5
        self.assertEqual(pdb._rfree, 1.5)



class PdbRfactorCountTests(TestCase):

    def test_can_get_pdb_rcount(self):
        pdb = Pdb()
        pdb._rcount = 1.2
        self.assertIs(pdb._rcount, pdb.rcount)


    def test_can_update_rcount(self):
        pdb = Pdb()
        pdb._rcount = 1.2
        pdb.rcount = 1.5
        self.assertEqual(pdb._rcount, 1.5)



class PdbKeywordsTests(TestCase):

    def test_can_get_pdb_keywords(self):
        pdb = Pdb()
        pdb._keywords = ["a", "b"]
        self.assertIs(pdb._keywords, pdb.keywords)



class PdbBiomoleculesTests(TestCase):

    def test_can_get_pdb_biomolecules(self):
        pdb = Pdb()
        pdb._biomolecules = ["a", "b"]
        self.assertIs(pdb._keywords, pdb.keywords)



class PdbAssemblyTests(TestCase):

    def setUp(self):
        self.pdb = Pdb()
        self.pdb._biomolecules = [{
         "delta_energy": -20, "id": 1, "transformations": [{
          "matrix": 1, "vector": 2, "chains": ["A", "B"]
         }, {
          "matrix": 3, "vector": 4, "chains": ["A", "C"]
         }]
        }, {
         "delta_energy": -60, "id": 2
        }, {
         "delta_energy": -7, "id": 3
        }]
        model = Mock()
        chains = [Mock(), Mock(), Mock()]
        new_chains = [Mock(), Mock(), Mock()]
        for chain, new_ in zip(chains, new_chains):
            chain.copy.return_value = new_
        model.chain.side_effect = [chains[0], chains[1], chains[0], chains[2]]
        self.pdb._models = [model, "MODEL"]
        self.chains, self.new_chains = chains, new_chains


    @patch("atomium.files.pdb.Model")
    def test_can_generate_assembly(self, mock_model):
        model = self.pdb.generate_assembly(1)
        for chain in self.chains: chain.copy.assert_any_call()
        self.new_chains[0].transform.assert_any_call(1)
        self.new_chains[0].transform.assert_any_call(3)
        self.new_chains[1].transform.assert_called_with(1)
        self.new_chains[2].transform.assert_called_with(3)
        self.new_chains[0].translate.assert_any_call(2)
        self.new_chains[0].translate.assert_any_call(4)
        self.new_chains[1].translate.assert_called_with(2)
        self.new_chains[2].translate.assert_called_with(4)
        mock_model.assert_called_with(
         self.new_chains[0], self.new_chains[1], self.new_chains[0], self.new_chains[2]
        )
        self.assertIs(model, mock_model.return_value)


    def test_id_must_be_valid(self):
        with self.assertRaises(ValueError):
            self.pdb.generate_assembly(4)



class BestAssemblyTests(TestCase):

    def test_can_get_best_assembly(self):
        pdb = Pdb()
        pdb._biomolecules = [
         {"delta_energy": -20, "id": 1,},
         {"delta_energy": -60, "id": 2},
         {"delta_energy": -7, "id": 3}
        ]
        model = pdb.best_assembly
        self.assertEqual(model, {"delta_energy": -60, "id": 2})


    def test_can_return_none_as_best_assembly(self):
        pdb = Pdb()
        model = pdb.best_assembly
        self.assertEqual(model, None)



class BestAssemblyGenerationTests(TestCase):

    @patch("atomium.files.pdb.Pdb.generate_assembly")
    @patch("atomium.files.pdb.Pdb.best_assembly", new_callable=PropertyMock)
    def test_can_generate_best_assembly(self, mock_best, mock_model):
        pdb = Pdb()
        mock_best.return_value = {"delta_energy": -60, "id": 2}
        model = pdb.generate_best_assembly()
        self.assertIs(model, mock_model.return_value)
        mock_model.assert_called_with(2)

    @patch("atomium.files.pdb.Pdb.best_assembly", new_callable=PropertyMock)
    def test_can_return_model_as_best_assembly(self, mock_best):
        mock_best.return_value = None
        pdb = Pdb()
        pdb._models = "ABCDEF"
        model = pdb.generate_best_assembly()
        self.assertEqual(model, "A")



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
