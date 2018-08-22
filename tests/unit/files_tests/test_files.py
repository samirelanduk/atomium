from datetime import date
from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.files.file import File

class FileCreationTests(TestCase):

    def test_can_create_files(self):
        file_ = File()
        self.assertIsNone(file_._filetype)
        self.assertIsNone(file_._code)
        self.assertIsNone(file_._title)
        self.assertIsNone(file_._deposition_date)
        self.assertIsNone(file_._classification)
        self.assertEqual(file_._keywords, [])
        self.assertEqual(file_._authors, [])
        self.assertIsNone(file_._technique)
        self.assertIsNone(file_._source_organism)
        self.assertIsNone(file_._expression_system)
        self.assertIsNone(file_._resolution)
        self.assertIsNone(file_._rvalue)
        self.assertIsNone(file_._rfree)
        self.assertEqual(file_._models, [])



class FileReprTests(TestCase):

    def test_generic_file_repr(self):
        file_ = File()
        self.assertEqual(repr(file_), "<File>")


    def test_file_file_repr(self):
        file_ = File()
        file_._filetype = "abc"
        self.assertEqual(repr(file_), "<.abc File>")


    def test_code_and_file_file_repr(self):
        file_ = File()
        file_._filetype = "abc"
        file_._code = "ABCD"
        self.assertEqual(repr(file_), "<ABCD.abc File>")



class FileFiletypeTests(TestCase):

    def test_can_get_filetype(self):
        file_ = File()
        file_._filetype = "1xxx"
        self.assertIs(file_._filetype, file_.filetype)


    def test_can_update_filetype(self):
        file_ = File()
        file_._filetype = "1xxx"
        file_.filetype = "2yyy"
        self.assertEqual(file_._filetype, "2yyy")



class FileCodeTests(TestCase):

    def test_can_get_code(self):
        file_ = File()
        file_._code = "1xxx"
        self.assertIs(file_._code, file_.code)


    def test_can_update_code(self):
        file_ = File()
        file_._code = "1xxx"
        file_.code = "2yyy"
        self.assertEqual(file_._code, "2yyy")



class FileTitleTests(TestCase):

    def test_can_get_title(self):
        file_ = File()
        file_._title = "TTT"
        self.assertIs(file_._title, file_.title)


    def test_can_update_title(self):
        file_ = File()
        file_._title = "TTT"
        file_.title = "TTTTTT"
        self.assertEqual(file_._title, "TTTTTT")



class FileDepositionDateTests(TestCase):

    def test_can_get_deposition_date(self):
        file_ = File()
        file_._deposition_date = date(1990, 9, 28)
        self.assertIs(file_._deposition_date, file_.deposition_date)


    def test_can_update_deposition_date(self):
        file_ = File()
        file_._deposition_date = date(1990, 9, 28)
        file_.deposition_date = date(1993, 9, 30)
        self.assertEqual(file_._deposition_date, date(1993, 9, 30))



class FileClassificationTests(TestCase):

    def test_can_get_classification(self):
        file_ = File()
        file_._classification = "CCC"
        self.assertIs(file_._classification, file_.classification)


    def test_can_update_classification(self):
        file_ = File()
        file_._classification = "CCC"
        file_.classification = "CCCCCC"
        self.assertEqual(file_._classification, "CCCCCC")



class FileKeywordsTests(TestCase):

    def test_can_get_keywords(self):
        file_ = File()
        file_._keywords = ["1", "2"]
        self.assertIs(file_._keywords, file_.keywords)



class FileAuthorsTests(TestCase):

    def test_can_get_authors(self):
        file_ = File()
        file_._authors = ["1", "2"]
        self.assertIs(file_._authors, file_.authors)



class FileTechniqueTests(TestCase):

    def test_can_get_technique(self):
        file_ = File()
        file_._technique = "ABC"
        self.assertIs(file_._technique, file_.technique)


    def test_can_update_technique(self):
        file_ = File()
        file_._technique = "ABC"
        file_.technique = "DEF"
        self.assertEqual(file_._technique, "DEF")



class FileSourceOrganismTests(TestCase):

    def test_can_get_source_organism(self):
        file_ = File()
        file_._source_organism = "ABC"
        self.assertIs(file_._source_organism, file_.source_organism)


    def test_can_update_source_organism(self):
        file_ = File()
        file_._source_organism = "ABC"
        file_.source_organism = "DEF"
        self.assertEqual(file_._source_organism, "DEF")



class FileExpressionSystemTests(TestCase):

    def test_can_get_expression_system(self):
        file_ = File()
        file_._expression_system = "ABC"
        self.assertIs(file_._expression_system, file_.expression_system)


    def test_can_update_expression_system(self):
        file_ = File()
        file_._expression_system = "ABC"
        file_.expression_system = "DEF"
        self.assertEqual(file_._expression_system, "DEF")



class FileResolutionTests(TestCase):

    def test_can_get_resolution(self):
        file_ = File()
        file_._resolution = 2.5
        self.assertIs(file_._resolution, file_.resolution)


    def test_can_update_resolution(self):
        file_ = File()
        file_._resolution = 2.5
        file_.resolution = 1.9
        self.assertEqual(file_._resolution, 1.9)



class FileRvalueTests(TestCase):

    def test_can_get_rvalue(self):
        file_ = File()
        file_._rvalue = 2.5
        self.assertIs(file_._rvalue, file_.rvalue)


    def test_can_update_rvalue(self):
        file_ = File()
        file_._rvalue = 2.5
        file_.rvalue = 1.9
        self.assertEqual(file_._rvalue, 1.9)



class FileRfreeTests(TestCase):

    def test_can_get_rfree(self):
        file_ = File()
        file_._rfree = 2.5
        self.assertIs(file_._rfree, file_.rfree)


    def test_can_update_rfree(self):
        file_ = File()
        file_._rfree = 2.5
        file_.rfree = 1.9
        self.assertEqual(file_._rfree, 1.9)



class FileAssembliesTests(TestCase):

    def test_can_get_assemblies(self):
        file_ = File()
        file_._assemblies = ["1", "2"]
        self.assertIs(file_._assemblies, file_.assemblies)



class FileModelsTests(TestCase):

    def test_can_get_models(self):
        file_ = File()
        file_._models = ["1", "2"]
        self.assertIs(file_._models, file_.models)



class FileModelTests(TestCase):

    def test_can_get_model(self):
        file_ = File()
        file_._models = ["1", "2"]
        self.assertIs(file_.model, file_._models[0])


    def test_can_get_no_model(self):
        file_ = File()
        self.assertIsNone(file_.model)



class PdbAssemblyTests(TestCase):

    def setUp(self):
        self.f = File()
        self.f._assemblies = [{
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
        self.f._models = [model, "MODEL"]
        self.chains, self.new_chains = chains, new_chains


    @patch("atomium.files.file.Model")
    def test_can_generate_assembly(self, mock_model):
        model = self.f.generate_assembly(1)
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
            self.f.generate_assembly(4)



class BestAssemblyTests(TestCase):

    def test_can_get_best_assembly(self):
        f = File()
        f._assemblies = [
         {"delta_energy": -20, "id": 1,},
         {"delta_energy": -60, "id": 2},
         {"delta_energy": -7, "id": 3},
         {"delta_energy": None, "id": 4}
        ]
        model = f.best_assembly
        self.assertEqual(model, {"delta_energy": -60, "id": 2})


    def test_can_return_none_as_best_assembly(self):
        f = File()
        model = f.best_assembly
        self.assertEqual(model, None)



class BestAssemblyGenerationTests(TestCase):

    @patch("atomium.files.file.File.generate_assembly")
    @patch("atomium.files.file.File.best_assembly", new_callable=PropertyMock)
    def test_can_generate_best_assembly(self, mock_best, mock_model):
        f = File()
        mock_best.return_value = {"delta_energy": -60, "id": 2}
        model = f.generate_best_assembly()
        self.assertIs(model, mock_model.return_value)
        mock_model.assert_called_with(2)


    @patch("atomium.files.file.File.best_assembly", new_callable=PropertyMock)
    def test_can_return_model_as_best_assembly(self, mock_best):
        mock_best.return_value = None
        f = File()
        f._models = "ABCDEF"
        model = f.generate_best_assembly()
        self.assertEqual(model, "A")



class FileSavingTests(TestCase):

    @patch("atomium.files.pdb.file_to_pdb_string")
    @patch("builtins.open")
    def test_can_save_as_pdb(self, mock_open, mock_string):
        open_return = MagicMock()
        mock_file = Mock()
        mock_write = MagicMock()
        mock_file.write = mock_write
        open_return.__enter__.return_value = mock_file
        mock_open.return_value = open_return
        f = File()
        f.save("/path/file.pdb")
        mock_string.assert_called_with(f)
        mock_open.assert_called_once_with("/path/file.pdb", "w")
        mock_write.assert_called_with(mock_string.return_value)


    @patch("atomium.files.xyz.file_to_xyz_string")
    @patch("builtins.open")
    def test_can_save_as_xyz(self, mock_open, mock_string):
        open_return = MagicMock()
        mock_file = Mock()
        mock_write = MagicMock()
        mock_file.write = mock_write
        open_return.__enter__.return_value = mock_file
        mock_open.return_value = open_return
        f = File()
        f.save("/path/file.xyz")
        mock_string.assert_called_with(f)
        mock_open.assert_called_once_with("/path/file.xyz", "w")
        mock_write.assert_called_with(mock_string.return_value)


    def test_unknown_file_extension(self):
        f = File()
        with self.assertRaises(ValueError):
            f.save("/path/file.abc")
