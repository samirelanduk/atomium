from datetime import date
from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.data import File

class FileCreationTests(TestCase):

    def test_can_create_file(self):
        f = File("abc")
        self.assertEqual(f._filetype, "abc")
        self.assertEqual(f._models, [])



class FileReprTests(TestCase):

    def test_generic_file_repr(self):
        f = File()
        self.assertEqual(repr(f), "<File>")


    def test_generic_file_repr(self):
        f = File("abc")
        f._code = None
        self.assertEqual(repr(f), "<.abc File>")


    def test_code_file_repr(self):
        f = File("abc")
        f._code = "ABCD"
        self.assertEqual(repr(f), "<ABCD.abc File>")



class FileTypePropertyTests(TestCase):

    def test_file_type(self):
        f = File("abc")
        self.assertIs(f.filetype, f._filetype)



class TitlePropertyTests(TestCase):

    def test_title(self):
        f = File("abc")
        f._title = "XXX"
        self.assertIs(f.title, f._title)



class DepositionDatePropertyTests(TestCase):

    def test_deposition_date(self):
        f = File("abc")
        f._deposition_date = "XXX"
        self.assertIs(f.deposition_date, f._deposition_date)



class ClassificationPropertyTests(TestCase):

    def test_classification(self):
        f = File("abc")
        f._classification = "XXX"
        self.assertIs(f.classification, f._classification)



class KeywordsPropertyTests(TestCase):

    def test_keywords(self):
        f = File("abc")
        f._keywords = "XXX"
        self.assertIs(f.keywords, f._keywords)



class AuthorsPropertyTests(TestCase):

    def test_authors(self):
        f = File("abc")
        f._authors = "XXX"
        self.assertIs(f.authors, f._authors)



class TechniquePropertyTests(TestCase):

    def test_technique(self):
        f = File("abc")
        f._technique = "XXX"
        self.assertIs(f.technique, f._technique)



class SourceOrganismPropertyTests(TestCase):

    def test_source_organism(self):
        f = File("abc")
        f._source_organism = "XXX"
        self.assertIs(f.source_organism, f._source_organism)



class ExpressionSystemPropertyTests(TestCase):

    def test_expression_system(self):
        f = File("abc")
        f._expression_system = "XXX"
        self.assertIs(f.expression_system, f._expression_system)



class MissingResiduesPropertyTests(TestCase):

    def test_missing_residues(self):
        f = File("abc")
        f._missing_residues = "XXX"
        self.assertIs(f.missing_residues, f._missing_residues)



class ResolutionPropertyTests(TestCase):

    def test_resolution(self):
        f = File("abc")
        f._resolution = "XXX"
        self.assertIs(f.resolution, f._resolution)



class RvaluePropertyTests(TestCase):

    def test_rvalue(self):
        f = File("abc")
        f._rvalue = "XXX"
        self.assertIs(f.rvalue, f._rvalue)



class RfreePropertyTests(TestCase):

    def test_rfree(self):
        f = File("abc")
        f._rfree = "XXX"
        self.assertIs(f.rfree, f._rfree)



class AssemblyPropertyTests(TestCase):

    def test_assemblies(self):
        f = File("abc")
        f._assemblies = "XXX"
        self.assertIs(f.assemblies, f._assemblies)



class ModelsPropertyTests(TestCase):

    def test_models(self):
        f = File("abc")
        f._models = "XXX"
        self.assertIs(f.models, f._models)



class ModelPropertyTests(TestCase):

    def test_model(self):
        f = File("abc")
        f._models = "XYZ"
        self.assertIs(f.model, f._models[0])
