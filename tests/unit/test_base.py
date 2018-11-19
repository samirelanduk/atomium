from collections import OrderedDict
from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.base import *

class ObjectFilteringTests(TestCase):

    def setUp(self):

        self.objects = Mock(structures=[
         Mock(x="A", y=1), Mock(x="B", y=3), Mock(x="B", y=3),
         Mock(x="C", y=2), Mock(x="D", y=4), Mock(x="D", y=4)
        ])


    @patch("atomium.base.StructureSet")
    def test_can_filter_by_property(self, mock_set):
        objects = filter_objects(self.objects, "x", "B")
        self.assertEqual(objects, mock_set.return_value)
        mock_set.assert_called_with(*self.objects.structures[1:3])
        objects = filter_objects(self.objects, "y", 2)
        self.assertEqual(objects, mock_set.return_value)
        mock_set.assert_called_with(self.objects.structures[3])
        objects = filter_objects(self.objects, "y", 20)
        self.assertEqual(objects, mock_set.return_value)
        mock_set.assert_called_with()



    @patch("atomium.base.StructureSet")
    def test_can_filter_by_property_regex(self, mock_set):
        objects = filter_objects(self.objects, "x__regex", "B|A")
        self.assertEqual(objects, mock_set.return_value)
        mock_set.assert_called_with(*self.objects.structures[:3])
        objects = filter_objects(self.objects, "x__regex", ".")
        self.assertEqual(objects, mock_set.return_value)
        mock_set.assert_called_with(*self.objects.structures)


    @patch("atomium.base.StructureSet")
    def test_can_filter_by_property_comparison(self, mock_set):
        objects = filter_objects(self.objects, "y__gt", 3)
        self.assertEqual(objects, mock_set.return_value)
        mock_set.assert_called_with(*self.objects.structures[4:])
        objects = filter_objects(self.objects, "y__ge", 3)
        self.assertEqual(objects, mock_set.return_value)
        mock_set.assert_called_with(*self.objects.structures[1:3] + self.objects.structures[4:])
        objects = filter_objects(self.objects, "y__ne", 4)
        self.assertEqual(objects, mock_set.return_value)
        mock_set.assert_called_with(*self.objects.structures[:4])



class QueryDecoratorTests(TestCase):

    def setUp(self):
        self.s = Mock(structures={2, 4, 6}, ids={1, 3, 5})
        self.f = lambda s: self.s


    def test_can_get_unfiltered_objects(self):
        f = query(self.f)
        self.assertEqual(f(self), {2, 4, 6})


    @patch("atomium.base.filter_objects")
    def test_can_get_filtered_objects(self, mock_filter):
        mock_filter.side_effect = [Mock(structures={20}, ids={10})]
        f = query(self.f)
        self.assertEqual(f(self, a=1), {20})
        mock_filter.assert_any_call(self.s, "a", 1)


    @patch("atomium.base.filter_objects")
    def test_can_get_filtered_objects_as_tuple(self, mock_filter):
        mock_filter.side_effect = [Mock(structures={2}, ids={1})]
        f = query(self.f, tuple_=True)
        self.assertEqual(f(self, a=1), (2,))
        mock_filter.assert_any_call(self.s, "a", 1)


    def test_can_get_objects_by_id(self):
        f = query(self.f)
        self.assertEqual(f(self, 3), {self.s.get.return_value})
        self.s.get.assert_called_with(3)
        self.assertEqual(f(self, 8), set())



class GetOneDecoratorTests(TestCase):

    def test_can_get_one(self):
        f = lambda s: [4, 6, 7]
        f = getone(f)
        self.assertEqual(f(self), 4)


    def test_can_get_mone(self):
        f = lambda s: []
        f = getone(f)
        self.assertEqual(f(self), None)



class StructureClassMetaclassTests(TestCase):

    @patch("atomium.base.query")
    @patch("atomium.base.getone")
    def test_structure_class_metaclass(self, mock_getone, mock_query):
        class TestClass(metaclass=StructureClass):
            def a(self): return 1000
            def chains(self): return {1: 2, 3: 4}
            def residues(self): return {10: 2, 30: 4}
            def ligands(self): return {11: 2, 31: 4}
            def waters(self): return {12: 2, 32: 4}
            def molecules(self): return {13: 2, 33: 4}
            def atoms(self): return {14: 2, 34: 4}
            def b(self): return 2000
        obj = TestClass()
        self.assertIs(obj.chains, mock_query.return_value)
        self.assertIs(obj.chain, mock_getone.return_value)
        self.assertIs(obj.residues, mock_query.return_value)
        self.assertIs(obj.residue, mock_getone.return_value)
        self.assertIs(obj.ligands, mock_query.return_value)
        self.assertIs(obj.ligand, mock_getone.return_value)
        self.assertIs(obj.waters, mock_query.return_value)
        self.assertIs(obj.water, mock_getone.return_value)
        self.assertIs(obj.molecules, mock_query.return_value)
        self.assertIs(obj.molecule, mock_getone.return_value)
        self.assertIs(obj.atoms, mock_query.return_value)
        self.assertIs(obj.atom, mock_getone.return_value)
        self.assertEqual(obj.a(), 1000)
        self.assertEqual(obj.b(), 2000)
