from collections import OrderedDict
from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.base import *

class ObjectFromFilterTests(TestCase):

    def test_can_get_same_object(self):
        obj = Mock()
        obj.__lt__ = MagicMock()
        obj2 = get_object_from_filter(obj, ["height"])
        self.assertIs(obj, obj2)
        obj2 = get_object_from_filter(obj, ["height", "regex"])
        self.assertIs(obj, obj2)
        obj2 = get_object_from_filter(obj, ["height", "lt"])
        self.assertIs(obj, obj2)
    

    def test_can_get_chained_object(self):
        obj = Mock()
        obj2 = get_object_from_filter(obj, ["o1", "o2", "o3", "height", "regex"])
        self.assertIs(obj2, obj.o1.o2.o3)
        obj2 = get_object_from_filter(obj, ["o1", "o2", "o3", "height"])
        self.assertIs(obj2, obj.o1.o2.o3)



class AttributeGettingTests(TestCase):

    def test_can_get_basic_attribute(self):
        obj = Mock(x=10)
        self.assertEqual(get_object_attribute_from_filter(obj, ["x"]), 10)
        self.assertEqual(get_object_attribute_from_filter(obj, ["y", "x"]), 10)
    

    def test_can_get_attribute_from_early_chain(self):
        obj = Mock(x=10)
        del obj.regex
        self.assertEqual(get_object_attribute_from_filter(obj, ["x", "regex"]), 10)
    

    def test_can_get_no_attribute(self):
        obj = Mock(x=10)
        del obj.y
        self.assertIsNone(get_object_attribute_from_filter(obj, ["y"]))



class AttributeMatchingTests(TestCase):

    def test_exact_match(self):
        self.assertTrue(attribute_matches_value(10, 10, ["height", "xy"]))
        self.assertFalse(attribute_matches_value(10, 11, ["height", "xy"]))
    

    def test_regex_match(self):
        self.assertTrue(attribute_matches_value("jon", "jon|joe", ["name", "regex"]))
        self.assertFalse(attribute_matches_value("jon", "jon|joe", ["name", "rogox"]))


    def test_magic_method_match(self):
        self.assertTrue(attribute_matches_value(12, 10, ["height", "gt"]))
        self.assertFalse(attribute_matches_value(10, 10, ["height", "gt"]))
        self.assertTrue(attribute_matches_value(10, 10, ["height", "gte"]))



class ObjectFilteringTests(TestCase):

    @patch("atomium.base.get_object_from_filter")
    @patch("atomium.base.get_object_attribute_from_filter")
    @patch("atomium.base.attribute_matches_value")
    @patch("atomium.base.StructureSet")
    def test_can_filter_objects(self, mock_s, mock_match, mock_getat, mock_getob):
        structures=[
         Mock(x="A", y=1), Mock(x="B", y=3), Mock(x="B", y=3),
         Mock(x="C", y=2), Mock(x="D", y=4), Mock(x="D", y=4)
        ]
        objects = Mock(structures=structures)
        mock_getob.side_effect = lambda s, c: s
        mock_getat.side_effect = lambda s, c: c[0]
        mock_match.side_effect = [False, True, False, True, False, False]
        filter_objects(objects, "key__key2__key_3", "value")
        for structure in structures:
            mock_getob.assert_any_call(structure, ["key", "key2", "key_3"])
            mock_getat.assert_any_call(structure, ["key", "key2", "key_3"])
            mock_match.assert_any_call("key", "value", ["key", "key2", "key_3"])
        mock_s.assert_called_with(structures[1], structures[3])



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



class StructureSetTests(TestCase):

    def test_can_make_structure_set(self):
        objects = [Mock(_id=n) for n in range(5)]
        s = StructureSet(*objects)
        self.assertEqual(s._d, {
         0: {objects[0]}, 1: {objects[1]}, 2: {objects[2]},
         3: {objects[3]}, 4: {objects[4]}
        })
        objects[2]._id = 0
        s = StructureSet(*objects)
        self.assertEqual(s._d, {
         0: {objects[0], objects[2]}, 1: {objects[1]},
         3: {objects[3]}, 4: {objects[4]}
        })
    

    def test_can_add_two_structure_sets(self):
        objects = [Mock(_id=n) for n in range(5)]
        objects[2]._id = 0
        s1 = StructureSet(*objects[:3])
        s2 = StructureSet(*objects[3:])
        self.assertEqual(s1._d, {
         0: {objects[0], objects[2]}, 1: {objects[1]},
        })
        self.assertEqual(s2._d, {3: {objects[3]}, 4: {objects[4]}})
        s3 = s1 + s2
        self.assertEqual(s3._d, {
         0: {objects[0], objects[2]}, 1: {objects[1]},
         3: {objects[3]}, 4: {objects[4]}
        })
    

    def test_can_get_length_of_structure_sets(self):
        objects = [Mock(_id=n) for n in range(5)]
        s = StructureSet(*objects)
        self.assertEqual(len(s), 5)
        objects[2]._id = 0
        s = StructureSet(*objects)
        self.assertEqual(len(s), 5)
    

    def test_can_get_structure_set_ids(self):
        objects = [Mock(_id=n) for n in range(5)]
        s = StructureSet(*objects)
        self.assertEqual(s.ids, {0, 1, 2, 3, 4})
    

    def test_can_get_structure_set_structures(self):
        objects = [Mock(_id=n) for n in range(5)]
        s = StructureSet(*objects)
        self.assertEqual(s.structures, objects)
        objects[2]._id = 0
        s = StructureSet(*objects)
        self.assertEqual(set(s.structures), set(objects))
    

    def test_can_get_structures_by_id(self):
        objects = [Mock(_id=n) for n in range(5)]
        s = StructureSet(*objects)
        self.assertEqual(s.get(0), objects[0])
        self.assertEqual(s.get(4), objects[4])
        self.assertEqual(s.get(5), None)
        objects[2]._id = 0
        s = StructureSet(*objects)
        self.assertIn(s.get(0), (objects[0], objects[2]))
        self.assertEqual(s.get(4), objects[4])
        self.assertEqual(s.get(2), None)