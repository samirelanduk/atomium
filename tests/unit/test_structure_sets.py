from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.base import StructureSet

class StructureSetTest(TestCase):

    def setUp(self):
        self.structures = [Mock(_id=1), Mock(_id=2), Mock(_id=2)]



class StructureSetCreationTests(StructureSetTest):

    def test_can_create_structure_set(self):
        s = StructureSet(*self.structures)
        self.assertEqual(s._d, {1: {self.structures[0]}, 2: set(self.structures[1:])})



class StructureSetAdditionTests(StructureSetTest):

    def test_can_add_structure_sets(self):
        s1 = StructureSet(*self.structures[:2])
        s2 = StructureSet(self.structures[2])
        s = s1 + s2
        self.assertEqual(s._d, {1: {self.structures[0]}, 2: set(self.structures[1:])})



class StructureSetLengthTests(StructureSetTest):

    def test_can_get_structure_set_len(self):
        s = StructureSet(*self.structures)
        self.assertEqual(len(s), 3)



class StructureSetAddingTests(StructureSetTest):

    def test_can_add_new_id_structure(self):
        s = StructureSet(*self.structures[1:])
        s.add(self.structures[0])
        self.assertEqual(s._d, {1: {self.structures[0]}, 2: set(self.structures[1:])})


    def test_can_add_existing_id_structure(self):
        s = StructureSet(*self.structures[:2])
        s.add(self.structures[2])
        self.assertEqual(s._d, {1: {self.structures[0]}, 2: set(self.structures[1:])})



class StructureSetRemovalTests(StructureSetTest):

    def test_can_remove_only_id_structure(self):
        s = StructureSet(*self.structures)
        s.remove(self.structures[0])
        self.assertEqual(s._d, {2: set(self.structures[1:])})


    def test_can_remove_duplicate_id_structure(self):
        s = StructureSet(*self.structures)
        s.remove(self.structures[2])
        self.assertEqual(s._d, {1: {self.structures[0]}, 2: {self.structures[1]}})



class StructureSetIdsTests(StructureSetTest):

    def test_can_get_ids(self):
        s = StructureSet(*self.structures)
        self.assertEqual(s.ids, {1, 2})



class StructureSetStructuresTests(StructureSetTest):

    def test_can_get_structures(self):
        s = StructureSet(*self.structures)
        self.assertEqual(set(s.structures), set(self.structures))



class StructureSetGettingTests(StructureSetTest):

    def test_can_get_structure_by_id(self):
        s = StructureSet(*self.structures)
        self.assertEqual(s.get(1), self.structures[0])
        self.assertIn(s.get(2), self.structures[1:])
