from unittest import TestCase
from unittest.mock import patch, Mock, PropertyMock
from atomium.models.molecules import Het, Residue

class HetCreationTests(TestCase):

    def test_can_create_residue(self):
        residue = Residue()
        self.assertIsInstance(residue, Het)


    @patch("atomium.models.molecules.Het.__init__")
    def test_can_create_residue(self, mock_init):
        residue = Residue("a", b="c")
        self.assertIsInstance(residue, Het)
        mock_init.assert_called_with(residue, "a", b="c")
        self.assertIsNone(residue._next)
        self.assertIsNone(residue._previous)



class ResidueFullNameTests(TestCase):

    def test_can_get_name(self):
        res = Residue()
        self.assertIsNone(res.full_name)
        res._name = "XMP"
        self.assertEqual(res.full_name, "XMP")


    def test_can_lookup_full_name(self):
        res = Residue()
        res._name = "met"
        self.assertEqual(res.full_name, "methionine")



class ResidueCodeTests(TestCase):

    def test_can_get_generic_code(self):
        res = Residue()
        self.assertIsNone(res.code)
        res._name = "PMP"
        self.assertEqual(res.code, "X")


    def test_can_lookup_code(self):
        res = Residue()
        res._name = "met"
        self.assertEqual(res.code, "M")



class ResidueNextTests(TestCase):

    def test_next_gets_next(self):
        res = Residue()
        res._next = Mock(Residue)
        self.assertIs(res.next, res._next)


    def test_can_assign_next(self):
        res = Residue()
        next_res = Mock(Residue)
        res.next = next_res
        self.assertIs(res._next, next_res)
        self.assertIs(next_res._previous, res)


    def test_next_res_cannot_be_self(self):
        res = Residue()
        with self.assertRaises(ValueError):
            res.next = res


    def test_can_remove_next_residue(self):
        res = Residue()
        next_res = Mock(Residue)
        res._next = next_res
        next_res._previous = res
        res.next = None
        self.assertIsNone(res._next)
        self.assertIsNone(next_res._previous)



class ResiduePreviousTests(TestCase):

    def test_previous_gets_previous(self):
        res = Residue()
        res._previous = Mock(Residue)
        self.assertIs(res.previous, res._previous)


    def test_can_assign_previous(self):
        res = Residue()
        previous_res = Mock(Residue)
        res.previous = previous_res
        self.assertIs(res._previous, previous_res)
        self.assertIs(previous_res._next, res)


    def test_previous_res_cannot_be_self(self):
        res = Residue()
        with self.assertRaises(ValueError):
            res.previous = res


    def test_can_remove_previous_residue(self):
        res = Residue()
        previous_res = Mock(Residue)
        res._previous = previous_res
        previous_res._next = res
        res.previous = None
        self.assertIsNone(res._previous)
        self.assertIsNone(previous_res._next)
