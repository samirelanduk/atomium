from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.structures import Residue, Atom, Molecule, Het, AtomStructure

class ResidueTest(TestCase):

    def setUp(self):
        self.atoms = [Mock(Atom), Mock(Atom)]



class ResidueCreationTests(ResidueTest):

    @patch("atomium.structures.Het.__init__")
    def test_can_create_residue(self, mock_init):
        residue = Residue(*self.atoms)
        self.assertIsNone(residue._id)
        self.assertIsNone(residue._name)
        self.assertIsNone(residue._next)
        self.assertIsNone(residue._previous)
        self.assertIsNone(residue._chain)
        mock_init.assert_called_with(residue, *self.atoms)


    @patch("atomium.structures.Het.__init__")
    def test_can_create_residue_with_attributes(self, mock_init):
        residue = Residue(*self.atoms, id="A1", name="VAL")
        self.assertEqual(residue._id, "A1")
        self.assertEqual(residue._name, "VAL")
        mock_init.assert_called_with(residue, *self.atoms)



class ResidueReprTests(ResidueTest):

    def test_residue_repr(self):
        residue = Residue(*self.atoms, id="A1", name="VAL")
        self.assertEqual(repr(residue), "<Residue VAL (A1)>")



class ResidueAtomsTests(ResidueTest):

    def test_can_get_atoms(self):
        residue = Residue(*self.atoms, id="A1", name="VAL")
        self.assertEqual(residue.atoms(), set(self.atoms))



class ResidueNextPropertyTests(ResidueTest):

    def test_can_get_next_residue(self):
        residue = Residue(*self.atoms, id="A1", name="VAL")
        residue._next = "NEXT"
        self.assertIs(residue.next, residue._next)


    def test_can_set_next_residue(self):
        other = Mock(Residue, _previous=None)
        residue = Residue(*self.atoms, id="A1", name="VAL")
        residue.next = other
        self.assertEqual(residue._next, other)
        self.assertEqual(other._previous, residue)


    def test_can_unset_next_residue(self):
        residue = Residue(*self.atoms, id="A1", name="VAL")
        other = Mock(Residue, _previous=residue)
        residue._next = other
        residue.next = None
        self.assertIsNone(residue._next)
        self.assertIsNone(other._previous)


    def test_cant_connect_residue_to_self(self):
        residue = Residue(*self.atoms, id="A1", name="VAL")
        with self.assertRaises(ValueError):
            residue.next = residue



class ResiduePreviousPropertyTests(ResidueTest):

    def test_can_get_previous_residue(self):
        residue = Residue(*self.atoms, id="A1", name="VAL")
        residue._previous = "PREV"
        self.assertIs(residue.previous, residue._previous)


    def test_can_set_previous_residue(self):
        other = Mock(Residue, _next=None)
        residue = Residue(*self.atoms, id="A1", name="VAL")
        residue.previous = other
        self.assertEqual(residue._previous, other)
        self.assertEqual(other._next, residue)


    def test_can_unset_previous_residue(self):
        residue = Residue(*self.atoms, id="A1", name="VAL")
        other = Mock(Residue, _next=residue)
        residue._previous = other
        residue.previous = None
        self.assertIsNone(residue._previous)
        self.assertIsNone(other._next)


    def test_cant_connect_residue_to_self(self):
        residue = Residue(*self.atoms, id="A1", name="VAL")
        with self.assertRaises(ValueError):
            residue.previous = residue



class ResidueCodeTests(ResidueTest):

    def test_can_get_residue_code(self):
        residue = Residue(*self.atoms, id="A1", name="VAL")
        self.assertEqual(residue.code, "V")


    def test_can_get_no_code(self):
        residue = Residue(*self.atoms, id="A1", name="UNKNOWN")
        self.assertEqual(residue.code, "X")



class ResidueFullNameTests(ResidueTest):

    def test_can_get_residue_full_name(self):
        residue = Residue(*self.atoms, id="A1", name="VAL")
        self.assertEqual(residue.full_name, "valine")


    def test_can_get_no_code(self):
        residue = Residue(*self.atoms, id="A1", name="UNKNOWN")
        self.assertEqual(residue.full_name, "UNKNOWN")



class ResidueCopyingTests(ResidueTest):

    @patch("atomium.structures.Residue.atoms")
    def test_can_copy_residue(self, mock_atoms):
        mock_atoms.return_value = self.atoms
        residue = Residue(*self.atoms, id="A1", name="VAL")
        p = patch("atomium.structures.Residue")
        mock_residue = p.start()
        try:
            copy = residue.copy()
            mock_residue.assert_called_with(
             self.atoms[0].copy.return_value, self.atoms[1].copy.return_value,
             id="A1", name="VAL"
            )
        finally: p.stop()
