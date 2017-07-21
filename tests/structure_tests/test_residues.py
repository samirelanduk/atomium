from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.structures.molecules import Molecule, Residue
from atomium.structures.atoms import Atom

class ResidueTest(TestCase):

    def setUp(self):
        self.atom1, self.atom2, self.atom3 = Mock(Atom), Mock(Atom), Mock(Atom)
        self.atoms = [self.atom1, self.atom2, self.atom3]



class ResidueCreationTests(ResidueTest):

    @patch("atomium.structures.molecules.Molecule.__init__")
    def test_can_create_residue(self, mock_init):
        mock_init.return_value = None
        res = Residue(self.atom1, self.atom2, self.atom3)
        self.assertIsInstance(res, Molecule)
        mock_init.assert_called_with(res, self.atom1, self.atom2, self.atom3)


    @patch("atomium.structures.molecules.Molecule.__init__")
    def test_can_create_residue_with_name(self, mock_init):
        mock_init.return_value = None
        res = Residue(self.atom1, self.atom2, self.atom3, name="GLY")
        mock_init.assert_called_with(res, self.atom1, self.atom2, self.atom3, name="GLY")


    @patch("atomium.structures.molecules.Molecule.__init__")
    def test_can_residue_init_doesnt_send_residue_id(self, mock_init):
        mock_init.return_value = None
        res = Residue(self.atom1, self.atom2, self.atom3, residue_id="B100")
        mock_init.assert_called_with(res, self.atom1, self.atom2, self.atom3)


    def test_can_create_residue_with_id(self):
        res = Residue(self.atom1, self.atom2, self.atom3, residue_id="A101")
        self.assertEqual(res._id, "A101")


    def test_residue_id_must_be_str(self):
        with self.assertRaises(TypeError):
            Residue(self.atom1, self.atom2, self.atom3, residue_id=100)


    def test_residue_id_must_be_unique(self):
        res = Residue(self.atom1, self.atom2, self.atom3, residue_id="A200")
        with self.assertRaises(ValueError):
            Residue(self.atom1, self.atom2, self.atom3, residue_id="A200")


    def test_atoms_are_linked_to_residue(self):
        res = Residue(self.atom1, self.atom2, self.atom3)
        self.assertIs(self.atom1._residue, res)
        self.assertIs(self.atom2._residue, res)
        self.assertIs(self.atom3._residue, res)



class ResiudeIdChangesTests(ResidueTest):

    def test_changing_id_updates_known_ids(self):
        res = Residue(self.atom1, self.atom2, self.atom3, residue_id="A220")
        self.assertIn("A220", Molecule.known_ids)
        res._id = "A221"
        self.assertIn("A221", Molecule.known_ids)
        self.assertNotIn("A220", Molecule.known_ids)


    def test_deleting_a_molecule_frees_up_its_id(self):
        res = Residue(self.atom1, self.atom2, self.atom3, residue_id="A225")
        self.assertIn("A225", Molecule.known_ids)
        res = Residue(self.atom1, self.atom2, self.atom3, residue_id="A226")
        self.assertNotIn("A225", Molecule.known_ids)
        self.assertIn("A226", Molecule.known_ids)



class ResidueReprTests(ResidueTest):

    def test_molecule_repr_no_id_or_name(self):
        res = Residue(self.atom1, self.atom2, self.atom3)
        self.assertEqual(str(res), "<Residue (3 atoms)>")


    def test_molecule_repr_id_no_name(self):
        res = Residue(self.atom1, self.atom2, self.atom3, residue_id="C10A")
        self.assertEqual(str(res), "<Residue C10A (3 atoms)>")


    def test_molecule_repr_name_no_id(self):
        res = Residue(self.atom1, self.atom2, self.atom3, name="GLY")
        self.assertEqual(str(res), "<Residue (GLY, 3 atoms)>")


    def test_molecule_repr_id_and_name(self):
        res = Residue(self.atom1, self.atom2, residue_id="C10B", name="GLY")
        self.assertEqual(str(res), "<Residue C10B (GLY, 2 atoms)>")



class ResidueIdTests(ResidueTest):

    def test_residue_id_property(self):
        res = Residue(self.atom1, self.atom2, self.atom3, residue_id="C10C")
        self.assertIs(res._id, res.residue_id())


    def test_can_update_residue_id(self):
        res = Residue(self.atom1, self.atom2, self.atom3, residue_id="C10D")
        res.residue_id("A20")
        self.assertEqual(res._id, "A20")


    def test_residue_id_must_be_str(self):
        res = Residue(self.atom1, self.atom2, self.atom3, residue_id="C10E")
        with self.assertRaises(TypeError):
            res.residue_id(10)
