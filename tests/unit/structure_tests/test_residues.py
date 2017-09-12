from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.structures.molecules import Molecule, Residue
from atomium.structures.atoms import Atom

class ResidueTest(TestCase):

    def setUp(self):
        self.atom1, self.atom2, self.atom3 = Mock(Atom), Mock(Atom), Mock(Atom)
        self.atoms = [self.atom1, self.atom2, self.atom3]
        def mock_init(obj, *args, **kwargs):
            obj._atoms = set()
            obj._id_atoms = {}
        self.mock_init = mock_init



class ResidueCreationTests(ResidueTest):

    @patch("atomium.structures.molecules.Molecule.__init__")
    def test_can_create_residue(self, mock_init):
        mock_init.side_effect = self.mock_init
        res = Residue(self.atom1, self.atom2, self.atom3)
        self.assertIsInstance(res, Molecule)
        mock_init.assert_called_with(res, self.atom1, self.atom2, self.atom3)
        self.assertIsNone(res._next)
        self.assertIsNone(res._previous)


    @patch("atomium.structures.molecules.Molecule.__init__")
    def test_can_create_residue_with_name(self, mock_init):
        mock_init.side_effect = self.mock_init
        res = Residue(self.atom1, self.atom2, self.atom3, name="GLY")
        mock_init.assert_called_with(res, self.atom1, self.atom2, self.atom3, name="GLY")


    @patch("atomium.structures.molecules.Molecule.__init__")
    def test_can_create_residue_with_id(self, mock_init):
        mock_init.side_effect = self.mock_init
        res = Residue(self.atom1, self.atom2, self.atom3, residue_id="A1")
        mock_init.assert_called_with(res, self.atom1, self.atom2, self.atom3, molecule_id="A1")


    def test_atoms_are_linked_to_residue(self):
        res = Residue(self.atom1, self.atom2, self.atom3)
        self.assertIs(self.atom1._residue, res)
        self.assertIs(self.atom2._residue, res)
        self.assertIs(self.atom3._residue, res)



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



class ResidueAtomAdditionTests(ResidueTest):

    def test_adding_atom_updates_atom(self):
        res = Residue(self.atom1, self.atom2)
        res.add_atom(self.atom3)
        self.assertIs(self.atom3._residue, res)



class ResidueAtomRemovalTests(ResidueTest):

    def test_removing_atom_updates_atom(self):
        res = Residue(self.atom1, self.atom2, self.atom3)
        res.remove_atom(self.atom3)
        self.assertIs(self.atom3._residue, None)



class ResidueNextTests(ResidueTest):

    def test_next_gets_next(self):
        res = Residue(self.atom1, self.atom2, self.atom3)
        res._next = Mock(Residue)
        self.assertIs(res.next(), res._next)


    def test_can_assign_next(self):
        res = Residue(self.atom1, self.atom2, self.atom3)
        next_res = Mock(Residue)
        res.next(next_res)
        self.assertIs(res._next, next_res)
        self.assertIs(next_res._previous, res)


    def test_next_residue_must_be_residue(self):
        res = Residue(self.atom1, self.atom2, self.atom3)
        mol = Mock(Molecule)
        with self.assertRaises(TypeError):
            res.next(mol)


    def test_next_res_cannot_be_self(self):
        res = Residue(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(ValueError):
            res.next(res)


    def test_can_remove_next_residue(self):
        res = Residue(self.atom1, self.atom2, self.atom3)
        next_res = Mock(Residue)
        res._next = next_res
        next_res._previous = res
        res.next(None)
        self.assertIsNone(res._next)
        self.assertIsNone(next_res._previous)



class ResiduePreviousTests(ResidueTest):

    def test_previous_gets_previous(self):
        res = Residue(self.atom1, self.atom2, self.atom3)
        res._previous = Mock(Residue)
        self.assertIs(res.previous(), res._previous)


    def test_can_assign_previous(self):
        res = Residue(self.atom1, self.atom2, self.atom3)
        previous_res = Mock(Residue)
        res.previous(previous_res)
        self.assertIs(res._previous, previous_res)
        self.assertIs(previous_res._next, res)


    def test_previous_residue_must_be_residue(self):
        res = Residue(self.atom1, self.atom2, self.atom3)
        mol = Mock(Molecule)
        with self.assertRaises(TypeError):
            res.previous(mol)


    def test_previous_res_cannot_be_self(self):
        res = Residue(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(ValueError):
            res.previous(res)


    def test_can_remove_previous_residue(self):
        res = Residue(self.atom1, self.atom2, self.atom3)
        previous_res = Mock(Residue)
        res._previous = previous_res
        previous_res._next = res
        res.previous(None)
        self.assertIsNone(res._previous)
        self.assertIsNone(previous_res._next)



class ResidueChainTests(ResidueTest):

    @patch("atomium.structures.molecules.Residue.atoms")
    def test_can_get_chain(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms)
        chain = Mock()
        self.atom1.chain.return_value = chain
        self.atom2.chain.return_value = chain
        self.atom3.chain.return_value = chain
        res = Residue(self.atom1, self.atom2, self.atom3)
        self.assertIs(res.chain(), chain)


    @patch("atomium.structures.molecules.Residue.atoms")
    def test_can_get_no_chain(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms)
        self.atom1.chain.return_value = None
        self.atom2.chain.return_value = None
        self.atom3.chain.return_value = None
        res = Residue(self.atom1, self.atom2, self.atom3)
        self.assertIs(res.chain(), None)
