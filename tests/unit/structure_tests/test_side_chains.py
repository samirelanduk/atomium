from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.structures.molecules import SideChain, AtomicStructure
from atomium.structures.atoms import Atom

class SideChainTest(TestCase):

    def setUp(self):
        self.atom1, self.atom2, self.atom3 = Mock(Atom), Mock(Atom), Mock(Atom)
        self.atoms = [self.atom1, self.atom2, self.atom3]
        def mock_init(obj, *args, **kwargs):
            obj._atoms = set()
            obj._id_atoms = {}
        self.mock_init = mock_init



class SideCreationTests(SideChainTest):

    @patch("atomium.structures.molecules.AtomicStructure.__init__")
    def test_can_create_side_chain(self, mock_init):
        mock_init.side_effect = self.mock_init
        side_chain = SideChain(self.atom1, self.atom2, self.atom3)
        self.assertIsInstance(side_chain, AtomicStructure)
        mock_init.assert_called_with(side_chain, self.atom1, self.atom2, self.atom3)
        self.assertEqual(side_chain._occupancy, 1)


    @patch("atomium.structures.molecules.AtomicStructure.__init__")
    def test_can_create_side_chain_with_occupancy(self, mock_init):
        mock_init.side_effect = self.mock_init
        side_chain = SideChain(self.atom1, self.atom2, self.atom3, occupancy=0.5)
        self.assertIsInstance(side_chain, AtomicStructure)
        mock_init.assert_called_with(side_chain, self.atom1, self.atom2, self.atom3)
        self.assertEqual(side_chain._occupancy, 0.5)


    def test_occupancy_must_be_number(self):
        with self.assertRaises(TypeError):
            SideChain(self.atom1, self.atom2, self.atom3, occupancy="o")


    def test_occupancy_must_be_in_range(self):
        with self.assertRaises(ValueError):
            SideChain(self.atom1, self.atom2, self.atom3, occupancy=0)
        with self.assertRaises(ValueError):
            SideChain(self.atom1, self.atom2, self.atom3, occupancy=1.001)



class SideChainReprTests(SideChainTest):

    def test_side_chain_repr(self):
        side_chain = SideChain(self.atom1, self.atom2, self.atom3)
        self.assertEqual(str(side_chain), "<SideChain (3 atoms)>")
