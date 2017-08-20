from unittest import TestCase
from unittest.mock import patch, Mock, MagicMock
from atomium.structures.models import ChainStructure
from atomium.structures.chains import Chain
from atomium.structures.atoms import Atom

class ChainStructureTest(TestCase):

    def setUp(self):
        self.structure = ChainStructure()
        self.atom1, self.atom2 = Mock(Atom), Mock(Atom)
        self.atom3, self.atom4 = Mock(Atom), Mock(Atom)
        self.chain1, self.chain2 = Mock(Chain), Mock(Chain)
        self.chain1.chain_id.return_value = "C"
        self.chain2.chain_id.return_value = "D"
        self.chain1.name.return_value = "III"
        self.chain2.name.return_value = "JJJ"
        self.atom1.chain.return_value = self.chain1
        self.atom2.chain.return_value = self.chain1
        self.atom3.chain.return_value = self.chain2
        self.atom4.chain.return_value = self.chain2
        self.chain1.atoms.return_value = set([self.atom1, self.atom2])
        self.chain2.atoms.return_value = set([self.atom3, self.atom4])
        self.structure.atoms = lambda: set([
         self.atom1, self.atom2, self.atom3, self.atom4
        ])
        self.structure.add_atom = MagicMock()
        self.structure.remove_atom = MagicMock()


class ChainStructureChainsTests(ChainStructureTest):

    def test_can_get_chains(self):
        self.assertEqual(
         self.structure.chains(),
         set([self.chain1, self.chain2])
        )


    def test_can_filter_none_from_chains(self):
        self.atom4.chain.return_value = None
        self.assertEqual(
         self.structure.chains(),
         set([self.chain1, self.chain2])
        )


    def test_can_get_chains_by_id(self):
        self.assertEqual(
         self.structure.chains(chain_id="C"), set([self.chain1])
        )
        self.assertEqual(
         self.structure.chains(chain_id="D"), set([self.chain2])
        )
        self.assertEqual(self.structure.chains(chain_id="E"), set())


    def test_can_get_chains_by_name(self):
        self.assertEqual(
         self.structure.chains(name="III"), set([self.chain1])
        )
        self.assertEqual(
         self.structure.chains(name="JJJ"), set([self.chain2])
        )
        self.assertEqual(self.structure.chains(name="GLY"), set())



class ChainStructureChainTests(ChainStructureTest):

    @patch("atomium.structures.models.ChainStructure.chains")
    def test_chain_calls_chains(self, mock_chains):
        mock_chains.return_value = set([self.chain2])
        chain = self.structure.chain(name="A")
        mock_chains.assert_called_with(name="A")
        self.assertIs(chain, self.chain2)


    @patch("atomium.structures.models.ChainStructure.chains")
    def test_chain_can_return_none(self, mock_chains):
        mock_chains.return_value = set()
        self.assertIs(self.structure.chain(name="E"), None)


    @patch("atomium.structures.models.ChainStructure.chains")
    def test_chain_can_get_chain_by_id_and_name(self, mock_chains):
        mock_chains.return_value = set([self.chain1])
        chain = self.structure.chain(chain_id="A1", name="A")
        mock_chains.assert_called_with(chain_id="A1", name="A")
        self.assertIs(chain, self.chain1)



class ChainStructureChainAdditionTests(ChainStructureTest):

    def test_can_add_chain(self):
        self.structure._atoms = set([
         self.atom1, self.atom2
        ])
        self.structure.add_chain(self.chain2)
        self.structure.add_atom.assert_any_call(self.atom3)
        self.structure.add_atom.assert_any_call(self.atom4)


    def test_can_only_add_chains(self):
        with self.assertRaises(TypeError):
            self.structure.add_chain("self.chain4")



class ChainStructureChainRemovalTests(ChainStructureTest):

    def test_can_remove_chain(self):
        self.structure.remove_chain(self.chain2)
        self.structure.remove_atom.assert_any_call(self.atom3)
        self.structure.remove_atom.assert_any_call(self.atom4)


    def test_can_only_remove_chains(self):
        with self.assertRaises(TypeError):
            self.structure.remove_chain("self.chain4")
