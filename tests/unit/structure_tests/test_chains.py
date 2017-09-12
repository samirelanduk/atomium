from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.structures.chains import Chain, ResidueSequence
from atomium.structures.molecules import AtomicStructure, Residue, Molecule
from atomium.structures.atoms import Atom

class ChainTest(TestCase):

    def setUp(self):
        self.atom1, self.atom2, self.atom3 = Mock(Atom), Mock(Atom), Mock(Atom)
        self.atom4, self.atom5, self.atom6= Mock(Atom), Mock(Atom), Mock(Atom)
        self.residue1, self.residue2 = Mock(Residue), Mock(Residue)
        self.patcher1 = patch("atomium.structures.chains.ResidueSequence.verify")
        self.mock_verify = self.patcher1.start()
        self.mock_verify.return_value = True
        self.patcher2 = patch("atomium.structures.chains.ResidueSequence.residues")
        self.mock_residues = self.patcher2.start()
        self.mock_residues.return_value = (self.residue1, self.residue2)
        def mock_init(obj, *args, **kwargs):
            obj._atoms = set()
            obj._id_atoms = {}
        self.mock_init = mock_init


    def tearDown(self):
        self.patcher1.stop()
        self.patcher2.stop()



class ChainCreationTests(ChainTest):

    @patch("atomium.structures.chains.Molecule.__init__")
    def test_can_create_chain(self, mock_init):
        mock_init.side_effect = self.mock_init
        chain = Chain(self.atom1, self.atom2, self.atom3)
        self.assertIsInstance(chain, Molecule)
        self.assertIsInstance(chain, ResidueSequence)
        mock_init.assert_called_with(chain, self.atom1, self.atom2, self.atom3)
        self.mock_verify.assert_called_with(chain)


    @patch("atomium.structures.chains.Molecule.__init__")
    def test_can_create_chain_with_name(self, mock_init):
        mock_init.side_effect = self.mock_init
        chain = Chain(self.atom1, self.atom2, name="BORG")
        mock_init.assert_called_with(chain, self.atom1, self.atom2, name="BORG")


    @patch("atomium.structures.chains.Molecule.__init__")
    def test_can_create_chain_with_id(self, mock_init):
        mock_init.side_effect = self.mock_init
        chain = Chain(self.atom1, self.atom2, chain_id="A")
        mock_init.assert_called_with(chain, self.atom1, self.atom2, molecule_id="A")


    def test_atoms_are_linked_to_chain(self):
        chain = Chain(self.atom1, self.atom2, self.atom3)
        self.assertIs(self.atom1._chain, chain)
        self.assertIs(self.atom2._chain, chain)
        self.assertIs(self.atom3._chain, chain)



class ChainReprTests(ChainTest):

    def test_chain_repr_no_id(self):
        chain = Chain(self.atom1, self.atom2, self.atom3)
        chain._name = None
        self.assertEqual(str(chain), "<Chain (2 residues)>")


    def test_chain_repr_no_id(self):
        chain = Chain(self.atom1, self.atom2, self.atom3)
        chain._name = "A"
        self.assertEqual(str(chain), "<Chain A (2 residues)>")



class ChainAtomAdditionTests(ChainTest):

    def test_adding_atom_updates_atom(self):
        chain = Chain(self.atom1, self.atom2)
        chain.add_atom(self.atom3)
        self.assertIs(self.atom3._chain, chain)



class ChainAtomRemovalTests(ChainTest):

    def test_removing_atom_updates_atom(self):
        chain = Chain(self.atom1, self.atom2, self.atom3)
        chain.remove_atom(self.atom3)
        self.assertIs(self.atom3._chain, None)



class ChainIdTests(ChainTest):

    def test_chain_id_property(self):
        chain = Chain(self.atom1, self.atom2, self.atom3, chain_id="A")
        self.assertIs(chain._id, chain.chain_id())
