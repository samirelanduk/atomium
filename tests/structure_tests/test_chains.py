from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.structures.chains import Chain
from atomium.structures.molecules import AtomicStructure, Residue, Molecule
from atomium.structures.atoms import Atom

class ChainTest(TestCase):

    def setUp(self):
        self.atom1, self.atom2, self.atom3 = Mock(Atom), Mock(Atom), Mock(Atom)
        self.atom4, self.atom5, self.atom6= Mock(Atom), Mock(Atom), Mock(Atom)
        self.residue1, self.residue2 = Mock(Residue), Mock(Residue)



class ChainCreationTests(ChainTest):

    @patch("atomium.structures.chains.Molecule.__init__")
    def test_can_create_chain(self, mock_init):
        mock_init.return_value = None
        chain = Chain(self.atom1, self.atom2, self.atom3)
        self.assertIsInstance(chain, Molecule)
        mock_init.assert_called_with(chain, self.atom1, self.atom2, self.atom3)


    @patch("atomium.structures.chains.Molecule.__init__")
    def test_can_create_chain_with_name(self, mock_init):
        mock_init.return_value = None
        chain = Chain(self.atom1, self.atom2, name="BORG")
        mock_init.assert_called_with(chain, self.atom1, self.atom2, name="BORG")


    @patch("atomium.structures.chains.Molecule.__init__")
    def test_can_create_chain_with_id(self, mock_init):
        mock_init.return_value = None
        chain = Chain(self.atom1, self.atom2, chain_id="A")
        mock_init.assert_called_with(chain, self.atom1, self.atom2, molecule_id="A")


    def test_atoms_are_linked_to_chain(self):
        chain = Chain(self.atom1, self.atom2, self.atom3)
        self.assertIs(self.atom1._chain, chain)
        self.assertIs(self.atom2._chain, chain)
        self.assertIs(self.atom3._chain, chain)
