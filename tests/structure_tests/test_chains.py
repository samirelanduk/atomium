from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.structures.chains import Chain, ResidueSequence
from atomium.structures.molecules import Residue, Molecule
from atomium.structures.atoms import Atom

class ChainTest(TestCase):

    def setUp(self):
        self.residue1 = Mock(Residue)
        self.residue2 = Mock(Residue)
        self.atom1, self.atom2 = Mock(Atom), Mock(Atom)
        self.atom3, self.atom4 = Mock(Atom), Mock(Atom)
        self.residue1.atoms.return_value = set([self.atom1, self.atom2])
        self.residue2.atoms.return_value = set([self.atom3, self.atom4])
        self.atom1.residue.return_value = self.residue1
        self.atom2.residue.return_value = self.residue1
        self.atom3.residue.return_value = self.residue2
        self.atom4.residue.return_value = self.residue2



class ChainCreationTests(ChainTest):

    @patch("atomium.structures.chains.ResidueSequence.__init__")
    def test_can_create_chain(self, mock_init):
        mock_init.return_value = None
        chain = Chain(self.residue1, self.residue2)
        self.assertIsInstance(chain, ResidueSequence)
        self.assertIsInstance(chain, Molecule)
        self.assertEqual(mock_init.call_count, 1)
        args, kwargs = mock_init.call_args_list[0]
        self.assertEqual(args[0], chain)
        self.assertEqual(set(args[1:]), set(
         [self.residue1, self.residue2]
        ))
        self.assertEqual(chain._chain_id, None)


    def test_can_create_chain_with_id(self):
        chain = Chain(self.residue1, self.residue2, chain_id="A")
        self.assertEqual(chain._chain_id, "A")


    def test_chain_id_must_be_str(self):
        with self.assertRaises(TypeError):
            Chain(self.residue1, self.residue2, chain_id=100)


    def test_chain_id_must_be_unique(self):
        chain = Chain(self.residue1, self.residue2, chain_id="B")
        with self.assertRaises(ValueError):
            Chain(self.residue1, self.residue2, chain_id="B")
