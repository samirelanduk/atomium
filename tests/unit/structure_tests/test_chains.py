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
        self.patch1 = patch("atomium.structures.chains.ResidueSequence.verify")
        self.mock_verify = self.patch1.start()
        self.mock_verify.return_value = True
        self.patch2 = patch("atomium.structures.chains.ResidueSequence.residues")
        self.mock_residues = self.patch2.start()
        self.mock_residues.return_value = (self.residue1, self.residue2)
        def mock_init(obj, *args, **kwargs):
            obj._atoms = set(args)
        self.patch3 = patch("atomium.structures.molecules.Molecule.__init__")
        self.mock_init = self.patch3.start()
        self.mock_init.side_effect = mock_init


    def tearDown(self):
        self.patch1.stop()
        self.patch2.stop()
        self.patch3.stop()



class ChainCreationTests(ChainTest):

    def test_can_create_chain(self):
        chain = Chain(self.atom1, self.atom2, self.atom3, a=1, b=2)
        self.assertIsInstance(chain, Molecule)
        self.assertIsInstance(chain, ResidueSequence)
        self.mock_init.assert_called_with(chain, self.atom1, self.atom2, self.atom3, a=1, b=2)
        self.mock_verify.assert_called_with(chain)


    def test_atoms_are_linked_to_chain(self):
        self.atom1.id = 100
        chain = Chain(self.atom1, self.atom2, self.atom3)
        self.assertIs(self.atom1._chain, chain)
        self.assertIs(self.atom2._chain, chain)
        self.assertIs(self.atom3._chain, chain)



class ChainReprTests(ChainTest):

    @patch("atomium.structures.chains.Chain.residues")
    def test_chain_repr_no_id(self, mock_res):
        mock_res.return_value = [1, 2]
        chain = Chain(self.atom1, self.atom2, self.atom3)
        chain._name = None
        self.assertEqual(str(chain), "<Chain (2 residues)>")


    @patch("atomium.structures.chains.Chain.residues")
    def test_chain_repr_with_id(self, mock_res):
        mock_res.return_value = [1, 2]
        chain = Chain(self.atom1, self.atom2, self.atom3)
        chain._name = "A"
        self.assertEqual(str(chain), "<Chain A (2 residues)>")
