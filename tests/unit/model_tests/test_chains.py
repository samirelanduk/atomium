from unittest import TestCase
from unittest.mock import patch, Mock, PropertyMock
from atomium.models.molecules import Chain
from atomium.models.structures import AtomStructure
from atomium.models.exceptions import SequenceConnectivityError

class ChainCreationTests(TestCase):

    @patch("atomium.models.structures.AtomStructure.__init__")
    @patch("atomium.models.molecules.Chain.verify")
    def test_can_create_chain(self, mock_verify, mock_init):
        chain = Chain("a", b="c")
        self.assertIsInstance(chain, AtomStructure)
        mock_init.assert_called_with(chain, "a", b="c")
        mock_verify.assert_called_with()
        self.assertEqual(chain._rep_sequence, "")


    @patch("atomium.models.structures.AtomStructure.__init__")
    @patch("atomium.models.molecules.Chain.verify")
    def test_can_create_chain_with_rep_sequence(self, mock_veriy, mock_init):
        chain = Chain("a", b="c", rep="ABC")
        self.assertIsInstance(chain, AtomStructure)
        mock_init.assert_called_with(chain, "a", b="c")
        mock_veriy.assert_called_with()
        self.assertEqual(chain._rep_sequence, "ABC")



class ChainMembershipPropertiesTests(TestCase):

    def test_chains_have_correct_properties(self):
        for attr in ["ligand", "residue"]:
            self.assertIn(attr, Chain.__dict__)
            self.assertIn(attr + "s", Chain.__dict__)
        self.assertIn("model", Chain.__dict__)



class ChainLenTests(TestCase):

    @patch("atomium.models.molecules.Chain.residues")
    def test_can_get_len(self, mock_residues):
        chain = Chain()
        mock_residues.return_value = [1, 2, 4]
        self.assertEqual(len(chain), 3)



class ChainIndexingTests(TestCase):

    @patch("atomium.models.molecules.Chain.residues")
    def test_can_get_len(self, mock_residues):
        chain = Chain()
        mock_residues.return_value = [1, 2, 4]
        self.assertEqual(chain[0], 1)
        self.assertEqual(chain[1], 2)
        self.assertEqual(chain[2], 4)



class ChainCorrectCheckingTests(TestCase):

    def setUp(self):
        self.atom1, self.atom2 = Mock(), Mock()
        self.atom3, self.atom4 = Mock(), Mock()
        self.atom5, self.atom6 = Mock(), Mock()
        self.atom7, self.atom8 = Mock(), Mock()
        self.residue1, self.residue2 = Mock(), Mock()
        self.residue3, self.residue4 = Mock(), Mock()
        self.residue1.next = self.residue2
        self.residue2.next = self.residue3
        self.residue3.next = self.residue4
        self.residue4.next = None
        self.residue1.previous = None
        self.residue2.previous = self.residue1
        self.residue3.previous = self.residue2
        self.residue4.previous = self.residue3
        self.atom1.residue = self.residue1
        self.atom2.residue = self.residue1
        self.atom3.residue = self.residue2
        self.atom4.residue = self.residue2
        self.atom5.residue = self.residue3
        self.atom6.residue = self.residue3
        self.atom7.residue = self.residue4
        self.atom8.residue = self.residue4
        self.chain = Chain()
        self.chain._atoms = set([
         self.atom1, self.atom2, self.atom3, self.atom4,
         self.atom5, self.atom6, self.atom7, self.atom8
        ])


    def test_can_verify_conected_residues(self):
        self.chain._atoms = set([self.atom1, self.atom2])
        self.assertTrue(self.chain.verify())
        self.chain._atoms = set([self.atom3, self.atom4])
        self.assertTrue(self.chain.verify())
        self.chain._atoms = set([self.atom5, self.atom6])
        self.assertTrue(self.chain.verify())
        self.chain._atoms = set([self.atom7, self.atom8])
        self.assertTrue(self.chain.verify())
        self.chain._atoms = set([
         self.atom1, self.atom2, self.atom3, self.atom4
        ])
        self.assertTrue(self.chain.verify())
        self.chain._atoms = set([
         self.atom3, self.atom4, self.atom5, self.atom6
        ])
        self.assertTrue(self.chain.verify())
        self.chain._atoms = set([
         self.atom5, self.atom6, self.atom7, self.atom8
        ])
        self.assertTrue(self.chain.verify())


    def test_can_reject_unconnected_residues(self):
        self.chain._atoms = set([
         self.atom3, self.atom4, self.atom7, self.atom8
        ])
        with self.assertRaises(SequenceConnectivityError):
            self.chain.verify()


    def test_empty_sequences_pass(self):
        self.chain._atoms = set()
        self.assertTrue(self.chain.verify())



class ChainLengthTests(TestCase):

    @patch("atomium.models.molecules.Chain.__len__")
    def test_can_get_len(self, mock_len):
        chain = Chain()
        mock_len.return_value = 100
        self.assertEqual(chain.length, 100)



class ResidueSequenceSequence(TestCase):

    @patch("atomium.models.molecules.Chain.residues")
    def test_can_get_sequence(self, mock_residues):
        chain = Chain()
        mock_residues.return_value = (Mock(code="Y"), Mock(code="T"))
        self.assertEqual(chain.sequence, "YT")



class ChainRepSequenceTests(TestCase):

    def test_chain_rep_sequence_property(self):
        chain = Chain(rep="ABC")
        self.assertIs(chain._rep_sequence, chain.rep_sequence)


    def test_can_update_chain_rep_sequence(self):
        chain = Chain(rep="ABC")
        chain.rep_sequence = "DEF"
        self.assertEqual(chain._rep_sequence, "DEF")



class ChainCopyingTests(TestCase):

    @patch("atomium.models.molecules.Chain.residues")
    @patch("atomium.models.molecules.Chain.ligands")
    def test_can_copy_chain(self, mock_ligands, mock_residues):
        residues = [Mock(), Mock()]
        mock_residues.return_value = residues
        ligands = [Mock(), Mock()]
        mock_ligands.return_value = ligands
        residues[0].copy.return_value, residues[1].copy.return_value = Mock(), Mock()
        ligands[0].copy.return_value, ligands[1].copy.return_value = Mock(), Mock()
        chain = Chain()
        patcher = patch("atomium.models.molecules.Chain")
        mock_chain = patcher.start()
        try:
            new_chain = chain.copy()
            mock_chain.assert_called_with(
             ligands[0].copy.return_value, ligands[1].copy.return_value,
             residues[0].copy.return_value, residues[1].copy.return_value
            )
        finally:
            patcher.stop()


    @patch("atomium.models.molecules.Chain.residues")
    @patch("atomium.models.molecules.Chain.ligands")
    def test_can_copy_chain_with_atoms(self, mock_ligands, mock_residues):
        residues = [Mock(), Mock()]
        mock_residues.return_value = residues
        ligands = [Mock(), Mock()]
        mock_ligands.return_value = ligands
        residues[0].copy.return_value, residues[1].copy.return_value = Mock(), Mock()
        ligands[0].copy.return_value, ligands[1].copy.return_value = Mock(), Mock()
        chain = Chain()
        atoms = [Mock(), Mock(), Mock()]
        chain._atoms = set(atoms)
        patcher = patch("atomium.models.molecules.Chain")
        mock_chain = patcher.start()
        mock_chain.return_value = Mock(_atoms=set(atoms[:-1]))
        try:
            new_chain = chain.copy()
            mock_chain.assert_called_with(
             ligands[0].copy.return_value, ligands[1].copy.return_value,
             residues[0].copy.return_value, residues[1].copy.return_value
            )
            atoms[-1].copy.assert_called_with()
            self.assertFalse(atoms[0].called)
            self.assertFalse(atoms[0].called)
            self.assertEqual(new_chain._atoms, set(atoms[:-1] + [atoms[-1].copy.return_value]))
        finally:
            patcher.stop()
