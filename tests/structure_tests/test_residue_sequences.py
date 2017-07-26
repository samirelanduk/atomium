from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.structures.chains import ResidueSequence
from atomium.structures.molecules import Residue, AtomicStructure
from atomium.structures.atoms import Atom

class ResidueSequenceTest(TestCase):

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



class ResidueSequenceCreationTests(ResidueSequenceTest):

    @patch("atomium.structures.molecules.AtomicStructure.__init__")
    def test_can_create_residue_sequence(self, mock_init):
        mock_init.return_value = None
        sequence = ResidueSequence(self.residue1, self.residue2)
        self.assertIsInstance(sequence, AtomicStructure)
        self.assertEqual(mock_init.call_count, 1)
        args, kwargs = mock_init.call_args_list[0]
        self.assertEqual(args[0], sequence)
        self.assertEqual(set(args[1:]), set(
         [self.atom1, self.atom2, self.atom3, self.atom4]
        ))
        self.assertEqual(sequence._residues, [self.residue1, self.residue2])


    def test_residue_sequence_needs_residues(self):
        with self.assertRaises(TypeError):
            ResidueSequence(self.residue1, "self.residue2")


    def test_residue_needs_at_least_one_residue(self):
        with self.assertRaises(ValueError):
            ResidueSequence()



class ResidueSequenceReprTests(ResidueSequenceTest):

    def test_residue_sequence_repr(self):
        sequence = ResidueSequence(self.residue1, self.residue2)
        self.assertEqual(str(sequence), "<ResidueSequence (2 residues)>")



class ResidueSequenceLenTests(ResidueSequenceTest):

    def test_residue_sequence_len(self):
        sequence = ResidueSequence(self.residue1, self.residue2)
        self.assertEqual(len(sequence), 2)



class ResidueSequenceLengthTests(ResidueSequenceTest):

    @patch("atomium.structures.chains.ResidueSequence.__len__")
    def test_residue_sequence_lenght_is_len(self, mock_len):
        mock_len.return_value = 100
        sequence = ResidueSequence(self.residue1, self.residue2)
        self.assertEqual(sequence.length(), 100)



class ResidueSequenceIterableTests(ResidueSequenceTest):

    def test_residue_sequence_is_iterable(self):
        sequence = ResidueSequence(self.residue1, self.residue2)
        for residue, correct_residue in zip(sequence, (self.residue1, self.residue2)):
            self.assertEqual(residue, correct_residue)
