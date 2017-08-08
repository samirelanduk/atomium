from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.structures.chains import ResidueStructure, ResidueSequence
from atomium.structures.molecules import Residue
from atomium.structures.atoms import Atom

class ResidueSequenceTest(TestCase):

    def setUp(self):
        self.sequence = ResidueSequence()
        self.atom1, self.atom2 = Mock(Atom), Mock(Atom)
        self.atom3, self.atom4 = Mock(Atom), Mock(Atom)
        self.atom5, self.atom6 = Mock(Atom), Mock(Atom)
        self.atom7, self.atom8 = Mock(Atom), Mock(Atom)
        self.residue1, self.residue2 = Mock(Residue), Mock(Residue)
        self.residue3, self.residue4 = Mock(Residue), Mock(Residue)
        self.residue1.residue_id.return_value = "A1"
        self.residue2.residue_id.return_value = "A2"
        self.residue3.residue_id.return_value = "A3"
        self.residue4.residue_id.return_value = "A4"
        self.residue1.name.return_value = "TYR"
        self.residue2.name.return_value = "VAL"
        self.residue3.name.return_value = "TYR"
        self.residue4.name.return_value = "MET"
        self.residue1.next.return_value = self.residue2
        self.residue2.next.return_value = self.residue3
        self.residue3.next.return_value = self.residue4
        self.residue4.next.return_value = None
        self.residue1.previous.return_value = None
        self.residue2.previous.return_value = self.residue1
        self.residue3.previous.return_value = self.residue2
        self.residue4.previous.return_value = self.residue3
        self.atom1.residue.return_value = self.residue1
        self.atom2.residue.return_value = self.residue1
        self.atom3.residue.return_value = self.residue2
        self.atom4.residue.return_value = self.residue2
        self.atom5.residue.return_value = self.residue3
        self.atom6.residue.return_value = self.residue3
        self.atom7.residue.return_value = self.residue4
        self.atom8.residue.return_value = self.residue4
        self.sequence.atoms = lambda: set([
         self.atom1, self.atom2, self.atom3, self.atom4,
         self.atom5, self.atom6, self.atom7, self.atom8
        ])



class ResidueSequenceResiduesTests(ResidueSequenceTest):

    @patch("atomium.structures.chains.ResidueStructure.residues")
    def test_can_get_residues(self, mock_residues):
        mock_residues.return_value = set(
         (self.residue1, self.residue2, self.residue3, self.residue4)
        )
        self.assertEqual(
         self.sequence.residues(),
         (self.residue1, self.residue2, self.residue3, self.residue4)
        )


    @patch("atomium.structures.chains.ResidueStructure.residues")
    def test_can_get_residues_after_filter(self, mock_residues):
        mock_residues.return_value = set(
         (self.residue2, self.residue4)
        )
        self.assertEqual(
         self.sequence.residues(a="a", b="b"),
         (self.residue2, self.residue4)
        )
        mock_residues.assert_called_with(self.sequence, a="a", b="b")
