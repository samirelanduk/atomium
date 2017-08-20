from unittest import TestCase
from unittest.mock import Mock, patch, MagicMock
from atomium.structures.chains import ResidueStructure
from atomium.structures.molecules import Residue
from atomium.structures.atoms import Atom

class ResidueStructureTest(TestCase):

    def setUp(self):
        self.structure = ResidueStructure()
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
        self.atom1.residue.return_value = self.residue1
        self.atom2.residue.return_value = self.residue1
        self.atom3.residue.return_value = self.residue2
        self.atom4.residue.return_value = self.residue2
        self.atom5.residue.return_value = self.residue3
        self.atom6.residue.return_value = self.residue3
        self.atom7.residue.return_value = self.residue4
        self.atom8.residue.return_value = self.residue4
        self.residue1.atoms.return_value = set([self.atom1, self.atom2])
        self.residue2.atoms.return_value = set([self.atom3, self.atom4])
        self.residue3.atoms.return_value = set([self.atom5, self.atom6])
        self.residue4.atoms.return_value = set([self.atom7, self.atom8])
        self.structure.atoms = lambda: set([
         self.atom1, self.atom2, self.atom3, self.atom4,
         self.atom5, self.atom6, self.atom7, self.atom8
        ])
        self.structure.add_atom = MagicMock()
        self.structure.remove_atom = MagicMock()



class ResidueStructureResiduesTests(ResidueStructureTest):

    def test_can_get_residues(self):
        self.assertEqual(
         self.structure.residues(),
         set([self.residue1, self.residue2, self.residue3, self.residue4])
        )


    def test_can_filter_none_from_residues(self):
        self.atom4.residue.return_value = None
        self.assertEqual(
         self.structure.residues(),
         set([self.residue1, self.residue2, self.residue3, self.residue4])
        )


    def test_can_get_residues_by_id(self):
        self.assertEqual(
         self.structure.residues(residue_id="A1"), set([self.residue1])
        )
        self.assertEqual(
         self.structure.residues(residue_id="A2"), set([self.residue2])
        )
        self.assertEqual(self.structure.residues(residue_id="A5"), set())


    def test_can_get_residues_by_name(self):
        self.assertEqual(
         self.structure.residues(name="VAL"), set([self.residue2])
        )
        self.assertEqual(
         self.structure.residues(name="TYR"), set([self.residue1, self.residue3])
        )
        self.assertEqual(self.structure.residues(name="GLY"), set())



class ResidueStructureResidueTests(ResidueStructureTest):

    @patch("atomium.structures.chains.ResidueStructure.residues")
    def test_residue_calls_residues(self, mock_residues):
        mock_residues.return_value = set([self.residue4])
        residue = self.structure.residue(name="A")
        mock_residues.assert_called_with(name="A")
        self.assertIs(residue, self.residue4)


    @patch("atomium.structures.chains.ResidueStructure.residues")
    def test_residue_can_return_none(self, mock_residues):
        mock_residues.return_value = set()
        self.assertIs(self.structure.residue(name="C"), None)


    @patch("atomium.structures.chains.ResidueStructure.residues")
    def test_residue_can_get_residue_by_id_and_name(self, mock_residues):
        mock_residues.return_value = set([self.residue1])
        residue = self.structure.residue(residue_id="A1", name="A")
        mock_residues.assert_called_with(residue_id="A1", name="A")
        self.assertIs(residue, self.residue1)



class ResidueStructureResidueAdditionTests(ResidueStructureTest):

    def test_can_add_residue(self):
        self.structure._atoms = set([
         self.atom1, self.atom2, self.atom3, self.atom4,
         self.atom5, self.atom6
        ])
        self.structure.add_residue(self.residue4)
        self.structure.add_atom.assert_any_call(self.atom7)
        self.structure.add_atom.assert_any_call(self.atom8)


    def test_can_only_add_residues(self):
        with self.assertRaises(TypeError):
            self.structure.add_residue("self.residue4")



class ResidueStructureResidueRemovalTests(ResidueStructureTest):

    def test_can_remove_residue(self):
        self.structure.remove_residue(self.residue4)
        self.structure.remove_atom.assert_any_call(self.atom7)
        self.structure.remove_atom.assert_any_call(self.atom8)


    def test_can_only_remove_residues(self):
        with self.assertRaises(TypeError):
            self.structure.remove_residue("self.residue4")
