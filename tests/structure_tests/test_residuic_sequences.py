from unittest import TestCase
import unittest.mock
from molecupy.structures import ResiduicSequence, ResiduicStructure, Residue, Atom

class ResiduicSequenceTest(TestCase):

    def setUp(self):
        self.residues = [unittest.mock.Mock(spec=Residue) for _ in range(10)]



class ResiduicSequenceCreationTests(ResiduicSequenceTest):

    def test_can_create_residuic_structure(self):
        residuic_sequence = ResiduicSequence(*self.residues)
        self.assertIsInstance(residuic_sequence, ResiduicStructure)
        self.assertEqual(residuic_sequence._residues, self.residues)


    def test_can_only_create_residuic_sequence_with_residues(self):
        with self.assertRaises(TypeError):
            ResiduicSequence("Atom1", "Atom2")


    def test_residuic_structure_repr(self):
        residuic_sequence = ResiduicSequence(*self.residues)
        self.assertEqual(str(residuic_sequence), "<ResiduicSequence (10 residues)>")


    def test_residuic_sequence_length(self):
        residuic_sequence = ResiduicSequence(*self.residues)
        self.assertEqual(len(residuic_sequence), 10)



class ResiduicSequencePropertyTests(ResiduicSequenceTest):

    def test_can_get_residues(self):
        residuic_sequence = ResiduicSequence(*self.residues)
        self.assertEqual(residuic_sequence.residues(), self.residues)


    def test_residuic_structure_residues_is_read_only(self):
        residue11 = unittest.mock.Mock(spec=Residue)
        residuic_sequence = ResiduicSequence(*self.residues)
        self.assertEqual(len(residuic_sequence.residues()), 10)
        residuic_sequence.residues().append(residue11)
        self.assertEqual(len(residuic_sequence.residues()), 10)


    def test_can_exclude_missing_residues(self):
        for residue in self.residues:
            residue.is_missing.return_value = False
        self.residues[6].is_missing.return_value = True
        residuic_sequence = ResiduicSequence(*self.residues)
        self.assertEqual(
         len(residuic_sequence.residues(include_missing=False)),
         9
        )
        self.assertNotIn(
         self.residues[6],
         residuic_sequence.residues(include_missing=False)
        )


    def test_can_add_residue(self):
        residue11 = unittest.mock.Mock(spec=Residue)
        residuic_sequence = ResiduicSequence(*self.residues)
        residuic_sequence.add_residue(residue11)
        self.assertEqual(len(residuic_sequence.residues()), 11)
        self.assertEqual(residue11, residuic_sequence.residues()[-1])


    def test_can_only_add_residues(self):
        residuic_sequence = ResiduicSequence(*self.residues)
        with self.assertRaises(TypeError):
            residuic_sequence.add_residue("atom21")


    def test_can_remove_residues(self):
        residuic_sequence = ResiduicSequence(*self.residues)
        residuic_sequence.remove_residue(self.residues[5])
        self.assertEqual(len(residuic_sequence.residues()), 9)
        self.assertNotIn(self.residues[5], residuic_sequence.residues())


    def test_can_get_atoms(self):
        atoms = [unittest.mock.Mock(spec=Atom) for _ in range(10)]
        for index, residue in enumerate(self.residues):
            residue.atoms.return_value = set([atoms[index]])
        residuic_sequence = ResiduicSequence(*self.residues)
        self.assertEqual(residuic_sequence._atoms, set(atoms))
        self.assertEqual(residuic_sequence.atoms(), set(atoms))



class SequenceGenerationTests(ResiduicSequenceTest):

    def setUp(self):
        ResiduicSequenceTest.setUp(self)
        names = ["VAL", "TYR", "TRP", "LEU", "ILE"]
        for index, residue in enumerate(self.residues[:5]):
            residue.residue_name.return_value = names[index]
            residue.is_missing.return_value = False if index else True
        self.residuic_sequence = ResiduicSequence(*self.residues[:5])


    def test_can_get_string_sequence(self):
        self.assertEqual(
         self.residuic_sequence.sequence_string(),
         "VYWLI"
        )


    def test_can_get_protein_sequence_with_unknown_residues(self):
        self.residues[2].residue_name.return_value = "ABC"
        self.assertEqual(
         self.residuic_sequence.sequence_string(),
         "VYXLI"
        )


    def test_can_ignore_missing_residues(self):
        self.assertEqual(
         self.residuic_sequence.sequence_string(include_missing=False),
         "YWLI"
        )
