from unittest import TestCase
import unittest.mock
from molecupy.structures import Residue, AtomicStructure, PdbAtom

class ResidueTest(TestCase):

    def setUp(self):
        self.atoms = [unittest.mock.Mock(spec=PdbAtom) for _ in range(10)]



class ResidueCreationTests(ResidueTest):

    def test_can_create_residue(self):
        residue = Residue("A5", "TYR", *self.atoms)
        self.assertIsInstance(residue, AtomicStructure)
        self.assertEqual(residue._residue_id, "A5")
        self.assertEqual(residue._residue_name, "TYR")
        self.assertEqual(residue._atoms, set(self.atoms))


    def test_residue_id_must_be_str(self):
        with self.assertRaises(TypeError):
            Residue(1.1, "TYR", *self.atoms)


    def test_residue_name_must_be_str(self):
        with self.assertRaises(TypeError):
            Residue("A5", 1, *self.atoms)


    def test_residue_repr(self):
        residue = Residue("A5", "TYR", *self.atoms)
        self.assertEqual(str(residue), "<Residue A5 (TYR)>")



class ResiduePropertyTests(ResidueTest):

    def test_residue_properties(self):
        residue = Residue("A5", "TYR", *self.atoms)
        self.assertEqual(residue.residue_id(), "A5")
        self.assertEqual(residue.residue_name(), "TYR")


    def test_can_change_residue_name(self):
        residue = Residue("A5", "TYR", *self.atoms)
        self.assertEqual(residue.residue_name(), "TYR")
        residue.residue_name("ASP")
        self.assertEqual(residue.residue_name(), "ASP")


    def test_residue_name_can_only_be_changed_to_str(self):
        residue = Residue("A5", "TYR", *self.atoms)
        with self.assertRaises(TypeError):
            residue.residue_name(100)
