from unittest import TestCase
import unittest.mock
from molecupy.structures import Residue, AtomicStructure, Atom, GhostAtom
from molecupy.exceptions import MultipleResidueConnectionError

class ResidueTest(TestCase):

    def setUp(self):
        self.atoms = [unittest.mock.Mock(Atom) for _ in range(10)]
        for atom in self.atoms:
            atom._molecule = None



class ResidueCreationTests(ResidueTest):

    def test_can_create_residue(self):
        residue = Residue("A5", "TYR", *self.atoms)
        self.assertIsInstance(residue, AtomicStructure)
        self.assertEqual(residue._residue_id, "A5")
        self.assertEqual(residue._residue_name, "TYR")
        self.assertEqual(residue._chain, None)
        self.assertEqual(residue._atoms, set(self.atoms))


    def test_residue_id_must_be_str(self):
        with self.assertRaises(TypeError):
            Residue(1.1, "TYR", *self.atoms)


    def test_molecule_id_must_be_valid(self):
        with self.assertRaises(ValueError):
            Residue("100", "TYR", *self.atoms)
        with self.assertRaises(ValueError):
            Residue("B_100", "TYR", *self.atoms)
        with self.assertRaises(ValueError):
            Residue("CA34", "TYR", *self.atoms)


    def test_residue_name_must_be_str(self):
        with self.assertRaises(TypeError):
            Residue("A5", 1, *self.atoms)


    def test_residue_repr(self):
        residue = Residue("A5", "TYR", *self.atoms)
        self.assertEqual(str(residue), "<Residue A5 (TYR)>")


    def test_residues_update_atoms(self):
        for atom in self.atoms:
            self.assertIs(atom._molecule, None)
        residue = Residue("A5", "TRP", *self.atoms)
        for atom in self.atoms:
            self.assertIs(atom._molecule, residue)



class ResiduePropertyTests(ResidueTest):

    def test_residue_properties(self):
        residue = Residue("A5", "TYR", *self.atoms)
        self.assertEqual(residue.residue_id(), "A5")
        self.assertEqual(residue.residue_name(), "TYR")
        self.assertEqual(residue.chain(), None)


    def test_can_change_residue_name(self):
        residue = Residue("A5", "TYR", *self.atoms)
        self.assertEqual(residue.residue_name(), "TYR")
        residue.residue_name("ASP")
        self.assertEqual(residue.residue_name(), "ASP")


    def test_residue_name_can_only_be_changed_to_str(self):
        residue = Residue("A5", "TYR", *self.atoms)
        with self.assertRaises(TypeError):
            residue.residue_name(100)


    def test_residue_integrates_added_atoms(self):
        residue = Residue("A5", "MET", *self.atoms[:-1])
        self.assertIs(self.atoms[-1]._molecule, None)
        residue.add_atom(self.atoms[-1])
        self.assertIs(self.atoms[-1]._molecule, residue)
        self.assertEqual(len(residue.atoms()), 10)


    def test_residue_deintegrates_removed_atoms(self):
        residue = Residue("A5", "MET", *self.atoms)
        self.assertIs(self.atoms[-1]._molecule, residue)
        residue.remove_atom(self.atoms[-1])
        self.assertIs(self.atoms[-1]._molecule, None)
        self.assertEqual(len(residue.atoms()), 9)



class MissingResidueTests(ResidueTest):

    def test_residue_knows_if_it_is_missing(self):
        generic_atoms = [unittest.mock.Mock(spec=GhostAtom) for _ in range(10)]
        whole_residue = Residue("A3", "MET", *self.atoms)
        part_residue = Residue("A2", "MET", self.atoms[0], self.atoms[1], generic_atoms[0])
        missing_residue = Residue("A1", "MET", *generic_atoms)
        self.assertFalse(whole_residue.is_missing())
        self.assertFalse(part_residue.is_missing())
        self.assertTrue(missing_residue.is_missing())



class ResidueConnectionTests(ResidueTest):

    def setUp(self):
        ResidueTest.setUp(self)
        self.residue1 = Residue("A1", "MET", *self.atoms[:2])
        self.residue2 = Residue("A2", "MET", *self.atoms[2:4])
        self.residue3 = Residue("A3", "MET", *self.atoms[4:6])
        self.residue4 = Residue("A4", "MET", *self.atoms[6:8])
        self.residue5 = Residue("A5", "MET", *self.atoms[8:])


    def test_can_connect_residues(self):
        self.assertIs(self.residue1.downstream_residue(), None)
        self.assertIs(self.residue1.upstream_residue(), None)
        self.assertIs(self.residue2.downstream_residue(), None)
        self.assertIs(self.residue2.upstream_residue(), None)
        self.assertIs(self.residue3.downstream_residue(), None)
        self.assertIs(self.residue3.upstream_residue(), None)
        self.assertIs(self.residue4.downstream_residue(), None)
        self.assertIs(self.residue4.upstream_residue(), None)
        self.assertIs(self.residue5.downstream_residue(), None)
        self.assertIs(self.residue5.upstream_residue(), None)
        self.residue1.connect_to(self.residue2)
        self.residue2.connect_to(self.residue3)
        self.residue3.connect_to(self.residue4)
        self.residue4.connect_to(self.residue5)
        self.assertIs(self.residue1.downstream_residue(), self.residue2)
        self.assertIs(self.residue1.upstream_residue(), None)
        self.assertIs(self.residue2.downstream_residue(), self.residue3)
        self.assertIs(self.residue2.upstream_residue(), self.residue1)
        self.assertIs(self.residue3.downstream_residue(), self.residue4)
        self.assertIs(self.residue3.upstream_residue(), self.residue2)
        self.assertIs(self.residue4.downstream_residue(), self.residue5)
        self.assertIs(self.residue4.upstream_residue(), self.residue3)
        self.assertIs(self.residue5.downstream_residue(), None)
        self.assertIs(self.residue5.upstream_residue(), self.residue4)


    def test_can_only_connect_residues(self):
        with self.assertRaises(TypeError):
            self.residue1.connect_to(self.atoms[0])


    def test_cannot_connect_residue_to_connected_residue(self):
        self.residue2.connect_to(self.residue3)
        with self.assertRaises(MultipleResidueConnectionError):
            self.residue2.connect_to(self.residue4)
        with self.assertRaises(MultipleResidueConnectionError):
            self.residue1.connect_to(self.residue3)


    def test_can_disconnect_residues(self):
        self.residue1.connect_to(self.residue2)
        self.residue2.connect_to(self.residue3)
        self.residue3.connect_to(self.residue4)
        self.residue4.connect_to(self.residue5)
        self.residue2.disconnect_from(self.residue3)
        self.assertIs(self.residue2.downstream_residue(), None)
        self.assertIs(self.residue3.upstream_residue(), None)
        self.residue4.disconnect_from(self.residue3)
        self.assertIs(self.residue3.downstream_residue(), None)
        self.assertIs(self.residue4.upstream_residue(), None)



class ResidueAtomRetrievalTest(ResidueTest):

    def setUp(self):
        ResidueTest.setUp(self)
        self.atom1 = self.atoms[0]
        self.atom1.atom_name.return_value = "CA"
        self.atom1.element.return_value = "C"
        self.atom2 = self.atoms[1]
        self.atom1.atom_name.return_value = "HX"
        self.atom2.element.return_value = "H"
        self.atom3 = self.atoms[2]
        self.atom1.atom_name.return_value = "N"
        self.atom3.element.return_value = "N"
        self.residue = Residue("A1", "RES", self.atom1, self.atom2, self.atom3)


    def test_can_get_named_alpha_carbon(self):
        self.assertIs(self.residue.alpha_carbon(), self.atom1)


    def test_can_get_unnamed_alpha_carbon(self):
        self.atom1.atom_name.return_value = "C"
        self.assertIs(self.residue.alpha_carbon(), self.atom1)


    def test_can_get_any_alpha_carbon(self):
        self.atom1.atom_name.return_value = "C"
        self.atom1.element.return_value = "P"
        self.assertIn(
         self.residue.alpha_carbon(),
         (self.atom1, self.atom2, self.atom3)
        )
