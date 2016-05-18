from unittest import TestCase
from molecupy import exceptions
from molecupy.structures import ResiduicStructure, AtomicStructure, PdbAtom, PdbResidue

class ResiduicStructureTest(TestCase):

    def setUp(self):
        self.atom1 = PdbAtom(1.0, 1.0, 1.0, "H", 1, "H1")
        self.atom2 = PdbAtom(1.0, 1.0, 2.0, "C", 2, "CA")
        self.atom3 = PdbAtom(1.0, 1.0, 3.0, "O", 3, "OX1")
        self.residue1 = PdbResidue("A1", "ARG", self.atom1, self.atom2, self.atom3)
        self.atom4 = PdbAtom(1.0, 1.0, 4.0, "H", 4, "H1")
        self.atom5 = PdbAtom(1.0, 1.0, 5.0, "C", 5, "CA")
        self.atom6 = PdbAtom(1.0, 1.0, 6.0, "O", 6, "OX1")
        self.residue2 = PdbResidue("A2", "HST", self.atom4, self.atom5, self.atom6)
        self.atom7 = PdbAtom(1.0, 1.0, 7.0, "H", 7, "H1")
        self.atom8 = PdbAtom(1.0, 1.0, 8.0, "C", 8, "CA")
        self.atom9 = PdbAtom(1.0, 1.0, 9.0, "O", 9, "OX1")
        self.residue3 = PdbResidue("A3", "TRP", self.atom7, self.atom8, self.atom9)


    def check_valid_residuic_structure(self, residuic_structure):
        self.assertIsInstance(residuic_structure, ResiduicStructure)
        self.assertIsInstance(residuic_structure, AtomicStructure)
        self.assertIsInstance(residuic_structure.residues, set)
        for residue in residuic_structure.residues:
            self.assertIsInstance(residue, PdbResidue)
        self.assertRegex(str(residuic_structure), r"<ResiduicStructure \((\d+) residues\)>")



class ResiduicStructureCreationTests(ResiduicStructureTest):

    def test_can_create_residuic_structure(self):
        residuic_structure = ResiduicStructure(self.residue1, self.residue2, self.residue3)
        self.check_valid_residuic_structure(residuic_structure)


    def test_residuic_structure_needs_residues(self):
        with self.assertRaises(TypeError):
            residuic_structure = ResiduicStructure(self.residue1, self.residue2, "res3")


    def test_residuic_structure_needs_at_least_one_residue(self):
        with self.assertRaises(exceptions.NoResiduesError):
            residuic_structure = ResiduicStructure()


    def test_atoms_in_residuic_structure(self):
        residuic_structure = ResiduicStructure(self.residue1, self.residue2, self.residue3)
        self.assertEqual(len(residuic_structure.atoms), 9)
        self.assertIn(self.atom1, residuic_structure)
        self.assertIn(self.atom5, residuic_structure)
        self.assertIn(self.atom9, residuic_structure)


    def test_residues_in_residuic_structure(self):
        residuic_structure = ResiduicStructure(self.residue1, self.residue2, self.residue3)
        self.assertEqual(len(residuic_structure.residues), 3)
        self.assertIn(self.residue1, residuic_structure)
        self.assertIn(self.residue2, residuic_structure)
        self.assertIn(self.residue3, residuic_structure)


    def test_residuic_structure_needs_unique_residue_ids(self):
        self.residue3.residue_id = "A2"
        with self.assertRaises(exceptions.DuplicateResidueIdError):
            residuic_structure = ResiduicStructure(
             self.residue1, self.residue2, self.residue3
            )
        self.residue3.residue_id = "A3"
        residuic_structure = ResiduicStructure(self.residue1, self.residue2, self.residue3)



class ResiduicStructureResidueRetrievalTests(ResiduicStructureTest):

    def test_can_get_residue_by_name(self):
        residuic_structure = ResiduicStructure(self.residue1, self.residue2, self.residue3)
        self.assertIs(residuic_structure.get_residue_by_name("ARG"), self.residue1)
        self.assertIs(residuic_structure.get_residue_by_name("xxx"), None)


    def test_can_get_multiple_residues_by_name(self):
        residuic_structure = ResiduicStructure(self.residue1, self.residue2, self.residue3)
        self.assertEqual(
         residuic_structure.get_residues_by_name("HST"),
         set([self.residue2])
        )
        self.assertEqual(residuic_structure.get_residues_by_name("XXX"), set())


    def test_can_only_search_by_string_name(self):
        residuic_structure = ResiduicStructure(self.residue1, self.residue2, self.residue3)
        with self.assertRaises(TypeError):
            residuic_structure.get_residue_by_name(None)
        with self.assertRaises(TypeError):
            residuic_structure.get_residues_by_name(None)


    def test_can_get_residue_by_id(self):
        residuic_structure = ResiduicStructure(self.residue1, self.residue2, self.residue3)
        self.assertIs(residuic_structure.get_residue_by_id("A1"), self.residue1)
        self.assertIs(residuic_structure.get_residue_by_id("A2"), self.residue2)
        self.assertIs(residuic_structure.get_residue_by_id("A3"), self.residue3)
        self.assertIs(residuic_structure.get_residue_by_id("A4"), None)


    def test_can_only_search_by_str_id(self):
        residuic_structure = ResiduicStructure(self.residue1, self.residue2, self.residue3)
        with self.assertRaises(TypeError):
            residuic_structure.get_residue_by_id(None)
