from unittest import TestCase
from molecupy import exceptions
from molecupy.molecules import Atom, AtomicStructure
from molecupy.macromolecules import Residue, ResiduicStructure

class ResiduicStructureTest(TestCase):

    def setUp(self):
        self.atom1 = Atom(1.0, 1.0, 1.0, "H", atom_id=1, atom_name="H1")
        self.atom2 = Atom(1.0, 1.0, 2.0, "C", atom_id=2, atom_name="CA")
        self.atom3 = Atom(1.0, 1.0, 3.0, "O", atom_id=3, atom_name="OX1")
        self.atom2.covalent_bond_to(self.atom1)
        self.atom2.covalent_bond_to(self.atom3)
        self.residue1 = Residue(1, "MON1", self.atom1, self.atom2, self.atom3)
        self.atom4 = Atom(1.0, 1.0, 4.0, "H", atom_id=4, atom_name="H1")
        self.atom5 = Atom(1.0, 1.0, 5.0, "C", atom_id=5, atom_name="CA")
        self.atom6 = Atom(1.0, 1.0, 6.0, "O", atom_id=6, atom_name="OX1")
        self.atom5.covalent_bond_to(self.atom4)
        self.atom5.covalent_bond_to(self.atom6)
        self.residue2 = Residue(2, "MON2", self.atom4, self.atom5, self.atom6)
        self.atom7 = Atom(1.0, 1.0, 7.0, "H", atom_id=7, atom_name="H1")
        self.atom8 = Atom(1.0, 1.0, 8.0, "C", atom_id=8, atom_name="CA")
        self.atom9 = Atom(1.0, 1.0, 9.0, "O", atom_id=9, atom_name="OX1")
        self.atom8.covalent_bond_to(self.atom7)
        self.atom8.covalent_bond_to(self.atom9)
        self.residue3 = Residue(3, "MON3", self.atom7, self.atom8, self.atom9)


    def check_valid_residuic_structure(self, residuic_structure):
        self.assertIsInstance(residuic_structure, ResiduicStructure)
        self.assertIsInstance(residuic_structure, AtomicStructure)
        self.assertIsInstance(residuic_structure.residues, set)
        for residue in residuic_structure.residues:
            self.assertIsInstance(residue, Residue)
            for atom in residue.atoms:
                self.assertIn(atom, residuic_structure.atoms)
        self.assertRegex(
         str(residuic_structure),
         r"<ResiduicStructure \((\d+) residues\)>"
        )



class ResiduicStructureCreationTests(ResiduicStructureTest):

    def test_can_create_residuic_structure(self):
        residuic_structure = ResiduicStructure(
         self.residue1,
         self.residue2,
         self.residue3
        )
        self.check_valid_residuic_structure(residuic_structure)


    def test_can_residuic_structure_needs_residues(self):
        with self.assertRaises(TypeError):
            residuic_structure = ResiduicStructure(
             self.residue1,
             self.residue2,
             self.atom9
            )


    def test_residuic_structure_needs_at_least_one_residue(self):
        with self.assertRaises(exceptions.NoResiduesError):
            residuic_structure = ResiduicStructure()


    def test_residuic_structure_needs_unique_residue_ids(self):
        self.residue3.residue_id = 2
        with self.assertRaises(exceptions.DuplicateResidueIdError):
            residuic_structure = ResiduicStructure(
             self.residue1, self.residue2, self.residue3
            )
        self.residue3.residue_id = 3
        residuic_structure = ResiduicStructure(
         self.residue1, self.residue2, self.residue3
        )


    def test_residues_in_residuic_structure(self):
        residuic_structure = ResiduicStructure(
         self.residue1,
         self.residue2,
         self.residue3
        )
        self.assertEqual(len(residuic_structure.residues), 3)
        self.assertIn(self.residue1, residuic_structure)
        self.assertIn(self.residue2, residuic_structure)
        self.assertIn(self.residue3, residuic_structure)


    def test_atoms_in_residuic_structure(self):
        residuic_structure = ResiduicStructure(
         self.residue1,
         self.residue2,
         self.residue3
        )
        self.assertEqual(len(residuic_structure.atoms), 9)
        self.assertIn(self.atom1, residuic_structure)
        self.assertIn(self.atom4, residuic_structure)
        self.assertIn(self.atom7, residuic_structure)
