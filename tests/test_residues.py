from unittest import TestCase
from molecupy import exceptions
from molecupy.structures import PdbResidue, PdbAtom, AtomicStructure

class ResidueTest(TestCase):

    def setUp(self):
        self.atom1 = PdbAtom(1.0, 1.0, 1.0, "H", atom_id=1, atom_name="H1")
        self.atom2 = PdbAtom(1.0, 1.0, 2.0, "C", atom_id=2, atom_name="CA")


    def check_valid_residue(self, residue):
        self.assertIsInstance(residue, PdbResidue)
        self.assertIsInstance(residue, AtomicStructure)
        for atom in residue.atoms:
            self.assertEqual(atom.molecule, residue)
        self.assertIsInstance(residue.residue_id, str)
        self.assertIsInstance(residue.residue_name, str)
        residue.chain
        residue.downstream_residue
        residue.upstream_residue
        self.assertRegex(str(residue), r"<Residue \((.+)\)>")



class ResidueCreationTest(ResidueTest):

    def test_can_create_residue(self):
        residue = PdbResidue("A1", "HET", self.atom1, self.atom2)
        self.check_valid_residue(residue)


    def test_residue_id_must_be_str(self):
        with self.assertRaises(TypeError):
            residue = PdbResidue(1.1, "HET", self.atom1, self.atom2)
        with self.assertRaises(TypeError):
            residue = PdbResidue(1, "HET", self.atom1, self.atom2)


    def test_residue_name_must_be_str(self):
        with self.assertRaises(TypeError):
            residue = PdbResidue(1, 1, self.atom1, self.atom2)



class ResidueConnectionTest(ResidueTest):

    def setUp(self):
        ResidueTest.setUp(self)
        self.atom3 = PdbAtom(1.0, 1.0, 3.0, "H", 3, "H1")
        self.atom4 = PdbAtom(1.0, 1.0, 4.0, "C", 4, "CA")
        self.atom5 = PdbAtom(1.0, 1.0, 5.0, "H", 5, "H1")
        self.atom6 = PdbAtom(1.0, 1.0, 6.0, "C", 6, "CA")
        self.residue1 = PdbResidue("A1", "HET", self.atom1, self.atom2)
        self.residue2 = PdbResidue("A2", "HET", self.atom3, self.atom4)
        self.residue3 = PdbResidue("A3", "HET", self.atom5, self.atom6)


    def test_can_connect_residues(self):
        self.residue1.connect_to(self.residue2, self.atom2, self.atom3)
        self.assertIn(self.atom3, self.atom2.get_covalent_bonded_atoms())
        self.assertIn(self.atom2, self.atom3.get_covalent_bonded_atoms())
        self.residue2.connect_to(self.residue3, self.atom4, self.atom5)
        self.assertIn(self.atom5, self.atom4.get_covalent_bonded_atoms())
        self.assertIn(self.atom4, self.atom5.get_covalent_bonded_atoms())


    def test_can_only_connect_residues_with_valid_atoms(self):
        with self.assertRaises(ValueError):
            self.residue1.connect_to(self.residue2, self.atom2, self.atom5)
        with self.assertRaises(ValueError):
            self.residue1.connect_to(self.residue2, self.atom3, self.atom4)


    def test_connected_residues_know_about_each_other(self):
        self.residue1.connect_to(self.residue2, self.atom2, self.atom3)
        self.residue2.connect_to(self.residue3, self.atom4, self.atom5)
        self.assertIs(self.residue1.upstream_residue, None)
        self.assertIs(self.residue1.downstream_residue, self.residue2)
        self.assertIs(self.residue2.upstream_residue, self.residue1)
        self.assertIs(self.residue2.downstream_residue, self.residue3)
        self.assertIs(self.residue3.upstream_residue, self.residue2)
        self.assertIs(self.residue3.downstream_residue, None)



class ResidueAtomRetrievalTest(ResidueTest):

    def test_can_get_named_alpha_carbon(self):
        residue = PdbResidue("A1", "HET", self.atom1, self.atom2)
        self.assertIs(residue.get_alpha_carbon(), self.atom2)


    def test_can_get_unnamed_alpha_carbon(self):
        self.atom2.atom_name = "C"
        residue = PdbResidue("A1", "HET", self.atom1, self.atom2)
        self.assertIs(residue.get_alpha_carbon(), self.atom2)


    def test_can_get_any_alpha_carbon(self):
        self.atom2.atom_name = "N"
        self.atom2.element = "N"
        residue = PdbResidue("A1", "HET", self.atom1, self.atom2)
        self.assertIn(residue.get_alpha_carbon(), (self.atom1, self.atom2))
