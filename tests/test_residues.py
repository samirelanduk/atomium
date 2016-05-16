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
        self.assertIsInstance(residue.residue_id, int)
        self.assertIsInstance(residue.residue_name, str)
        residue.chain
        self.assertRegex(str(residue), r"<Residue \((.+)\)>")



class ResidueCreationTest(ResidueTest):

    def test_can_create_residue(self):
        residue = PdbResidue(1, "HET", self.atom1, self.atom2)
        self.check_valid_residue(residue)


    def test_residue_id_must_be_int(self):
        with self.assertRaises(TypeError):
            residue = PdbResidue(1.1, "HET", self.atom1, self.atom2)
        with self.assertRaises(TypeError):
            residue = PdbResidue("1", "HET", self.atom1, self.atom2)


    def test_residue_name_must_be_str(self):
        with self.assertRaises(TypeError):
            residue = PdbResidue(1, 1, self.atom1, self.atom2)
