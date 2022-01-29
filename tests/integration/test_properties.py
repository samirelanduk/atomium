from unittest import TestCase
import atomium

class SpecificStructureProperties(TestCase):

    def test_atom_properties(self):
        pdb = atomium.open("tests/integration/files/1lol.cif")
        atom = pdb.model.atom(1)
        self.assertEqual(atom.mass, 14.0067)
        self.assertEqual(atom.covalent_radius, 0.71)
        self.assertFalse(atom.is_metal)
        self.assertTrue(atom.is_backbone)
        self.assertFalse(atom.is_side_chain)
        atom = pdb.model.atom(5)
        self.assertEqual(atom.mass, 12.0107)
        self.assertEqual(atom.covalent_radius, 0.76)
        self.assertFalse(atom.is_metal)
        self.assertFalse(atom.is_backbone)
        self.assertTrue(atom.is_side_chain)

        pdb = atomium.open("tests/integration/files/12ca.cif")
        atom = pdb.model.atom(2028)
        self.assertEqual(atom.mass, 65.39)
        self.assertEqual(atom.covalent_radius, 1.22)
        self.assertTrue(atom.is_metal)
        self.assertFalse(atom.is_backbone)
        self.assertFalse(atom.is_side_chain)
    

    def test_residue_properties(self):
        pdb = atomium.open("tests/integration/files/1lol.cif")
        residue = pdb.model.residue("A.20")
        self.assertEqual(residue.full_name, "aspartic acid")
        self.assertIs(residue.next, pdb.model.residue("A.21"))
        self.assertIs(residue.previous, pdb.model.residue("A.19"))
        self.assertFalse(residue.in_helix)
        self.assertFalse(residue.in_strand)
        residue = pdb.model.residue("A.11")
        self.assertEqual(residue.full_name, "valine")
        self.assertIsNone(residue.previous)
        self.assertIs(residue.next, pdb.model.residue("A.12"))
        self.assertTrue(residue.in_helix)
        self.assertFalse(residue.in_strand)
        residue = pdb.model.residue("A.222")
        self.assertEqual(residue.full_name, "isoleucine")
        self.assertIsNone(residue.next)
        self.assertIs(residue.previous, pdb.model.residue("A.221"))
        residue = pdb.model.residue("A.15")
        self.assertEqual(residue.full_name, "leucine")
        self.assertTrue(residue.in_strand)


    def test_polymer_properties(self):
        pdb = atomium.open("tests/integration/files/1lol.cif")
        polymer = pdb.model.polymer("A")
        self.assertEqual(polymer.length, 204)
        self.assertTrue(polymer.present_sequence.startswith("VMNRLILAMDLMNRDDALRVTGEVREY"))
        self.assertTrue(polymer.present_sequence.endswith("IVGRSIYLADNPAAAAAGIIESI"))
        self.assertEqual(len(polymer.present_sequence), len(polymer))
        self.assertEqual(len(polymer.helices), 11)
        self.assertEqual(len(polymer.helices[0]), 3)
        self.assertIs(polymer.helices[0][0], pdb.model.residue("A.11"))
        self.assertIs(polymer.helices[0][-1], pdb.model.residue("A.13"))
        self.assertEqual(len(polymer.strands), 9)
        self.assertEqual(len(polymer.strands[0]), 5)
        self.assertIs(polymer.strands[0][0], pdb.model.residue("A.15"))
        self.assertIs(polymer.strands[0][-1], pdb.model.residue("A.19"))


    def test_model_properties(self):
        pdb = atomium.open("tests/integration/files/1lol.cif")
        self.assertEqual(pdb.model.file, pdb)



class AtomStructurePropertyTests(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.pdb = atomium.open("tests/integration/files/1lol.cif")


    def test_atom_structure_mass(self):
        self.assertEqual(round(self.pdb.model.mass), 46018)
        self.assertEqual(round(self.pdb.model.residue("A.11").mass, 1), 90.1)
        self.assertEqual(round(self.pdb.model.polymer("A").mass, 1), 20630.9)
        self.assertEqual(round(self.pdb.model.non_polymer("E").mass, 1), 80)
    

    def test_atom_structure_charge(self):
        self.pdb.model.atom(1).charge = -1
        self.pdb.model.atom(3222).charge = 0.5
        self.assertEqual(self.pdb.model.charge, -0.5)
        self.assertEqual(self.pdb.model.residue("A.11").charge, -1)
        self.assertEqual(self.pdb.model.polymer("A").charge, -1)
        self.assertEqual(self.pdb.model.non_polymer("E").charge, 0.5)