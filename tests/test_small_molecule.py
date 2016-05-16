from unittest import TestCase
from molecupy import exceptions
from molecupy.structures import PdbSmallMolecule, PdbAtom, AtomicStructure

class SmallMoleculeTest(TestCase):

    def setUp(self):
        self.atom1 = PdbAtom(1.0, 1.0, 1.0, "H", atom_id=1, atom_name="H1")
        self.atom2 = PdbAtom(1.0, 1.0, 2.0, "C", atom_id=2, atom_name="CA")


    def check_valid_small_molecule(self, small_molecule):
        self.assertIsInstance(small_molecule, PdbSmallMolecule)
        self.assertIsInstance(small_molecule, AtomicStructure)
        for atom in small_molecule.atoms:
            self.assertEqual(atom.molecule, small_molecule)
        self.assertIsInstance(small_molecule.molecule_id, str)
        self.assertIsInstance(small_molecule.molecule_name, str)
        self.assertRegex(str(small_molecule), r"<SmallMolecule \((.+)\)>")



class SmallMoleculeCreationTest(SmallMoleculeTest):

    def test_can_create_small_molecule(self):
        small_molecule = PdbSmallMolecule("A1", "HET", self.atom1, self.atom2)
        self.check_valid_small_molecule(small_molecule)


    def test_molecule_id_must_be_str(self):
        with self.assertRaises(TypeError):
            small_molecule = PdbSmallMolecule(1.1, "HET", self.atom1, self.atom2)
        with self.assertRaises(TypeError):
            small_molecule = PdbSmallMolecule(1, "HET", self.atom1, self.atom2)


    def test_molecule_name_must_be_str(self):
        with self.assertRaises(TypeError):
            small_molecule = PdbSmallMolecule("A1", 1, self.atom1, self.atom2)
