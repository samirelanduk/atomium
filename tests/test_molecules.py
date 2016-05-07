from unittest import TestCase
from molecupy import exceptions
from molecupy.molecules import Molecule, AtomicStructure, Atom

class MoleculeTest(TestCase):

    def setUp(self):
        self.atom1 = Atom(1.0, 1.0, 1.0, "H", atom_id=1, atom_name="H1")
        self.atom2 = Atom(1.0, 1.0, 2.0, "C", atom_id=2, atom_name="CA")
        self.atom3 = Atom(1.0, 1.0, 3.0, "O", atom_id=3, atom_name="OX1")
        self.atom1.covalent_bond_to(self.atom2)
        self.atom3.covalent_bond_to(self.atom2)


    def check_valid_molecule(self, molecule, check_molecule_id=False, check_molecule_name=False):
        self.assertIsInstance(molecule, Molecule)
        self.assertIsInstance(molecule, Molecule)
        for atom in molecule.atoms:
            self.assertEqual(atom.molecule, molecule)
        if check_molecule_id:
            self.assertIsInstance(molecule.molecule_id, int)
        if check_molecule_name:
            self.assertIsInstance(molecule.molecule_name, str)
        self.assertRegex(str(molecule), r"<Molecule \((\d+) atoms\)>")



class MoleculeCreationTest(MoleculeTest):

    def test_can_create_molecule(self):
        molecule = Molecule(self.atom1, self.atom2, self.atom3)
        self.check_valid_molecule(molecule)


    def test_can_create_molecule_with_id(self):
        molecule = Molecule(self.atom1, self.atom2, self.atom3, molecule_id=10)
        self.check_valid_molecule(molecule, check_molecule_id=True)


    def test_molecule_id_must_be_int(self):
        with self.assertRaises(TypeError):
            molecule = Molecule(self.atom1, self.atom2, self.atom3, molecule_id=1.1)
        with self.assertRaises(TypeError):
            molecule = Molecule(self.atom1, self.atom2, self.atom3, molecule_id="10")


    def test_can_create_molecule_with_name(self):
        molecule = Molecule(self.atom1, self.atom2, self.atom3, molecule_name="MOL")
        self.check_valid_molecule(molecule, check_molecule_name=True)


    def test_molecule_name_must_be_str(self):
        with self.assertRaises(TypeError):
            molecule = Molecule(self.atom1, self.atom2, self.atom3, molecule_name=1)


    def test_molecule_requires_atoms_to_be_bonded(self):
        self.atom2.break_covalent_bond_with(self.atom3)
        with self.assertRaises(exceptions.BrokenMoleculeError):
            molecule = Molecule(self.atom1, self.atom2, self.atom3)
