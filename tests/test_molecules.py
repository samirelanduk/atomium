from collections import Counter
from unittest import TestCase
from molecupy import exceptions
from molecupy.molecules import Molecule, AtomicStructure, Atom, Model

class MoleculeTest(TestCase):

    def setUp(self):
        self.atom1 = Atom(1.0, 1.0, 1.0, "H", atom_id=1, atom_name="H1")
        self.atom2 = Atom(1.0, 1.0, 2.0, "C", atom_id=2, atom_name="CA")
        self.atom3 = Atom(1.0, 1.0, 3.0, "O", atom_id=3, atom_name="OX1")
        self.atom1.covalent_bond_to(self.atom2)
        self.atom3.covalent_bond_to(self.atom2)


    def check_valid_molecule(self, molecule):
        self.assertIsInstance(molecule, Molecule)
        self.assertIsInstance(molecule, AtomicStructure)
        for atom in molecule.atoms:
            self.assertEqual(atom.molecule, molecule)
        if molecule.molecule_id is not None:
            self.assertIsInstance(molecule.molecule_id, int)
        if molecule.molecule_name is not None:
            self.assertIsInstance(molecule.molecule_name, str)
        if molecule.model is not None:
            self.assertIsInstance(molecule.model, Model)
        self.assertRegex(str(molecule), r"<Molecule \((\d+) atoms\)>")



class MoleculeCreationTest(MoleculeTest):

    def test_can_create_molecule(self):
        molecule = Molecule(self.atom1, self.atom2, self.atom3)
        self.check_valid_molecule(molecule)


    def test_can_create_molecule_with_id(self):
        molecule = Molecule(self.atom1, self.atom2, self.atom3, molecule_id=10)
        self.check_valid_molecule(molecule)


    def test_molecule_id_must_be_int(self):
        with self.assertRaises(TypeError):
            molecule = Molecule(self.atom1, self.atom2, self.atom3, molecule_id=1.1)
        with self.assertRaises(TypeError):
            molecule = Molecule(self.atom1, self.atom2, self.atom3, molecule_id="10")


    def test_can_create_molecule_with_name(self):
        molecule = Molecule(self.atom1, self.atom2, self.atom3, molecule_name="MOL")
        self.check_valid_molecule(molecule)


    def test_molecule_name_must_be_str(self):
        with self.assertRaises(TypeError):
            molecule = Molecule(self.atom1, self.atom2, self.atom3, molecule_name=1)


    def test_molecule_requires_atoms_to_be_bonded(self):
        bond = list(self.atom3.covalent_bonds)[0]
        self.atom2.covalent_bonds.remove(bond)
        self.atom3.covalent_bonds.remove(bond)
        with self.assertRaises(exceptions.BrokenMoleculeError):
            molecule = Molecule(self.atom1, self.atom2, self.atom3)


    def test_molecule_needs_unique_atom_ids(self):
        atom4 = Atom(1.0, 1.0, 4.0, "H", atom_id=1, atom_name="H1")
        atom5 = Atom(1.0, 1.0, 5.0, "H", atom_name="H1")
        atom6 = Atom(1.0, 1.0, 6.0, "H", atom_name="H1")
        self.atom3.covalent_bond_to(atom4)
        atom4.covalent_bond_to(atom5)
        atom5.covalent_bond_to(atom6)
        with self.assertRaises(exceptions.DuplicateAtomIdError):
            molecule = Molecule(
             self.atom1, self.atom2, self.atom3, atom4, atom5, atom6
            )
        atom4.atom_id = 4
        molecule = Molecule(
         self.atom1, self.atom2, self.atom3, atom4, atom5, atom6
        )


class MoleculeAtomRetrievalTests(MoleculeTest):

    def test_can_get_atom_by_id(self):
        moleclue = Molecule(self.atom1, self.atom2, self.atom3)
        self.assertIs(moleclue.get_atom_by_id(1), self.atom1)
        self.assertIs(moleclue.get_atom_by_id(2), self.atom2)
        self.assertIs(moleclue.get_atom_by_id(3), self.atom3)
        self.assertIs(moleclue.get_atom_by_id(4), None)


    def test_can_only_search_by_numeric_id(self):
        moleclue = Molecule(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(TypeError):
            moleclue.get_atom_by_id(None)



class MoleculeEmpiricalFormulaTests(MoleculeTest):

    def test_can_get_empirical_formula(self):
        atom4 = Atom(1.0, 1.0, 4.0, "O", atom_id=4, atom_name="O")
        atom5 = Atom(1.0, 1.0, 5.0, "N", atom_id=5, atom_name="NA")
        self.atom3.covalent_bond_to(atom4)
        atom4.covalent_bond_to(atom5)
        molecule = Molecule(self.atom1, self.atom2, self.atom3, atom4, atom5)
        self.assertEqual(
         molecule.get_empirical_formula(),
         Counter({"C": 1, "N": 1, "O": 2})
        )
